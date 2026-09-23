import numpy as np
import lmfit
import logging

from .base import Fitter, FitResult
from ..models.base import ParameterSet
from ..components import NoDataError
from ..timeseries.misc import robust_std


class LmfitWrapper(Fitter):
    """A wrapper for the LM-fitting of the chunks
    
    :param model: A :class:`SimpleModel`, which contains all the required 
        submodels that are used in the fitting procedure.
    :type model: :class:`SimpleModel`
    """
    
    def __init__(self, model):
        self.model = model

    @staticmethod
    def convert_params(params, to_lmfit=False, from_lmfit=False):
        """
        Convert between :class:`ParameterSet` and :class:`lmfit.Parameters`
        
        :param params: An object containing fitting parameters.
        :type params: :class:`ParameterSet` or :class:`lmfit.Parameters`
        :param to_lmfit: Convert input to :class:`lmfit.Parameters`.
        :type to_lmfit: bool
        :param from_lmfit: Convert input to :class:`ParameterSet`.
        :type from_lmfit: bool
        
        :return: The converted parameter object.
        :rtype: :class:`ParameterSet` or :class:`lmfit.Parameters`
        """
        if to_lmfit:
            lmfit_params = lmfit.Parameters()
            lmfit_params.add_many(*params.items())
            return lmfit_params
        elif from_lmfit:
            pyodine_params = ParameterSet()
            for k in params:
                pyodine_params[k] = params[k].value
            return pyodine_params
        else:
            raise ValueError('Either `to_lmfit` or `from_lmfit` must be true!')
    
    
    def fit_ostar(self, chunk, weight=None, chunk_ind=None):
        """Convenience function for fitting with fixed velocity and no template
        
        :param chunk: The chunk to be modelled.
        :type chunk: :class:`Chunk`
        :param weight: Pixel weights to use in the model evaluation. Defaults
            to None.
        :type weight: ndarray[nr_pix], or None
        :param chunk_ind: Index of the chunk, to grab the respective chunk 
            from the template. Defaults to None.
        :type chunk_ind: int, or None
        
        :return: The best-fit result.
        :rtype: :class:`LmfitResult`
        """
        params = self.model.guess_params(chunk)
        lmfit_params = self.convert_params(params, to_lmfit=True)
        # Fix parameters
        lmfit_params['velocity'].vary = False
        lmfit_params['tem_depth'].vary = False
        return self.fit(chunk, lmfit_params, weight=weight, chunk_ind=chunk_ind)
    

    def fit(self, chunk, lmfit_params, weight=None, chunk_ind=None, **kwargs):
        """Fit the chunk and return the best-fit result
        
        :param chunk: The chunk to be modelled.
        :type chunk: :class:`Chunk`
        :param lmfit_params: The parameter object, defining starting values, 
            limits, etc.
        :type lmfit_params: :class:`lmfit.Parameters`
        :param weight: Pixel weights to use in the model evaluation. Defaults
            to None.
        :type weight: ndarray[nr_pix], or None
        :param chunk_ind: Index of the chunk, to grab the respective chunk 
            from the template. Defaults to None.
        :type chunk_ind: int, or None
        
        :return: The best-fit result.
        :rtype: :class:`LmfitResult`
        """
        
        # Add-on: pixel weights in the fitting function, as used in dop code
        def func(lmfit_params, x, weight, chunk_ind):
            params = self.convert_params(lmfit_params, from_lmfit=True)
            if isinstance(weight, (list, np.ndarray)):
                return (self.model.eval(chunk, params, require=None, chunk_ind=chunk_ind) - chunk.flux) * np.sqrt(np.abs(weight))
            else:
                return self.model.eval(chunk, params, require=None, chunk_ind=chunk_ind) - chunk.flux

        # Carry out the fit
        try:
            # Make sure that the initial parameter guesses are consistent with
            # template and iodine atlas coverage
            params = self.convert_params(lmfit_params, from_lmfit=True)
            self.model.eval(chunk, params, require='full', chunk_ind=chunk_ind)
            # Carry out the fit
            lmfit_result = lmfit.minimize(func, lmfit_params, args=[chunk.pix, weight, chunk_ind]) #, xtol=1.e-7)
            # Make sure that the fitted parameters are consistent with
            # template and iodine atlas coverage
            new_params = self.convert_params(lmfit_result.params, from_lmfit=True)
            self.model.eval(chunk, new_params, require='full', chunk_ind=chunk_ind)
            # Return output as LmfitResult object
            return self.LmfitResult(chunk, self.model, lmfit_result, chunk_ind=chunk_ind)
        except NoDataError:
            logging.error('Chunk {}: No Data! Returning LmfitResult with None.'.format(chunk_ind))
            return self.LmfitResult(chunk, self.model, None, chunk_ind=chunk_ind)


    class LmfitResult(FitResult):
        """Results from LmfitWrapper
        
        :param chunk: The chunk which was modelled.
        :type chunk: :class:`Chunk`
        :param model: The model instance used in the fitting procedure.
        :type model: :class:`SimpleModel`
        :param lmfit_result: The best-fit results from the modelling. None if 
            it failed.
        :type lmfit_result: :class:`lmfit.MinimizerResult`
        :param chunk_ind: The index of the modelled chunk. Defaults to None.
        :type chunk_ind: int, or None
        """
        def __init__(self, chunk, model, lmfit_result, chunk_ind=None):
            self.chunk = chunk
            self.model = model
            self.lmfit_result = lmfit_result
            self.chunk_ind = chunk_ind
        
        def rel_residuals_rms(self, robust=True):
            """Compute the mean relative residuals of the chunk
            
            :param robust: If True, return the robust mean (default). Otherwise
                the normal mean.
            :type robust: bool
            
            :return: The (robust) mean of the relative residuals between data
                and model.
            :rtype: float
            """
            if self.lmfit_result is not None:
                if robust:
                    return robust_std(self.residuals/self.fitted_spectrum.flux)
                else:
                    return np.nanstd(self.residuals/self.fitted_spectrum.flux)
            else:
                return np.NaN

        @property
        def params(self):
            """Return a :class:'ParameterSet' with the fitted parameters"""
            if self.lmfit_result is not None:
                return LmfitWrapper.convert_params(self.lmfit_result.params, from_lmfit=True)
            else:
                return ParameterSet({p: np.NaN for p in self.model.all_param_names})

        @property
        def errors(self):
            """Return a dictionary of standard errors for the fitted parameters"""
            if self.lmfit_result is not None:
                lp = self.lmfit_result.params
                return {p: lp[p].stderr for p in lp}
            else:
                return {p: np.NaN for p in self.model.all_param_names}

        @property
        def init_params(self):
            """Return a dictionary of initial values"""
            if self.lmfit_result is not None:
                lp = self.lmfit_result.params
                params = ParameterSet()
                for n in list(self.lmfit_result.params.keys()):
                    if lp[n].vary is False:
                        params[n] = lp[n].value
                    else:
                        ii = self.lmfit_result.var_names.index(n)
                        params[n] = self.lmfit_result.init_vals[ii]
                return params
            else:
                return ParameterSet(
                    {p: np.NaN for p in self.model.all_param_names}
                )

        @property
        def report(self):
            """Return a fit report"""
            if self.lmfit_result is not None:
                return lmfit.fit_report(self.lmfit_result)
            else:
                return 'Chunk failed...'

        @property
        def redchi(self):
            """Return the red. Chi**2 of the fit"""
            if self.lmfit_result is not None:
                return self.lmfit_result.redchi
            else:
                return np.NaN

        @property
        def neval(self):
            """Return the number of evaluations of the fit"""
            if self.lmfit_result is not None:
                return self.lmfit_result.nfev
            else:
                return 0
        
        @property
        def medcnts(self):
            """Return the median counts of the chunk"""
            return np.median(self.chunk.flux)
    
    
    def fit_lsfs(self, lsf_model, params):
        """Fit the lsf model of this initialized fitter object to another
        lsf, defined by the input arguments.
        
        :param lsf_model: The LSF model to fit to.
        :type lsf_model: :class:`LSFModel`
        :param params: The LSF parameters to evaluate the supplied LSF model.
        :type params: :class:`ParameterSet`
        
        :return: The best-fit parameters of the fit.
        :rtype: :class:`ParameterSet`
        """
        
        def fit_func(lmpars, x):
            pars = self.convert_params(lmpars, from_lmfit=True)
            return lsf_y - self.model.lsf_model.eval(x, pars)
        
        x = self.model.lsf_model.generate_x(self.model.osample_factor, self.model.conv_width)
        
        lsf_y = lsf_model.eval(x, params)        
        
        lmpars = self.convert_params(self.model.lsf_model.guess_params(0), to_lmfit=True)
    
        lmfit_result = lmfit.minimize(fit_func, lmpars, args=[x])
        
        return self.convert_params(lmfit_result.params, from_lmfit=True)
        
        
