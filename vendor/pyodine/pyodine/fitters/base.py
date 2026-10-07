class Fitter:
    """Abstract base class"""
    def fit(self, *args, **kwargs):
        raise NotImplementedError


class FitResult:
    """Abstract base class"""
    chunk = None
    model = None

    @property
    def params(self):
        """Return a pyodine ParameterSet with the fitted parameters"""
        raise NotImplementedError

    @property
    def errors(self):
        """Return a pyodine ParameterSet with the fitted parameters"""
        raise NotImplementedError

    @property
    def report(self):
        """Return a (byte)string with a report of the fit"""
        raise NotImplementedError

    @property
    def redchi(self):
        """Reduced chi-square"""
        raise NotImplementedError

    @property
    def neval(self):
        """Number of function evaluations"""
        raise NotImplementedError

    @property
    def fitted_spectrum(self):
        """Return the fitted spectrum with modelled wavelength solution and
        continuum
        
        :return: The fitted spectrum.
        :rtype: :class:`Spectrum`
        """
        return self.model.eval_spectrum(self.chunk, self.params, chunk_ind=self.chunk_ind)

    @property
    def fitted_lsf(self):
        """Return x- and y-vector for the fitted LSF
        
        :return: The pixels over which the LSF is sampled.
        :rtype: ndarray[nr_lsf_pix]
        :return: The evaluated LSF values.
        :rtype: ndarray[nr_lsf_pix]
        """
        #osample_factor = 4  # FIXME: This should be dynamic somehow
        return self.model.eval_lsf(self.params)

    @property
    def residuals(self):
        """Calculate residuals
        
        :return: The residuals between the observation and model.
        :rtype: ndarray[nr_pix]
        """
        return self.chunk.flux - self.model.eval(self.chunk, self.params, chunk_ind=self.chunk_ind)
    
    @property
    def medcounts(self):
        """Median counts of the original spectrum chunk"""
        return NotImplementedError
