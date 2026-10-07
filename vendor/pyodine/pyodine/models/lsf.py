import numpy as np
import logging
import sys

from .base import StaticModel, ParameterSet


class LSFModel:
    """Abstract base class for the LSF models
    """
    @staticmethod
    def generate_x(osample_factor, conv_width=10.):
        """Generate a pixel vector to sample the LSF over
        
        :param osample_factor: The oversample factor of the model.
        :type osample_factor: float
        :param conv_width: The number of pixels on either side. Defaults to 10.
        :type conv_width: float, or None
        
        :return: The evaluated pixel vector.
        :rtype: ndarray
        """
        return np.linspace(-conv_width, conv_width, 
                           int(2*conv_width) * int(osample_factor) + 1)


class SingleGaussian(LSFModel, StaticModel):
    """The LSF model of a Single Gaussian
    
    This model has 1 free parameter: The FWHM of the Gaussian.
    
    The class cannot be initialized and all its methods are static or class 
    methods.
    """
    
    param_names = ['fwhm']
    
    param_guess = np.array([2.0])
    
    pars_dict = {}
    
    @staticmethod
    def adapt_LSF(pars_dict):
        """A dummy method implemented to ensure acccordance with high-level
        routines.
        """
        pass

    @classmethod
    def eval(cls, x, params):
        """Evaluate the LSF
        
        :param x: The pixel vector over which to evaluate the LSF.
        :type x: ndarray
        :param params: The LSF parameters.
        :type params: :class:`ParameterSet`
        
        :return: The normalized LSF.
        :rtype: ndarray
        """
        # A single gaussian defined by its FWHM
        y = np.exp(-2.77258872223978123768 * x**2. / params[cls.param_names[0]]**2)
        # Make sure that the sum equals one
        return y / np.sum(y)  # FIXME: Normalize to unit area?

    @classmethod
    def guess_params(cls, chunk):
        """Guess the LSF parameters for a given chunk
        
        This method returns just fixed values, independent of the chunk.
        
        :param chunk: The chunk for which to make the guess.
        :type: :class:`Chunk`
        
        :return: The guessed LSF parameters.
        :rtype: :class:`ParameterSet`
        """
        return ParameterSet(
                {name: guess for name, guess in zip(cls.param_names, cls.param_guess)})
    
    @staticmethod
    def name():
        """The name of the LSF as a string
        
        :return: The LSF name.
        :rtype: str
        """
        return __class__.__name__


class SuperGaussian(LSFModel, StaticModel):
    """The LSF model of a Super Gaussian
    
    The model consists of a central Gaussian-like function, whose exponent and 
    sigma are variable, as well as a satellite Gaussian on either side with 
    variable amplitude but fixed position and sigma.
    
    This class cannot be initialized and all its methods are static or class 
    methods. The positions of all Gaussians and sigmas of the satellites are 
    build in as class variables and can be changed through a dedicated method. 
    This allows the model to be adapted to different instruments.
    """
    param_names = ['sigma', 'exponent', 'left', 'right']
    
    param_guess = np.array([
        1.0, 1.9, 0.2, 0.2
    ])
    
    pars_dict = {
            'positions': np.array([
                    -1.0, 0.0, 1.0]),    
            'satellite_sigmas': np.array([
                    2.5, 2.5])
            }
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    @classmethod
    def adapt_LSF(cls, pars_dict):
        """Adapt the LSF setup to different instruments
        
        This method allows to change the positions of all Gaussians and sigmas 
        of the satellites by handing in a dictionary with desired values.
        
        :param pars_dict: A dictionary with keys 'positions' (length-3 list,
            tuple, or ndarray) and/or 'satellite_sigmas' (length-2 list, tuple, 
            or ndarray).
        :type pars_dict: dict
        """
        if isinstance(pars_dict, dict):
            
            # Update the positions in the parameter dictionary
            if 'positions' in pars_dict.keys() and \
            isinstance(pars_dict['positions'], (list,tuple,np.ndarray)) and \
            len(pars_dict['positions']) == 3:
                if isinstance(pars_dict['positions'], (list,tuple)):
                    cls.pars_dict['positions'] = np.array(pars_dict['positions'])
                else:
                    cls.pars_dict['positions'] = pars_dict['positions']
            
            # Update the satellite sigmas in the parameter dictionary
            if 'satellite_sigmas' in pars_dict.keys() and \
            isinstance(pars_dict['satellite_sigmas'], (list,tuple,np.ndarray)) and \
            len(pars_dict['satellite_sigmas']) == 2:
                if isinstance(pars_dict['satellite_sigmas'], (list,tuple)):
                    cls.pars_dict['satellite_sigmas'] = np.array(pars_dict['satellite_sigmas'])
                else:
                    cls.pars_dict['satellite_sigmas'] = pars_dict['satellite_sigmas']
        else:
            logging.error(pars_dict)
            raise ValueError('Make sure you supply a dictionary with positions' + 
                             ' and/or sigmas, in list, tuple or ndarray of length 3 and/or 2!')

    @classmethod
    def eval(cls, x, params):
        """Evaluate the LSF
        
        :param x: The pixel vector over which to evaluate the LSF.
        :type x: ndarray
        :param params: The LSF parameters.
        :type params: :class:`ParameterSet`
        
        :return: The normalized LSF.
        :rtype: ndarray
        """
        # Include fixed satellite parameters
        a = np.array([params['left'], 1.0, params['right']])
        c = np.array([cls.pars_dict['satellite_sigmas'][0], 
                      params['sigma'], 
                      cls.pars_dict['satellite_sigmas'][1]])
        n = np.array([2.0, params['exponent'], 2.0])

        # Supergauss function
        def func(x):
            xarr = np.transpose([x] * 3)
            f = np.sum(a * np.exp(-0.5 * (np.abs(xarr - cls.pars_dict['positions']) / c) ** n), axis=1)
            # Replace negative values with zero
            f[np.where(f < 0.0)] = 0.0
            return f
        
        try:
            # Evaluate function
            y = func(x)
            
            # Calculate centroid and re-center the LSF
            offset = np.sum(x * y) / np.sum(y)
            y = func(x + offset)
            
            # Make sure that the sum equals one
            y_sum = np.sum(y)
            y = y / y_sum
            
            return y
        except Exception as e:
            logging.error('LSF evaluation failed. Parameters:')
            logging.error(params)
            logging.error('Sum of un-normalized LSF: {}'.format(y_sum))
            raise e

    @classmethod
    def guess_params(cls, chunk):
        """Guess the LSF parameters for a given chunk
        
        This method returns just fixed values, independent of the chunk.
        
        :param chunk: The chunk for which to make the guess.
        :type: :class:`Chunk`
        
        :return: The guessed LSF parameters.
        :rtype: :class:`ParameterSet`
        """
        return ParameterSet(
                {name: guess for name, guess in zip(cls.param_names, cls.param_guess)})
    
    @staticmethod
    def name():
        """The name of the LSF as a string
        
        :return: The LSF name.
        :rtype: str
        """
        return __class__.__name__


class MultiGaussian(LSFModel, StaticModel):
    """The LSF model of a Multi Gaussian
    
    The model consists of a central, completely fixed Gaussian, and 5 satellite
    Gaussian both to the left and right. The positions and sigmas of the 
    satellites are fixed, but their amplitudes are the 10 free parameters.
    
    This class cannot be initialized and all its methods are static or class 
    methods. The positions and sigmas of all Gaussians are build in as class
    variables and can be changed through dedicated methods. This allows the 
    model to be adapted to different instruments.
    """
    param_names = [
        'left_5', 'left_4', 'left_3', 'left_2', 'left_1',
        'right1', 'right2', 'right3', 'right4', 'right5',
    ]
    
    param_guess = np.array([
        0.1, 0.2, 0.3, 0.5, 0.7,
        0.7, 0.5, 0.3, 0.2, 0.1
    ])
    
    pars_dict = {
            'positions': np.array([
                    -2.9, -2.5, -1.9, -1.4, -1.0, 
                    0.0, 
                    1.0, 1.4, 1.9, 2.5, 2.9]),
            'sigmas': np.array([
                    0.9, 0.9, 0.9, 0.9, 0.9, 
                    0.6, 
                    0.9, 0.9, 0.9, 0.9, 0.9])
            }
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    @classmethod
    def adapt_LSF(cls, pars_dict):
        """Adapt the LSF setup to different instruments
        
        This method allows to change the positions and sigmas of all Gaussians
        by handing in a dictionary with desired values.
        
        :param pars_dict: A dictionary with keys 'positions' and/or 'sigmas' 
            (for both: length-11 list, tuple, or ndarray).
        :type pars_dict: dict
        """
        if isinstance(pars_dict, dict):
            
            # Update the positions in the parameter dictionary
            if 'positions' in pars_dict.keys() and \
            isinstance(pars_dict['positions'], (list,tuple,np.ndarray)) and \
            len(pars_dict['positions']) == 11:
                if isinstance(pars_dict['positions'], (list,tuple)):
                    cls.pars_dict['positions'] = np.array(pars_dict['positions'])
                else:
                    cls.pars_dict['positions'] = pars_dict['positions']
            
            # Update the sigmas in the parameter dictionary
            if 'sigmas' in pars_dict.keys() and \
            isinstance(pars_dict['sigmas'], (list,tuple,np.ndarray)) and \
            len(pars_dict['sigmas']) == 11:
                if isinstance(pars_dict['sigmas'], (list,tuple)):
                    cls.pars_dict['sigmas'] = np.array(pars_dict['sigmas'])
                else:
                    cls.pars_dict['sigmas'] = pars_dict['sigmas']
        else:
            logging.error(pars_dict)
            raise ValueError('Make sure you supply a dictionary with positions' + 
                             ' and/or sigmas, in list, tuple or ndarray of length 11!')

    @classmethod
    def eval(cls, x, params):
        """Evaluate the LSF
        
        :param x: The pixel vector over which to evaluate the LSF.
        :type x: ndarray
        :param params: The LSF parameters.
        :type params: :class:`ParameterSet`
        
        :return: The normalized LSF.
        :rtype: ndarray
        """
        # Convert input dict to list
        params = np.array([params[k] for k in cls.param_names])

        # Set up parameter vectors, including central gaussian
        a = np.array([
            params[0], params[1], params[2], params[3], params[4],
            1.0,
            params[5], params[6], params[7], params[8], params[9],
        ])

        # Multigauss function
        def func(x):
            xarr = np.repeat([x], len(a), axis=0)
            f = np.sum(a * np.exp(-0.5 * ((np.transpose(xarr) - cls.pars_dict['positions']) 
                                            / cls.pars_dict['sigmas'])**2.), axis=1)
            f[np.where(f < 0.0)] = 0.0
            return f
        
        try:
            # Evaluate function and find centroid
            y = func(x)
            
            # Calculate centroid and re-center the LSF
            offset = np.sum(x * y) / np.sum(y)
            y = func(x + offset)
            
            # Make sure that the sum equals one
            y_sum = np.sum(y)
            y_norm = y / y_sum
            
            if len(np.where(np.isnan(y))[0]) > 0:
                logging.debug('NaNs detected in LSF. Parameters:')
                logging.debug(params)
                logging.debug('Sum of un-normalized LSF: {}'.format(y_sum))
                logging.debug('Unnormalized LSF:')
                logging.debug(y)
            
            return y_norm
        except Exception as e:
            logging.error('LSF evaluation failed. Parameters:')
            logging.error(params)
            logging.error('Sum of un-normalized LSF: {}'.format(y_sum))
            raise e
            

    @classmethod
    def guess_params(cls, chunk):
        """Guess the LSF parameters for a given chunk
        
        This method returns just fixed values, independent of the chunk.
        
        :param chunk: The chunk for which to make the guess.
        :type: :class:`Chunk`
        
        :return: The guessed LSF parameters.
        :rtype: :class:`ParameterSet`
        """
        return ParameterSet(
                {name: guess for name, guess in zip(cls.param_names, cls.param_guess)})
    
    @staticmethod
    def name():
        """The name of the LSF as a string
        
        :return: The LSF name.
        :rtype: str
        """
        return __class__.__name__


class MultiGaussian_Lick(LSFModel, StaticModel):
    """The LSF model of a Multi Gaussian as used in Lick (algorithm employed
    as in dop-code by D. Fisher, Yale University, no re-centering of the LSF)
    
    The model consists of a central, completely fixed Gaussian, and 5 satellite
    Gaussian both to the left and right. The positions and sigmas of the 
    satellites are fixed, but their amplitudes are the 10 free parameters.
    
    This class cannot be initialized and all its methods are static or class 
    methods. The positions and sigmas of all Gaussians are build in as class
    variables and can be changed through dedicated methods. This allows the 
    model to be adapted to different instruments.
    """
    
    param_names = [
        'left_5', 'left_4', 'left_3', 'left_2', 'left_1', 
        'right1', 'right2', 'right3', 'right4', 'right5'
        ]
    
    param_guess = np.array([
        0.1, 0.2, 0.3, 0.4, 0.5,
        0.5, 0.4, 0.3, 0.2, 0.1
    ])
    
    pars_dict = {
            'positions': np.array([
                    -2.4, -2.1, -1.6, -1.1, -0.6,
                    0.0,
                    0.6, 1.1, 1.6, 2.1, 2.4]),
            'sigmas': np.array([
                    0.3, 0.3, 0.3, 0.3, 0.3,
                    0.4,
                    0.3, 0.3, 0.3, 0.3, 0.3])
            }
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    @classmethod
    def adapt_LSF(cls, pars_dict):
        """Adapt the LSF setup to different instruments
        
        This method allows to change the positions and sigmas of all Gaussians
        by handing in a dictionary with desired values.
        
        :param pars_dict: A dictionary with keys 'positions' and/or 'sigmas' 
            (for both: length-11 list, tuple, or ndarray).
        :type pars_dict: dict
        """
        if isinstance(pars_dict, dict):
            
            # Update the positions in the parameter dictionary
            if 'positions' in pars_dict.keys() and \
            isinstance(pars_dict['positions'], (list,tuple,np.ndarray)) and \
            len(pars_dict['positions']) == 11:
                if isinstance(pars_dict['positions'], (list,tuple)):
                    cls.pars_dict['positions'] = np.array(pars_dict['positions'])
                else:
                    cls.pars_dict['positions'] = pars_dict['positions']
            
            # Update the sigmas in the parameter dictionary
            if 'sigmas' in pars_dict.keys() and \
            isinstance(pars_dict['sigmas'], (list,tuple,np.ndarray)) and \
            len(pars_dict['sigmas']) == 11:
                if isinstance(pars_dict['sigmas'], (list,tuple)):
                    cls.pars_dict['sigmas'] = np.array(pars_dict['sigmas'])
                else:
                    cls.pars_dict['sigmas'] = pars_dict['sigmas']
        else:
            logging.error(pars_dict)
            raise ValueError('Make sure you supply a dictionary with positions' + 
                             ' and/or sigmas, in list, tuple or ndarray of length 11!')
    
    @classmethod
    def eval(cls, x, params):
        """Evaluate the LSF
        
        :param x: The pixel vector over which to evaluate the LSF.
        :type x: ndarray
        :param params: The LSF parameters.
        :type params: :class:`ParameterSet`
        
        :return: The normalized LSF.
        :rtype: ndarray
        """
        # Convert input dict to list
        params = np.array([params[k] for k in cls.param_names])

        # Set up parameter vectors, including central gaussian
        a = np.array([
            params[0], params[1], params[2], params[3], params[4],
            1.0,
            params[5], params[6], params[7], params[8], params[9],
        ])
        # This is from the cf's
        b = cls.pars_dict['positions']
        c = cls.pars_dict['sigmas']

        # Gaussian function, computed as in Lick dop code
        def func(x, cntr=None):
            y = np.zeros(len(x))
            if cntr is None:
                cntr = 0.0
            cen = 0. - cntr
            # First the central Gaussian
            cent_wid = c[5] * 5.
            if cent_wid < 2.:
                cent_wid = 2.  # define the central gaussian over this restricted pixel interval
            xx = np.where((x >= cen-cent_wid) & (x <= cen+cent_wid))
            y[xx] = np.exp(-0.5 * ((x[xx] - cen) / c[5])**2.)
            
            # Now the satellite Gaussians
            for i in range(len(a)):
                if i != 5:
                    cen = b[i] - cntr
                    gd_range = 5. * c[i]
                    xx = np.where((x >= cen-gd_range) & (x <= cen+gd_range))
                    y[xx] += a[i] * np.exp(-0.5 * ((x[xx] - cen) / c[i])**2.)
            # Normalize
            #dx = (x[-1] - x[0]) / len(x)
            #y[np.where(y < 0.0)] = 0.0
            
            return y / np.sum(y) #(dx * np.sum(y))

        # Evaluate the function
        y = func(x)
        # Shift to Center? Not doing it gives better results for Lick spectra!
        """
        fwhm = 0.5 * np.max(y)
        x2 = np.where((y >= fwhm) & (np.abs(x) < 6.)) # peak points +/- from center
        
        if len(x2[0]) >= 3:
            #print('>=3')
            dd = np.where(y[x2] == np.max(y[x2]))
            ndd = len(dd[0])
            if ndd <= 2:
                #print('ndd<=2')
                dd = dd[0][0]
            if ndd > 2:
                #print('ndd>2')
                dd = dd[0][int(ndd/2.)]
            cntr = x[x2[0][dd]]
            #print(cntr)
            if abs(cntr) >= 0.1 and abs(cntr) < 1.2:
                #print('Shift LSF.')
                y = func(x, cntr=cntr)
                #y = y / np.sum(np.abs(y))#(dx * np.sum(y))
                #print('New: ', x[np.argmax(y)])
        """
        return y
    
    @classmethod
    def guess_params(cls, chunk):
        """Guess the LSF parameters for a given chunk
        
        This method returns just fixed values, independent of the chunk.
        
        :param chunk: The chunk for which to make the guess.
        :type: :class:`Chunk`
        
        :return: The guessed LSF parameters.
        :rtype: :class:`ParameterSet`
        """
        return ParameterSet(
                {name: guess for name, guess in zip(cls.param_names, cls.param_guess)})
    
    @staticmethod
    def name():
        """The name of the LSF as a string
        
        :return: The LSF name.
        :rtype: str
        """
        return __class__.__name__


class HermiteGaussian(LSFModel, StaticModel):
    """The LSF model of a Single Gaussian perturbed by Hermite polynomials
    
    This model has up to 10 free parameters: The FWHM of the Gaussian, which 
    also regulates the width of the Hermite polynomials, and 9 weights for the
    9 Hermite polynomials from order 1 to 9. (Order 0 is always 1 so that the
    product of all does not become 0!)
    
    The class cannot be initialized and all its methods are static or class 
    methods.
    
    Contributed by Ayk Jessen (Landessternwarte Heidelberg, 2022).
    """
    
    param_names = [
        'fwhm',
        'weight_1', 'weight_2', 'weight_3', 'weight_4', 'weight_5',
        'weight_6', 'weight_7', 'weight_8', 'weight_9'
        ]
 
    # Default parameter guess: Only Hermite orders 3 - 8 are used
    param_guess = np.array([
            1.07, 
            0.0, 0.0, 0.000000001, -0.0000000001, -0.00000000001, 
            0.000000000001, -0.000000000001, -0.000000000001, 0.0
            ])
    
    # Parameter and whether to use it. If set to 0, that parameter will always
    # be 0 in the evaluation.
    pars_dict = {'fwhm': 1, 'weight_1': 0, 'weight_2': 0, 'weight_3': 1, 
                 'weight_4': 1, 'weight_5': 1, 'weight_6': 1, 'weight_7': 1, 
                 'weight_8': 1, 'weight_9': 0}
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                            format='%(message)s')
    
    @classmethod
    def adapt_LSF(cls, pars_dict):
        """Adapt the LSF setup to different instruments
        
        This method allows to specify which parameters (specifically Hermite
        degrees) actually to use. Those set to 0 through the pars_dict will
        always be 0 in the evaluation!
        
        :param pars_dict: A dictionary with keys as in param_names, and either
            0 for non-usage of that parameter, or 1.
        :type pars_dict: dict
        """
        if isinstance(pars_dict, dict):
            # Update which parameters will actually be used, and which will
            # be held at 0
            for key in pars_dict.keys():
                if key in cls.pars_dict and pars_dict[key] in (0,1):
                    cls.pars_dict[key] = pars_dict[key]
        else:
            logging.error(pars_dict)
            raise ValueError('Make sure you supply a dictionary with keys as in' + 
                             ' param_names, and values either 0 or 1!')

    @classmethod
    def eval(cls, x, params):
        """Evaluate the LSF
        
        :param x: The pixel vector over which to evaluate the LSF.
        :type x: ndarray
        :param params: The LSF parameters.
        :type params: :class:`ParameterSet`
        
        :return: The normalized LSF.
        :rtype: ndarray
        """
        # Convert input dict to numpy array: Multiply each parameter by the
        # pars_dict value, non-activated parameters always become 0 then.
        params = np.array([params[k]*cls.pars_dict[k] for k in cls.param_names])
        
        # Using sigma instead of FWHM
        params[0] *= 1/(2.355)

        # HermiteGauss function
        def func(x):
            x_hermitearg = np.divide(x, params[0]) 
            Hermite_polynomial = np.polynomial.hermite.Hermite(
                    (1, params[1], params[2], params[3], params[4], 
                     params[5], params[6], params[7], params[8], params[9])
                    )
            # Do I even need to divide by the norm? LSF is normalized by sum
            # later anyway
            Norm_hermite =  np.sum(params[1:]) + 1
            exp = np.exp(-np.power(x, 2) / (2 * np.power(params[0], 2))) / Norm_hermite
            f = np.multiply(Hermite_polynomial(x_hermitearg), exp)
            
            #TODO: Constrain not to become 0?
            return f
        
        try:
            # Evaluate function and find centroid
            y = func(x)
            
            # Calculate centroid and re-center the LSF
            offset = np.sum(x * y) / np.sum(y)
            y = func(x + offset)
            
            # Make sure that the sum equals one
            y_sum = np.sum(y)
            y_norm = y / y_sum
            
            if len(np.where(np.isnan(y))[0]) > 0:
                logging.debug('NaNs detected in LSF. Parameters:')
                logging.debug(params)
                logging.debug('Sum of un-normalized LSF: {}'.format(y_sum))
                logging.debug('Unnormalized LSF:')
                logging.debug(y)
            
            return y_norm
        
        except Exception as e:
            logging.error('LSF evaluation failed. Parameters:')
            logging.error(params)
            logging.error('Sum of un-normalized LSF: {}'.format(y_sum))
            raise e
            
    @classmethod
    def guess_params(cls, chunk):
        """Guess the LSF parameters for a given chunk
        
        This method returns just fixed values, independent of the chunk.
        
        :param chunk: The chunk for which to make the guess.
        :type: :class:`Chunk`
        
        :return: The guessed LSF parameters.
        :rtype: :class:`ParameterSet`
        """
        return ParameterSet(
                {name: guess for name, guess in zip(cls.param_names, cls.param_guess)})
    
    @staticmethod
    def name():
        """The name of the LSF as a string
        
        :return: The LSF name.
        :rtype: str
        """
        return __class__.__name__  


class FixedLSF(LSFModel, StaticModel):
    """The LSF model for a fixed LSF
    
    3 free parameters: Amplitude of the LSF, order and pixel0 of the respective
        LSF chunk.
    """
    param_names = ['amplitude', 'order', 'pixel0']   
    
    pars_dict = {}
    
    @staticmethod
    def adapt_LSF(pars_dict):
        """A dummy method implemented to ensure acccordance with high-level
        routines.
        """
        pass
    
    @staticmethod
    def eval(x, params):
        """Evaluate the LSF
        
        The parameter x serves here not as the x-vector over which to compute 
        the LSF, but instead is the LSF itself! Small hack for the time being.
        
        :param x: The evaluated LSF vector for the given chunk.
        :type x: ndarray
        :param params: The LSF parameters.
        :type params: :class:`ParameterSet`
        
        :return: The normalized LSF.
        :rtype: ndarray
        """
        
        lsf_fixed = x[params['order'], params['pixel0']]
        return lsf_fixed * params['amplitude']
        

    @staticmethod
    def guess_params(chunk):
        """Guess the LSF parameters for a given chunk
        
        :param chunk: The chunk for which to make the guess.
        :type: :class:`Chunk`
        
        :return: The guessed LSF parameters.
        :rtype: :class:`ParameterSet`
        """
        return ParameterSet(
                amplitude=1., order=chunk.order, pixel0=chunk.abspix[0]
                )
    
    @staticmethod
    def name():
        """The name of the LSF as a string
        
        :return: The LSF name.
        :rtype: str
        """
        return __class__.__name__


model_index = {
        'SingleGaussian': SingleGaussian,
        'SuperGaussian': SuperGaussian,
        'MultiGaussian': MultiGaussian,
        'MultiGaussian_Lick': MultiGaussian_Lick,
        'HermiteGaussian': HermiteGaussian,
        'FixedLSF': FixedLSF
        }


class LSF_Array:
    """A convenience class to enable modeling of a fixed LSF
    
    Needs to be supplied with a full lsf_array, and respective orders
    and pixels arrays describing the position of each LSF chunk.
    
    Args:
        lsf_array (ndarray[nr_chunks_total,nr_lsf_pix]): An array of 
            evaluated LSFs for all chunks in all orders.
        orders (ndarray[nr_chunks_total]): An array of order numbers for all
            chunks.
        pixels (ndarray[nr_chunks_total]): An array of pixel numbers for all
            chunks.
    """
    def __init__(self, lsf_array, orders, pixels):
        self.lsf_array = lsf_array
        self.orders = orders
        self.pixels = pixels
        
        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                                format='%(message)s')
    
    def __getitem__(self, args):
        """The dedicated get-method
        
        Return the LSF of a desired chunk.
        
        Args:
            args (tuple): Should be a 2-entry tuple with desired order and
                pixel number.
        
        Return:
            ndarray[nr_lsf_pix]: An array containing the LSF of the desired
                chunk.
        """
        if len(args) != 2:
            raise IndexError('Two indices expected, but got {}.'.format(len(args)))
        order, pixel = args
        try:
            return self.lsf_array[np.where((self.orders==order) & (self.pixels==pixel))[0][0]]
        except Exception as e:
            logging.error('LSF retrieval failed. Arguments:')
            logging.error(args)
            raise e