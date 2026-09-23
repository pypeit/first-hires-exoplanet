from ..components import Spectrum

class Model:
    """Abstract base model"""
    param_names = None

    def make_param_list(self, parameter_set):
        """Return a list of parameters
        
        :param parameter_set: The parameter set to transform.
        :type parameter_set: dict, :class:`ParameterSet`
        
        :return: A list of the parameter values from the dictionary.
        :rtype: list
        """
        return [parameter_set[k] for k in self.param_names]


class DynamicModel(Model):
    """A class to incorporate submodels and methods
    
    A DynamicModel must run as an instance.
    
    :param lsf_model: The LSF model to be used.
    :type lsf_model: :class:`LSFModel`
    :param wave_model: The wavelength model to be used.
    :type wave_model: :class:`WavelengthModel`
    :param cont_model: The continuum model to be used.
    :type cont_model: :class:`ContinuumModel`
    :param iodine_atlas: The I2 atlas to be used.
    :type iodine_atlas: :class:`IodineAtlas`
    :param stellar_template: The stellar template to be used. Leave out for the 
        hot-star modelling in the template creation process.
    :type stellar_template: :class:`StellarTemplate`, 
        :class:`StellarTemplate_Chunked`, or None
    :param lsf_array: If fitting with a fixed LSF is desired, pass it here.
    :type lsf_array: :class:`LSF_Array`, or None
    :param osample_factor: Oversampling factor for the model evaluation. 
        Default is an oversampling of 4.
    :type osample_factor: int, or None
    :param conv_width: Number of pixels to evaluate the LSF on (towards either 
        side). Default is 10.
    :type conv_width: int, or None
    """
    def __init__(self, lsf_model, wave_model, cont_model, iodine_atlas, 
                 stellar_template=None, lsf_array=None, 
                 osample_factor=4, conv_width=10.):
        # Sub-models
        self.lsf_model = lsf_model
        self.wave_model = wave_model
        self.cont_model = cont_model
        # Supplementary data
        self.iodine_atlas = iodine_atlas
        self.stellar_template = stellar_template
        # The oversampling factor used
        self.osample_factor = osample_factor
        # The width over which the lsf is evaluated
        self.conv_width=conv_width
        # FIXME: Make sure that if lsf_fixed is given,
        # lsf_model needs to be FixedLSF (how do I do that?!)
        self.lsf_array = lsf_array

    def eval(self, *args, **kwargs):
        raise NotImplementedError

    def eval_spectrum(self, chunk, params, chunk_ind=None, **kwargs):
        """Evaluate the model spectrum for a chunk and set of parameters
        
        :param chunk: The chunk to evaluate over.
        :type chunk: :class:`Chunk`
        :param params: The parameters to use.
        :type params: :class:`ParameterSet`
        :param chunk_ind: The index of the chunk to evaluate.
        :type chunk_ind: int, or None
        
        :return: A spectrum object containing the model.
        :rtype: :class:`Spectrum`
        """
        flux = self.eval(chunk, params, chunk_ind=chunk_ind, **kwargs)
        wave = self.wave_model.eval(chunk.pix, params.filter('wave'), **kwargs)
        cont = self.cont_model.eval(chunk.pix, params.filter('cont'), **kwargs)
        return Spectrum(flux, wave, cont)

    def eval_lsf(self, params, osample_factor=None, conv_width=None):
        """Evaluate and return the LSF for a set of parameters
        
        :param params: The parameters to use.
        :type params: :class:`ParameterSet`
        :param osample_factor: The oversampling factor to use. If None, the 
            model value is used.
        :type osample_factor: int, or None
        :param conv_width: Number of pixels to evaluate the LSF on (towards 
            either side). If None, the model value is used.
        :type conv_width: float, or None
        
        :return: The pixel vector that the LSF is evaluated over.
        :rtype: ndarray[nr_pix]
        :return: The evaluated LSF vector.
        :rtype: ndarray[nr_pix]
        """
        if osample_factor is None:
            osample_factor = self.osample_factor
        if conv_width is None:
            conv_width = self.conv_width
        x = self.lsf_model.generate_x(osample_factor, conv_width)
        
        if self.lsf_array is None:
            lsf = self.lsf_model.eval(x, params.filter('lsf'))
        else:
            lsf = self.lsf_model.eval(self.lsf_array, params.filter('lsf'))
        return x, lsf
    
    
    def guess_params(self, chunk):
        raise NotImplementedError

    @property
    def all_param_names(self):
        """A list of the parameter names of the model
        """
        names = self.param_names
        names = names + ['lsf_' + name for name in self.lsf_model.param_names]
        names = names + ['wave_' + name for name in self.wave_model.param_names]
        names = names + ['cont_' + name for name in self.cont_model.param_names]
        return names


class StaticModel(Model):
    """A class to act as parent class for submodels
    
    It should not be initialized; its methods are static."""
    def __new__(cls, *args, **kwargs):
        raise TypeError('Intended for static use only')

    @classmethod
    def eval(cls, x, params):
        raise NotImplementedError

    @classmethod
    def guess_params(cls, chunk):
        raise NotImplementedError


class ParameterSet(dict):
    """A general set of parameters for a model (a dict with extra methods)"""

    def __getitem__(self, item):
        """The dedicated get-method
        
        :param item: A string corresponding to a parameter key or a key prefix.
        :type item: str
        
        :return: Either a set of parameters corresponding to the prefix, or the 
            parameter value corresponding to the key name.
        :rtype: :class:`ParameterSet`, or value
        """
        if item in self.keys():
            return super().__getitem__(item)
        else:
            return self.filter(prefix=item)

    def filter(self, prefix): #=None):
        """Return a subset of parameters, defined by prefix
        (or something else in the future)
        
        :param prefix: A prefix to filter the parameter keys by (either of 
            'lsf', 'wave' or 'cont' at the moment).
        :type prefix: str
        
        :return: The parameters corresponding to the prefix.
        :rtype: :class:`ParameterSet`
        """
        if prefix is not None:
            new = {k[len(prefix) + 1:]: self[k] for k in self.keys() if k.startswith(prefix + '_')}
            return ParameterSet(new)
        else:
            raise ValueError('No filter keywords set')
    

    def add(self, parameter_set, prefix=''):
        """Add the parameters of another ParameterSet, adding a prefix if set.
        This will override existing items with the same key, if any.
        
        :param parameter_set: A set of parameters to add.
        :type parameter_set: :class:`ParameterSet`
        :param prefix: A prefix to set in front of the parameter keys.
        :type prefix: str
        """
        if prefix:
            prefix += '_'
        for k in parameter_set:
            self[prefix + k] = parameter_set[k]
    
    
    def __str__(self):
        """Return information about the contained parameters
        
        :return: The names and values of the parameters nicely formatted.
        :rtype: str
        """
        string = "<ParameterSet (values: {})>".format(len(self))
        if len(self) > 0:
            fill = max([len(k) for k in self])
            for k in sorted(self):
                string += "\n    {name:<{fill}}  =  {value}".format(
                    fill=fill,
                    name=k,
                    value=self[k]
                )
        return string
