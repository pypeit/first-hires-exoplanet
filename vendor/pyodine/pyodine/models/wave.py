import numpy as np
from .base import ParameterSet
from .shapes import LinearStaticModel, QuadraticStaticModel


class LinearWaveModel(LinearStaticModel):
    """A linear wavelength model"""
    @staticmethod
    def guess_params(chunk):
        """Make an educated guess of the wavelength parameters for a given chunk
        
        :param chunk: The chunk for which to guess the parameters.
        :type chunk: :class:`Chunk`
        
        :return: The guessed parameters (wavelength zero point and slope).
        :rtype: :class:`ParameterSet`
        """
        p = np.polyfit(chunk.pix, chunk.wave, 1)
        return ParameterSet(intercept=p[1], slope=p[0])
    
    @staticmethod
    def name():
        """The name of the wave model as a string
        
        :return: The wave model name.
        :rtype: str
        """
        return __class__.__name__


class QuadraticWaveModel(QuadraticStaticModel):
    """A 2nd degree polynomial wavelength model"""
    @staticmethod
    def guess_params(chunk):
        """Make an educated guess of the wavelength parameters for a given chunk
        
        :param chunk: The chunk for which to guess the parameters.
        :type chunk: :class:`Chunk`
        
        :return: The guessed parameters (wavelength zero point and slope).
        :rtype: :class:`ParameterSet`
        """
        p = np.polyfit(chunk.pix, chunk.wave, 2)
        return ParameterSet(intercept=p[2], slope=p[1], curvature=p[0])
    
    @staticmethod
    def name():
        """The name of the wave model as a string
        
        :return: The wave model name.
        :rtype: str
        """
        return __class__.__name__


model_index = {
        'LinearWaveModel': LinearWaveModel,
        'QuadraticWaveModel': QuadraticWaveModel
        }