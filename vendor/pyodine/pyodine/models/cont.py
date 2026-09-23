import numpy as np
from .base import ParameterSet
from .shapes import LinearStaticModel, QuadraticStaticModel


class LinearContinuumModel(LinearStaticModel):
    """A linear continuum model"""
    @staticmethod
    def guess_params(chunk):
        """Make an educated guess of the continuum parameters for a given chunk
        
        :param chunk: The chunk for which to guess the parameters.
        :type chunk: :class:`Chunk`
        
        :return: The guessed parameters (continuum zero point and slope).
        :rtype: :class:`ParameterSet`
        """
        if chunk.cont is not None:
            # Fit a straight line to the continuum and scale to fit the flux
            p = np.polyfit(chunk.pix, chunk.cont, 1)
            intercept = np.median(chunk.flux)
            slope = p[0] / p[1] * intercept
            return ParameterSet(intercept=intercept, slope=slope)
        else:
            # Fit a straight line to the spectral flux
            p = np.polyfit(chunk.pix, chunk.flux, 1)
            return ParameterSet(intercept=p[1], slope=p[0])
    
    @staticmethod
    def name():
        """The name of the continuum model as a string
        
        :return: The continuum model name.
        :rtype: str
        """
        return __class__.__name__


class QuadraticContinuumModel(QuadraticStaticModel):
    """A 2nd degree polynomial continuum model"""
    @staticmethod
    def guess_params(chunk):
        """Make an educated guess of the continuum parameters for a given chunk
        
        :param chunk: The chunk for which to guess the parameters.
        :type chunk: :class:`Chunk`
        
        :return: The guessed polynomial parameters.
        :rtype: :class:`ParameterSet`
        """
        if chunk.cont is not None:
            # Fit a parabola to the continuum and scale to fit the flux
            p = np.polyfit(chunk.pix, chunk.cont, 2)
            return ParameterSet(intercept=p[2], slope=p[1], curvature=p[0])
        else:
            # Fit a parabola to the spectral flux
            p = np.polyfit(chunk.pix, chunk.flux, 2)
            return ParameterSet(intercept=p[2], slope=p[1], curvature=p[0])
    
    @staticmethod
    def name():
        """The name of the continuum model as a string
        
        :return: The continuum model name.
        :rtype: str
        """
        return __class__.__name__


model_index = {
        'LinearContinuumModel': LinearContinuumModel,
        'QuadraticContinuumModel': QuadraticContinuumModel
        }