from .base import StaticModel


class LinearStaticModel(StaticModel):
    param_names = ['intercept', 'slope']

    @staticmethod
    def eval(x, params):
        return params['intercept'] + params['slope'] * x

    @staticmethod
    def guess_params(chunk):
        raise NotImplementedError


class QuadraticStaticModel(StaticModel):
    param_names = ['intercept', 'slope', 'curvature']
    
    @staticmethod
    def eval(x, params):
        # Is the order of the parameters correct?
        return params['intercept'] + params['slope'] * x + params['curvature'] * x**2
    
    @staticmethod
    def guess_params(chunk):
        raise NotImplementedError