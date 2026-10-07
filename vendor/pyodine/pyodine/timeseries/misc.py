import numpy as np
from ..lib import robust


def reweight(epsilon, alfa, bet, sigma):
    """
    Reweight function for velocity sigmas.
    For this function FWHM = 2.0 * alpha * sigma.
    """
    f = 1.0 / ( 1.0 + (np.abs(epsilon) / (alfa*sigma))**bet )
    return f


def robust_mean(x, axis=None):
    shape = x.shape
    ndim = len(shape)
    # If only one dimension
    if axis is None or ndim == 1:
        if len(shape) > 1:
            x = x.flatten()
        ii = np.where(np.isfinite(x))
        if len(ii[0]) > 0:
            return robust.mean(x[ii])
        else:
            return np.nan
    # Else take the mean along the specified axis
    elif axis == 0 and ndim == 2:
        y = np.zeros(shape[1])
        for k in range(shape[1]):
            ii = np.where(np.isfinite(x[:,k]))
            if len(ii[0]) > 0:
                y[k] = robust.mean(x[ii,k])
            else:
                y[k] = np.nan
        return y
    elif axis == 1 and ndim == 2:
        y = np.zeros(shape[0])
        for k in range(shape[0]):
            ii = np.where(np.isfinite(x[k,:]))
            if len(ii[0]) > 0:
                y[k] = robust.mean(x[k,ii])
            else:
                y[k] = np.nan
        return y
    else:
        # axis > 1 not yet implemented...
        raise


def robust_std(x):
    if len(x.shape) > 1:
        x = x.flatten()
    ii = np.where(np.isfinite(x))
    if len(ii[0]) > 0:
        return robust.std(x[ii])
    else:
        return np.nan