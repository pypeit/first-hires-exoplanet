import numpy as np
from scipy.interpolate import interp1d
from scipy.special import erfinv
import os
import json
import logging
import logging.config

_c = 299792458  # m/s


def setup_logging(config_file=None, level=logging.INFO, error_log=None,
                  info_log=None, quiet=False):
    """Setup logging configuration, following
    https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/
    
    :param config_file: Pathname to a logging configuration file, in json-
        format. If None is given, fall back on the basic logging configuration.
    :type config_file: str, or None
    :param level: The logging level, e.g. logging.DEBUG or logging.WARNING. 
        Defaults to logging.INFO.
    :type level: logging level
    :param error_log: Pathname to a file where to log errors. Defaults to None
        (no logging to file).
    :type error_log: str, or None
    :param info_log:Pathname to a file where to log info. Defaults to None
        (no logging to file).
    :type info_log: str, or None
    :param quiet: Whether or not to print info to the terminal (errors will 
        always be printed). Defaults to False.
    :type quiet: bool
    """
    
    log_handlers = []
    
    # Is the configuration file (json format) there?
    if isinstance(config_file, str) and os.path.exists(config_file):
        try:
            # Try and load the configuration dictionary
            with open(config_file, 'rt') as f:
                config = json.load(f)
            
            # If you want errors logged
            if isinstance(error_log, str):
                
                # Create directory structure if non-existent yet
                error_log_dir = os.path.dirname(error_log)
                if error_log_dir != '' and not os.path.exists(error_log_dir):
                    os.makedirs(error_log_dir)
                
                config['handlers']['error_file_handler']['filename'] = error_log
                log_handlers.append('error_file_handler')
            else:
                del config['handlers']['error_file_handler']
            
            # If you want info logged
            if isinstance(info_log, str):
                
                # Create directory structure if non-existent yet
                info_log_dir = os.path.dirname(info_log)
                if info_log_dir != '' and not os.path.exists(info_log_dir):
                    os.makedirs(info_log_dir)
                
                config['handlers']['info_file_handler']['filename'] = info_log
                config['handlers']['info_file_handler']['level'] = level
                log_handlers.append('info_file_handler')
            else:
                del config['handlers']['info_file_handler']
            
            # If you want info printed or not 
            # (if not, errors will still be printed!)
            if quiet:
                #del config['handlers']['console']
                config['handlers']['console']['level'] = logging.ERROR
            else:
                config['handlers']['console']['level'] = level
            log_handlers.append('console')
            
            config['root']['handlers'] = log_handlers
            config['root']['level'] = level
            
            logging.config.dictConfig(config)
            
        except Exception as e:
            logging.basicConfig(level=level)
            logging.error('Logger could not be configured from config file', 
                          exc_info=True)
    else:
        if not quiet:
            logging.basicConfig(level=level)
        else:
            logging.basicConfig(level=logging.ERROR)
        logging.warning('No config file supplied or config file not found.')


def return_existing_files(filenames):
    """Check the input filenames and only return those that actually exist
    
    :param filenames: The full pathname(s) of the input file(s).
    :type filenames: str, list, ndarray, tuple
    
    :return: A list of existing files.
    :rtype: list
    :return: A list of non-existing files.
    :rtype: list
    """
    
    if isinstance(filenames, str):
        filenames = [filenames]
    
    if isinstance(filenames, (list,tuple,np.ndarray)):
        existing_files = [f for f in filenames if os.path.isfile(f)]
        bad_files = [f for f in filenames if not os.path.isfile(f)]
    else:
        raise ValueError('No files supplied! (Either of str, list, ndarray, tuple)')
    
    return existing_files, bad_files


def findwave(multiorder, wavelength):
    """I think this is also not required anymore.
    """
    for i, spec in enumerate(multiorder):
        if spec.wave[0] < wavelength < spec.wave[-1]:
            return i
    raise ValueError('Could not find a corresponding order for wavelength %s Å' % wavelength)


def osample(x, factor):
    """Linear oversampling of a vector
    
    :param x: Input vector.
    :type x: ndarray[nr_pix_in]
    :param factor: Oversampling factor.
    :type factor: int
    
    :return: Resulting oversampled vector.
    :rtype: ndarray[nr_pix_out]
    """
    n = int(np.round(factor * (len(x) - 1) + 1))
    return np.linspace(x[0], x[-1], n)


def normgauss(x, fwhm):
    """A normalized Gaussian, defined by its FWHM
    
    :param x: Input vector to sample over.
    :type x: ndarray[nr_pix]
    :param fwhm: Full-width half-maximum of the Gaussian.
    :type fwhm: float
    
    :return: The normalized Gaussian.
    :rtype: ndarray[nr_pix]
    """
    y = np.exp(-2.77258872223978123768 * x**2 / fwhm**2)
    # Make sure that the sum equals one
    return y / np.sum(y)


def rebin(wold, sold, wnew):
    """Interpolates OR integrates a spectrum onto a new wavelength scale, depending
    on whether number of pixels per angstrom increases or decreases. Integration
    is effectively done analytically under a cubic spline fit to old spectrum.

    Ported to from rebin.pro (IDL) to Python by Frank Grundahl (FG).
    Original program written by Jeff Valenti.

    IDL Edit History:
    ; 10-Oct-90 JAV Create.
    ; 22-Sep-91 JAV Translated from IDL to ANA.
    ; 27-Aug-93 JAV Fixed bug in endpoint check: the "or" was essentially an "and".
    ; 26-Aug-94 JAV Made endpoint check less restrictive so that identical old and
    ;       new endpoints are now allowed. Switched to new Solaris library
    ;       in call_external.
    ; Nov01 DAF eliminated call_external code; now use internal idl fspline
    ; 2008: FG replaced fspline with spline
    
    :param wold: Input wavelength vector.
    :type wold: ndarray[nr_pix_in]
    :param sold: Input spectrum to be binned.
    :type sold: ndarray[nr_pix_in]
    :param wnew: New wavelength vector to bin to.
    :type wnew: ndarray[nr_pix_out]
    
    :return: Newly binned spectrum.
    :rtype: ndarray[nr_pix_out]
    """

    def idl_rebin(a, shape):
        sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
        return a.reshape(sh).mean(-1).mean(1)

    # Determine spectrum attributes.
    nold   = np.int32(len(wold))            # Number of old points
    nnew   = np.int32(len(wnew))            # Number of new points
    psold  = (wold[nold - 1] - wold[0]) / (nold - 1)  # Old pixel scale
    psnew  = (wnew[nnew - 1] - wnew[0]) / (nnew - 1)  # New pixel scale

    # Verify that new wavelength scale is a subset of old wavelength scale.
    if (wnew[0] < wold[0]) or (wnew[nnew - 1] > wold[nold - 1]):
        logging.warning('New wavelength scale not subset of old.')

    # Select integration or interpolation depending on change in dispersion.

    if psnew < psold:

        # Pixel scale decreased ie, finer pixels
        # Interpolating onto new wavelength scale.
        dum  = interp1d(wold, sold)  # dum  = interp1d( wold, sold, kind='cubic' ) # Very slow it seems.
        snew = dum(wnew)

    else:

        # Pixel scale increased ie more coarse
        # Integration under cubic spline - changed to interpolation.

        xfac = np.int32(psnew / psold + 0.5) 	# pixel scale expansion factor

        # Construct another wavelength scale (W) with a pixel scale close to that of
        # the old wavelength scale (Wold), but with the additional constraint that
        # every XFac pixels in W will exactly fill a pixel in the new wavelength
        # scale (Wnew). Optimized for XFac < Nnew.

        dw   = 0.5 * (wnew[2:] - wnew[:-2])        # Local pixel scale

        pre  = np.float(2.0 * dw[0] - dw[1])
        post = np.float(2.0 * dw[nnew - 3] - dw[nnew - 4])

        dw   = np.append(dw[::-1], pre)[::-1]
        dw   = np.append(dw, post)
        w    = np.zeros((nnew, xfac), dtype='float')

        # Loop thru subpixels
        for i in range(0, xfac):
            w[:, i] = wnew + dw * (np.float(2 * i + 1) / (2 * xfac) - 0.5)   # pixel centers in W

        nig  = nnew * xfac			# Elements in interpolation grid
        w    = np.reshape(w, nig)   # Make into 1-D

        # Interpolate old spectrum (Sold) onto wavelength scale W to make S. Then
        # sum every XFac pixels in S to make a single pixel in the new spectrum
        # (Snew). Equivalent to integrating under cubic spline through Sold.

        # dum    = interp1d( wold, sold, kind='cubic' ) # Very slow!
        # fill_value in interp1d added to deal with w-values just outside the interpolation range
        dum    = interp1d(wold, sold, fill_value="extrapolate")
        s      = dum(w)
        s      = s / xfac				# take average in each pixel
        sdummy = s.reshape(nnew, xfac)
        snew   = xfac * idl_rebin(sdummy, [nnew, 1])
        snew   = np.reshape(snew, nnew)

    return snew


# Chauvenet criterion, as implemented in dop code
def chauvenet_criterion(residuals, iterate=True):
    """Chauvenet's criterion, as implemented in dop-code by D. Fisher (Yale
    University):
    Find elements that lie too far away from the others.
    Updated: nan-values in the residuals are immediately marked as bad ones.
    
    :param residuals: Input vector with residuals to analyze.
    :type residuals: ndarray[nr_pix]
    :param iterate: Iteratively continue to throw out elements? Defaults to 
        True.
    :type iterate: bool
    
    :return: A mask of same length as input residuals array, with ones where 
        the criterion was passed, and zeros where it failed.
    :rtype: ndarray[nr_pix]
    :return:  An array with the indices where mask is True.
    :rtype: ndarray[nr_good]
    :return:  An array with the indices where mask is False.
    :rtype: ndarray[nr_bad]
    """
    # Find the finite values
    fin = np.where(np.isfinite(residuals))
    N = len(fin[0])
    if N == 0:
        raise ValueError('No finite values in the residuals!')
    # This only makes sense if there's at least two data points
    if N == 1:
        raise ValueError('Chauvenets criterion cant be applied to only one datapoint!')
    
    mean = np.nanmean(residuals)
    rms  = np.nanstd(residuals)
    
    # Create a mask of points that pass Chauvenet's criterion
    # (inverse error function of (1 - 0.5/N))
    mask = np.abs(residuals-mean) < 1.41421356 * rms * erfinv(1. - 0.5/N)

    # Iterate if desired, and throw out points until criterion is satisfied
    if iterate == True:
        indx = np.where(mask == True)
        if len(indx[0]) == 0:
            raise ValueError('All data have failed Chauvenets criterion - ' + \
                             'check that your model is appropriate for these data.')

        if len(indx[0]) < N:
            iter_mask, iter_mask_true, iter_mask_false = chauvenet_criterion(residuals[indx])
            mask[indx] = mask[indx] & iter_mask
            
    return mask, np.where(mask == True), np.where(mask == False)


# Smooth over given chunk LSFs
def smooth_lsf(chunk_arr, pixel_avg, order_avg, order_dist, fit_results, redchi2=None, 
               osample=None, lsf_conv_width=None):
    """Smooth over all chunk LSFs, with given 'radii' in orders & order pixels. 
    LSFs are weighted by red. Chi2 values of modeled chunks.
    Implemented with great parallels to the dop-code routine dop_psf_smooth.pro
    (D. Fisher, Yale University).
    
    NOTE: This currently only works if chunks are evenly distributed, i.e. same
    number of chunks within each order!!!
    
    :param chunk_arr: An array of chunks.
    :type chunk_arr: :class:`ChunkArray`
    :param pixel_avg: Pixels to smooth over in dispersion direction.
    :type pixel_avg: int
    :param order_avg: Orders to smooth over in cross-dispersion direction.
    :type order_avg: int
    :param order_dist: Approximate pixel distance between orders.
    :type order_dist: int
    :param fit_results: The fit result which holds the LSFs to smooth over and 
        red. Chi2s.
    :type fit_results: :class:`LmfitResult`
    :param redchi2: An array of red. Chi2 values to use instead of the ones 
        from the fit result.
    :type redchi2: ndarray[nr_chunks], or None
    :param osample: Oversampling to use; if not given, use the one from the fit 
        results model (default).
    :type osample: int, or None
    :param lsf_conv_width: Number of pixels to evaluate the LSF on (towards 
        either side). If None, use the one from the fit results model (default).
    :type lsf_conf_width: int, or None
    
    :return: An array with the smoothed LSFs for all chunks.
    :rtype: ndarray[nr_chunks, nr_pix]
    """
    # First we set up an easy chunk numpy-array with orders and center pixels
    # to find relevant smoothing indices for each chunk later
    chunk_arr_simple = []
    chunk_lsfs = []
    for i, chunk in enumerate(chunk_arr):
        chunk_arr_simple += [[chunk.order, chunk.abspix[int(len(chunk)/2.)]]]
        if lsf_conv_width is None:
            lsf_conv_width = fit_results[0].model.conv_width
        # Also we already evaluate all individual LSFs here, then we can later
        # simply pick them
        if osample is None:
            x_lsf, lsf = fit_results[i].fitted_lsf
        else:
            x_lsf = fit_results[i].model.lsf_model.generate_x(osample, conv_width=lsf_conv_width)
            lsf = fit_results[i].model.lsf_model.eval(x_lsf, fit_results[i].params.filter('lsf'))
        chunk_lsfs.append(lsf)
        
    chunk_arr_simple = np.array(chunk_arr_simple)
    chunk_lsfs = np.array(chunk_lsfs)
    
    # Now check if redchi2 is given. If not, use the redchi2 from fit_results
    if redchi2 is None or len(redchi2) is not len(chunk_arr_simple):
        redchi2 = np.zeros(len(chunk_arr_simple))
        for i, result in enumerate(fit_results):
            redchi2[i] = fit_results[i].redchi
    else:
        redchi2 = np.array(redchi2)
    
    # Now loop over chunks
    lsfs_smoothed = []
    for i, chunk in enumerate(chunk_arr):
        # These are the chunk indices to smooth over
        xsm = np.where( (chunk_arr_simple[:,0] <= chunk.order + order_avg) &
                        (chunk_arr_simple[:,0] >= chunk.order - order_avg) &
                        (chunk_arr_simple[:,1] <= chunk.abspix[int(len(chunk)/2.)] + pixel_avg) &
                        (chunk_arr_simple[:,1] >= chunk.abspix[int(len(chunk)/2.)] - pixel_avg))
        
        if len(xsm[0]) == 0:
            logging.warning('Chunk {}: No chunks to smooth over. Using input chunk.'.format(i))
            # Use lsf from fit result for that chunk
            lsfs_smoothed.append(chunk_lsfs[i])
        
        else:
            # Compute chunk weights for order averaging
            ord_sep = np.abs(chunk_arr_simple[xsm[0],0] - chunk.order) * order_dist
            ord_wt = 1. / np.sqrt(1. * ord_sep/len(chunk) + 1.)
            pix_sep = np.abs(chunk_arr_simple[xsm[0],1] - chunk.abspix[int(len(chunk)/2.)])
            pix_wt = 1. / np.sqrt(pix_sep/len(chunk) + 1.)
            
            wt = 1. / redchi2[xsm[0]]
            wt = ord_wt * pix_wt * wt
            
            # Check for nans
            nan_ind = np.argwhere(np.isnan(wt))
            if len(nan_ind) > 0:
                wt[nan_ind] = 0.
            
            # Now get all the lsfs
            lsf_array = []
            for j, xs in enumerate(xsm[0]):
                lsf = chunk_lsfs[xs]
                xneg = np.where(lsf < 0.0)
                if len(xneg[0]) > 0:
                    lsf[xneg[0]] = 0.0
                xnan = np.argwhere(np.isnan(lsf))
                if len(xnan) > 0:
                    print('Chunk {}, lsf {}: nan'.format(i, xs))
                else:
                    # No extra shifting to center here
                    lsf_array.append(lsf * wt[j])
            
            lsf_array = np.array(lsf_array)        
            lsf_av = np.median(lsf_array, axis=0)
            lsf_av = lsf_av / np.sum(lsf_av) # we use a different normalization -> change that?!
            lsfs_smoothed.append(lsf_av)\
    
    lsfs_smoothed = np.array(lsfs_smoothed)
    
    return lsfs_smoothed


def smooth_parameters_over_orders(parameters, par_name, chunks, deg=2):
    """For a given parameter (par_name) in a list of :class:`ParameterSet`
    objects, fit a polynomial of given degree over the central chunk pixels 
    within each order and return the evaluated results.
    
    :param parameters: A list of :class:`ParameterSet` objects.
    :type parameters: list
    :param par_name: Parameter key to smooth.
    :type par_name: str
    :param chunks: The chunks of the observation.
    :type chunks: :class:`ChunkArray`
    :param deg: Degree of the polynomial (default: 2).
    :type deg: int
    
    :return: A flattened array with the fit results of all chunks.
    :rtype: ndarray[nr_chunks]
    """
    
    pfits = np.zeros((len(chunks)))
    # Loop over orders that contain chunks
    for o in chunks.orders:
        chunk_ind = chunks.get_order_indices(o)
        # It only makes sense if there are at least two chunks within that order
        if len(chunk_ind) > 1:
            par_data = np.array([parameters[i][par_name] for i in chunk_ind])
            ch_pix = np.array([chunks[i].abspix[0] + (len(chunks[i]) // 2) \
                               for i in chunk_ind])
            
            # fit_polynomial returns fitted y-values and coefficients - only use the first
            pfits[chunk_ind] = fit_polynomial(ch_pix, par_data, deg=deg)[0]
        # If only one chunk within the order, use the original value
        else:
            pfits[chunk_ind[0]] = parameters[chunk_ind[0]][par_name]
        
    return pfits


def smooth_fitresult_over_orders(fit_result, par_name, deg=2):
    """For a given parameter (par_name) in a fit_result, fit a polynomial of 
    given degree over the central chunk pixels within each order and 
    return the evaluated results.
    
    :param fit_result: A list of :class:`LmfitResult` objects from the fit.
    :type fit_result: list
    :param par_name: Parameter key to smooth.
    :type par_name: str
    :param deg: Degree of the polynomial (default: 2).
    :type deg: int
    
    :return: An array with the smoothed fit results of all chunks.
    :rtype: ndarray[nr_chunks]
    """
    
    pfits = np.zeros((len(fit_result)))
    # Loop over orders that contain chunks
    chunk_orders = np.unique([res.chunk.order for res in fit_result])
    for o in chunk_orders:
        par_result_o = []
        ch_pix_o = []
        chunk_ind = []
        for i, res in enumerate(fit_result):
            if res.chunk.order == o:
                par_result_o.append(res.params[par_name])
                ch_pix_o.append(res.chunk.abspix[0] + (len(res.chunk) // 2))
                chunk_ind.append(i)
        
        par_result_o = np.array(par_result_o)
        ch_pix_o = np.array(ch_pix_o)
        
        # It only makes sense if there are at least two chunks within that order
        if len(par_result_o) > 1:
            # fit_polynomial returns fitted y-values and coefficients - only use the first
            pfits[chunk_ind] = fit_polynomial(ch_pix_o, par_result_o, deg=deg)[0]
        # If only one chunk within the order, use the original value
        else:
            pfits[chunk_ind[0]] = par_result_o[0]
    
    return pfits


def fit_polynomial(x_data, y_data, deg=2):
    """Fit a polynomial to data
    
    This routine masks out NaNs from the fit, and returns the input data
    if something goes wrong.
    
    :param x_data: The input x-data.
    :type x_data: ndarray, list
    :param y_data: The input y-data to be modelled.
    :type y_data: ndarray, list
    :param deg: Degree of the polynomial (default: 2).
    :type deg: int
    
    :return: The fitted polynomial data, or the input y-data if something went 
        wrong.
    :rtype: ndarray
    """
    try:
        if isinstance(x_data, list):
            x_data = np.array(x_data)
        elif isinstance(y_data, list):
            y_data = np.array(y_data)
        
        idx = np.where(np.isfinite(y_data))
        n = len(x_data[idx])
        
        pfit, stats = np.polynomial.Polynomial.fit(
                x_data[idx], y_data[idx], deg, full=True, window=(0, n), domain=(0, n))
        return pfit(x_data), pfit.coef
    except Exception as e:
        logging.error('Polynomial fitting failed. Falling back to input values', 
                      exc_info=True)
        return y_data, None


# Weights/error calculation
def analytic_chunk_weights(chunk_flux, wave_intercept, wave_slope):
    """Calculate analytic weight of a chunk, based on its spectral content
    
    Implemented with great parallels to the dop-code routine dsst_m.pro 
    (D. Fisher, Yale University), based on equations in Bulter&Marcy (1996): 
    sigma(mean) = 1./sqrt(sum(1/sig^2)).
    
    :param chunk_flux: Flux values of the pixels in a chunk.
    :type chunk_flux: ndarray[nr_pix]
    :param wave_intercept: The wavelength intercept w0 of that chunk.
    :type wave_intercept: float
    :param wave_slope: The wavelength dispersion w1 of that chunk.
    :type wave_slope: float
    
    :return: The weight of the chunk, representing its RV content.
    :rtype: float
    """
    eps = np.sqrt(chunk_flux[:-1])
    # slope: dI/d(pix)
    slope_int_per_pix = chunk_flux[1:] - chunk_flux[0:-1]
    # slope in real intensity per m/s
    slope_int_per_ms  = slope_int_per_pix * (wave_intercept / (_c * wave_slope))
    # Eqn. 5 in error write up
    weight = np.sum((slope_int_per_ms / eps)**2)
    
    return weight

