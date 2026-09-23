import h5py
import numpy as np
import logging
import dill
import os

from ..lib import h5quick
from ..components import ChunkArray, Chunk, SummedObservation #, Star
from ..template.base import StellarTemplate_Chunked
from ..models import lsf, wave, cont, spectrum
from ..fitters.lmfit_wrapper import LmfitWrapper

_group_keys = ('observation', 'chunks', 'params', 'errors', 'model')
_array_keys = ('reports', 'redchi2', 'residuals', 'medcnts')

_fileformats = {
        'h5py': ('.h5',),
        'dill': ('.pkl',)
        }


def create_results_dict(fit_results):
    """Pack most important fit results of an observation into a dictionary
    
    :param fit_results: The results of an observation.
    :type fit_results: list[:class:`LmfitWrapper.LmfitResult`]
    
    :return: The dictionary with the results.
    :rtype: dict
    """
    
    res_dict = {k: None for k in _group_keys}
    res_dict.update( {k: None for k in _array_keys} )
    
    # Collect observation info
    # Note: Unicode strings are currently not supported by h5py, so we need to
    # convert to bytestring (ascii). Special characters are replaced with
    # question marks
    obs = fit_results[0].chunk.observation
    res_dict['observation'] = {
            'instrument_name': obs.instrument.name.encode('utf8', 'replace'),
            'instrument_long': obs.instrument.longitude,
            'instrument_lat': obs.instrument.latitude,
            'instrument_alt': obs.instrument.altitude,
            'star_name': obs.star.name.encode('utf8', 'replace'),
            'orig_header': obs.orig_header.tostring(sep='\n').encode('utf8', 'replace'),
            'time_start': obs.time_start.isot.encode('utf8', 'replace'),
            'bary_date': obs.bary_date,
            'bary_vel_corr': obs.bary_vel_corr
            }
    
    # Additional star information if available
    if obs.star.coordinates is not None:
        res_dict['observation']['star_ra']  = obs.star.coordinates.ra.deg
        res_dict['observation']['star_dec'] = obs.star.coordinates.dec.deg
    if obs.star.proper_motion[0] is not None:
        res_dict['observation']['star_pmra']  = obs.star.proper_motion[0]
        res_dict['observation']['star_pmdec'] = obs.star.proper_motion[1]
    
    # Collect modelling info, and the original filename(s) of the modelled 
    # observation(s)
    res_dict['model'] = {
            'lsf_model': fit_results[0].model.lsf_model.name().encode('utf8', 'replace'),
            'lsf_pars_dict': fit_results[0].model.lsf_model.pars_dict,
            'wave_model': fit_results[0].model.wave_model.name().encode('utf8', 'replace'),
            'cont_model': fit_results[0].model.cont_model.name().encode('utf8', 'replace'),
            'iodine_file': os.path.abspath(fit_results[0].model.iodine_atlas.orig_filename).encode('utf8', 'replace'),
            'osample_factor': fit_results[0].model.osample_factor,
            'lsf_conv_width': fit_results[0].model.conv_width
            }
    # If not a fit result from O-star modelling
    if fit_results[0].model.stellar_template is not None:
        # Include the template info and the original filename of the observation
        res_dict['model']['stellar_template'] = os.path.abspath(fit_results[0].model.stellar_template.orig_filename).encode('utf8', 'replace')
        res_dict['observation']['orig_filename'] = os.path.abspath(obs.orig_filename).encode('utf8', 'replace')
        res_dict['observation']['temp_velocity'] = fit_results[0].model.stellar_template.velocity_offset
    else:
        # Include the original filenames of all modelled O-star observations (if more than one)
        if hasattr(obs, 'all_filenames'):
            res_dict['observation']['orig_filename'] = [os.path.abspath(f).encode('utf8', 'replace') for f in obs.all_filenames]
        else:
            res_dict['observation']['orig_filename'] = os.path.abspath(obs.orig_filename).encode('utf8', 'replace')

    # Assume parameter names are the same in all elements of input array
    param_names = list(fit_results[0].params.keys())
    nchunk = len(fit_results)

    # Collect info from all chunks
    res_dict['reports']   = np.array([res.report for res in fit_results], dtype='S')
    res_dict['redchi2']   = np.array([res.redchi for res in fit_results])
    res_dict['residuals'] = np.array([res.rel_residuals_rms() for res in fit_results])
    res_dict['medcnts']   = np.array([res.medcnts for res in fit_results])
    
    res_dict['chunks'] = {k: np.zeros(nchunk, dtype='int') for k in ('order', 'firstpix', 'lastpix', 'padding')}
    res_dict['params'] = {k: np.zeros(nchunk, dtype='float64') for k in param_names}
    res_dict['errors'] = {k: np.zeros(nchunk, dtype='float64') for k in param_names}
    
    for i in range(nchunk):
        res = fit_results[i]
        # Get chunk info
        res_dict['chunks']['order'][i]    = res.chunk.order
        res_dict['chunks']['firstpix'][i] = res.chunk.abspix[0]
        res_dict['chunks']['lastpix'][i]  = res.chunk.abspix[-1]
        res_dict['chunks']['padding'][i]  = res.chunk.padding
        # Get parameter values and errors
        for p in param_names:
            res_dict['params'][p][i] = res.params[p]
            res_dict['errors'][p][i] = res.errors[p]
    
    return res_dict


def filetype_from_ext(filename):
    """Determine the filetype from the filename extension
    
    :param filename: The pathname of the file.
    :type filename: str
    
    :return: The filetype matching the filename extension.
    :rtype: str
    """
    
    # Split the filename and check the extension
    ext = os.path.splitext(filename)[1]
    
    file_format = None
    for key in _fileformats.keys():
        if ext in _fileformats[key]:
            file_format = key
    
    if isinstance(file_format, str):
        return file_format
    else:
        raise KeyError('File extension {} does not correspond to any known format!'.format(ext))
    
    
def check_filename_format(filename, filetype, correct=True):
    """Check whether the extension of the filename matches the chosen filetype
    
    If the filename does not match and keyword 'correct' is True, a corrected 
    filename is returned. Otherwise just the old one.
    
    :param filename: The chosen filename.
    :type filename: str
    :param filetype: The chosen filetype.
    :type filetype: str
    :param correct: Whether to correct the filename if it does not match the
        filetype. Defaults to True.
    :type correct: bool
    
    :return: Whether the filename matches the filetype or not.
    :rtype: bool
    :return: The filename (corrected if it does not match and 'correct==True').
    :rtype: str
    """
    # Check whether the filetype is known
    if filetype in _fileformats.keys():
        # Split the filename and check the extension
        file_ext = os.path.splitext(filename)
        if file_ext[1] not in _fileformats[filetype]:
            logging.warning('The extension {} does not match the filetype {}.'.format(
                    file_ext[1], filetype))
            logging.warning('It should be one of: {}'.format(_fileformats[filetype]))
            
            # Possibly correct
            if correct:
                logging.warning('Correcting it to: {}'.format(_fileformats[filetype][0]))
                new_filename = file_ext[0] + _fileformats[filetype][0]
            else:
                new_filename = filename
                
            return False, new_filename
        else:
            return True, filename
    
    else:
        raise KeyError('Filetype {} is not known! Must be either of: '.format(
                filetype), _fileformats.keys())
    

def save_results(filename, fit_results, filetype='h5py'):
    """Preliminary function to save a set of fit results
    
    The results are either saved in 'h5py' format, which is basically a
    dictionary of the most important fit results, or in 'dill' format, where
    the whole object structure is saved and can be recovered later for in-depth 
    analysis. Note that this requires a lot more memory!
    If 'filename' exists, it will be overwritten.
    
    :param filename: Output path for the results file.
    :type filename: str
    :param fit_results: The list of :class:`LmfitWrapper.LmfitResult` objects 
        to save.
    :type fit_results: list[:class:`LmfitWrapper.LmfitResult`]
    :param filetype: In which format should the results be written? Default is 
        'h5py', which saves the most important data in a dictionary format to a 
        compact hdf5-file. If 'dill' is specified instead, the full object 
        structure is saved.
    :type filetype: str
    """
    
    # First check whether the filename matches the chosen type, and correct if
    # it does not
    match, new_filename = check_filename_format(filename, filetype, correct=True)
    
    # If savetype is dill, save the whole object structure
    if filetype == 'dill':
        with open(new_filename, 'wb') as f:
            dill.dump(fit_results, f)
    
    # Else create a results dictionary and save that as hdf5
    elif filetype == 'h5py':
        res_dict = create_results_dict(fit_results)
        
        with h5py.File(new_filename, 'w') as h:
            for key in _group_keys:
                h5quick.dict_to_group(res_dict[key], h, key)
            
            for key in _array_keys:
                h[key] = res_dict[key]
    
    else:
        raise KeyError('The savetype must be either of "h5py" or "dill",' + \
                       'but {} was supplied!'.format(filetype))        
    

def load_results(filename, filetype='h5py', force=True):
    """Function to load a set of results
    
    Returns them either as dictionary, if the filetype is 'h5py', or as the
    original object structure if filetype is 'dill'. If 'filetype' does not 
    match the extension of the filename and keyword 'force' is True, then 
    the routine attempts to still load the results as corresponding to the
    filename extension.
    
    :param filename: Path to the file.
    :type filename: str
    :param filetype: If 'h5py', the results are returned as dictionary 
        (default). If 'dill', the original object structure is recovered.
    :type filetype: str
    :param force: Whether or not to force the loading of the results, even if
        the filename does not match the filetype.
    :type force: bool
    
    :return: The fit results.
    :rtype: dict or list[:class:`LmfitResult`]
    """
    
    # First check whether the filename matches the chosen type
    match, filename = check_filename_format(filename, filetype, correct=False)
    
    # If it does not match but loading is forced, attempt to determine the
    # format from the filename
    if not match and force:
        filetype = filetype_from_ext(filename)
        match = True
    
    if match:
        # For dill: recover the whole object structure
        if filetype == 'dill':
            try:
                with open(filename, 'rb') as f:
                    fit_results = dill.load(f)
                return fit_results
            except Exception as e:
                raise(e)
            
        # For h5py: load the results as dictionary
        elif filetype == 'h5py':
            try:
                fit_results = {}
                with h5py.File(filename, 'r') as h:
                    for key in _group_keys + _array_keys:
                        try:
                            fit_results[key] = h5quick.h5data(h[key])
                        except:
                            fit_results[key] = None
                            logging.warning('Key {} not in result file!'.format(key))
                return fit_results
            except Exception as e:
                raise(e)        
    
    else:
        raise KeyError('The filename {} does not match the fileformat {}!'.format(
                filename, filetype))
    

def restore_results_object(utilities, filename, temp_path=None, obs_path=None,
                           iod_path=None):
    """A wrapper function to restore saved results as 
    :class:`LmfitWrapper.LmfitResult` objects
    
    If the results are saved as '.pkl' (dill format), this is straight forward.
    If the save format however is '.h5' (HDF5/h5py), it's more complicated as
    the objects need to be created from scratch, using the information 
    contained in the saved dictionary. IMPORTANT: This only works if all 
    pathnames of the input files (observation spectra, I2 atlas, templates) 
    have not changed!
    UPGRADE: Now the template and/or observation and/or I2 atlas pathnames
    can be handed as optional arguments, in case they DO have changed.
    
    :param utilities: The utilities module for the instrument used.
    :type utilities: library
    :param filename: The filename of the saved results.
    :type filename: str
    :param temp_path: An alternative pathname to the stellar template (optional,
        in case the pathname has changed).
    :type temp_path: str, or None
    :param obs_path: An alternative pathname to the observation (optional,
        in case the pathname has changed). If more than one observation,
        this should be a list of pathnames.
    :type obs_path: str, list, or None
    :param iod_path: An alternative pathname to the I2 atlas (optional,
        in case the pathname has changed).
    :type iod_path: str, or None
    
    :return: The chunks of the modelled observation.
    :rtype: :class:`ChunkArray`
    :return: The list of fit results of the modelled observation.
    :rtype: list
    """
    
    # First check the filetype of the filename
    filetype = filetype_from_ext(filename)
    results = load_results(filename, filetype)
    
    # If it was a 'dill' file, the result should already be recovered as
    # object structure and we only need to create the chunk array
    if not isinstance(results, dict):
        obs_chunks = ChunkArray()
        for r in results:
            obs_chunks.append(r.chunk)
        
        return obs_chunks, results
    
    # For a 'h5py' file in contrast, we need to try and manually rebuild the
    # object structure from the information in the saved dictionary
    else:
        
        # First get the names of the observation(s) (this could be more than
        # one, in case several observation spectra were summed up)
        # (first check if pathnames(s) were handed as argument)
        if isinstance(obs_path, (list, tuple, str)):
            if isinstance(obs_path, str):
                obs_path = [obs_path]
            elif isinstance(obs_path, tuple):
                obs_path = list(obs_path)
        else:
            if isinstance(results['observation']['orig_filename'], (list,np.ndarray)):
                obs_path = [f.decode() for f in results['observation']['orig_filename']]
            else:
                obs_path = [results['observation']['orig_filename'].decode()]
        
        # Now the stellar template name (if any - for hot star modelling during
        # template creation this is None)
        # (only if pathname was NOT handed as argument)
        if not isinstance(temp_path, str):
            if 'stellar_template' in results['model'].keys():
                temp_path = results['model']['stellar_template'].decode()
            else:
                temp_path = None
        
        # Now other important information: The pathname to the I2 atlas, the
        # used oversampling factor, the LSF convolution width, the name of the 
        # used LSF, and info about the chunks and best-fit parameters
        # (I2 atlas: check if it was handed as argument)
        if not isinstance(iod_path, str):
            iod_path   = results['model']['iodine_file'].decode()
        osample        = results['model']['osample_factor']
        lsf_conv_width = results['model']['lsf_conv_width']
        lsf_name       = results['model']['lsf_model'].decode()
        lsf_pars_dict  = results['model']['lsf_pars_dict']
        # For backward compatability: Check
        if 'wave_model' in results['model'].keys():
            wave_name      = results['model']['wave_model'].decode()
            cont_name      = results['model']['cont_model'].decode()
        else:
            wave_name = 'LinearWaveModel'
            cont_name = 'LinearContinuumModel'
        res_chunks     = results['chunks']
        res_params     = results['params']
        res_errors     = results['errors']
        
        # The orders covered by the chunks, the width of the chunks, and the
        # padding of the chunks (these are constructed from implicit 
        # information in the res_chunks dictionary)
        orders  = np.unique(np.array([o for o in res_chunks['order']]))
        
        # Load the observation data (either single one or multiple summed up)
        all_obs = [utilities.load_pyodine.ObservationWrapper(f) for f in obs_path]
        if len(obs_path) > 1:
            obs = SummedObservation(*all_obs)
        else:
            obs = all_obs[0]
        
        # If a template was used, load it also
        # ToDo: Allow also for 'StellarTemplate' object?
        if temp_path:
            temp = StellarTemplate_Chunked(temp_path)
        else:
            temp = None
        
        # Load the I2 atlas
        iod = utilities.load_pyodine.IodineTemplate(iod_path)
        
        # Build the model and fitter
        lsf_model  = lsf.model_index[lsf_name]
        # Adapt the LSF setup to the instrument
        lsf_model.adapt_LSF(lsf_pars_dict)
        
        wave_model = wave.model_index[wave_name]
        cont_model = cont.model_index[cont_name]
        model      = spectrum.SimpleModel(
                lsf_model, wave_model, cont_model, iod, stellar_template=temp, 
                osample_factor=osample, conv_width=lsf_conv_width)
        
        # Initialize the fitter object
        fitter = LmfitWrapper(model)
        """
        # Build the chunk array: If a template is given, use the wave_defined
        # algorithm, otherwise the user_defined function (default functions
        # in the standard pyodine distribution)
        if temp:
            # Compute possible order shifts between template and observation,
            # by searching for the best coverage of first template order in
            # observation
            obs_order_min, min_coverage = obs.check_wavelength_range(
                    temp[0].w0, temp[len(temp.get_order_indices(temp.orders_unique[0]))-1].w0)
            order_correction = obs_order_min - temp.orders_unique[0]
            logging.info('Order correction: {}\n'.format(order_correction))
            
            obs_chunks = chunks.wave_defined(obs, temp, width=width, orders=orders-order_correction, 
                                             padding=padding, order_correction=order_correction)
        else:
            chunks_per_order = len(np.where(res_chunks['order'] == orders[0])[0])
            pix_offset0      = res_chunks['firstpix'][0]
            
            obs_chunks = chunks.user_defined(obs, width=width, orders=orders, 
                                             padding=padding, chunks_per_order=chunks_per_order, 
                                             pix_offset0=pix_offset0)
        """
        obs_chunks = build_chunk_array(obs, res_chunks)
        
        # The total number of chunks, and the number of chunks in order 0
        nr_chunks_total  = len(obs_chunks)
        nr_chunks_order0 = len(obs_chunks.get_order(orders[0]))
        
        logging.info('Total number of created chunks: {} (in result file: {})'.format(
                nr_chunks_total, len(res_chunks['order'])))
        logging.info('Number of created chunks in order 0: {}'.format(nr_chunks_order0))
        
        # Loop over the chunks to build the fit_result object
        fit_results = []
        for i, chunk in enumerate(obs_chunks):
            # First guess the chunk parameters to create the ParameterSet object
            pars = model.guess_params(chunk)
            
            # Then fill it with the best-fit results for that chunk
            for key in res_params.keys():
                pars[key] = res_params[key][i]
            
            # Convert it to an LM-fit parameter object, fix all the parameters
            lmfit_pars = fitter.convert_params(pars, to_lmfit=True)
            for key in lmfit_pars.keys():
                lmfit_pars[key].set(vary=False)
            
            # Finally 'fit' it (work-around to create a fully functional 
            # FitResults object)
            fit_results.append(fitter.fit(chunk, lmfit_pars, chunk_ind=i))
            
            # And add the correct errors from the dictionary
            for key in res_errors.keys():
                if fit_results[-1].lmfit_result:
                    fit_results[-1].lmfit_result.params[key].stderr = res_errors[key][i]
        
        return obs_chunks, fit_results


def build_chunk_array(obs, chunk_dict):
    """Build a chunk array, using the order, pixel and padding information
    contained in a dictionary that was loaded from a saved results file
    
    :param obs: The observation to chunk.
    :type obs: :class:`ObservationWrapper` or :class:`SummedObservation`
    :param chunk_dict: A dictionary with orders ('order'), pixels ('firstpix'
        and 'lastpix'), and padding ('padding') information for each chunk.
    :type chunk_dict: dict
    
    :return: The built chunk array.
    :rtype: :class:`ChunkArray`
    """
    
    chunks = ChunkArray()
    
    for i in range(len(chunk_dict['order'])):
        
        order   = chunk_dict['order'][i]
        pixels  = np.arange(chunk_dict['firstpix'][i], chunk_dict['lastpix'][i]+1, dtype='int')
        padding = chunk_dict['padding'][i]
        
        chunk = Chunk(obs, order, pixels, padding)
        chunks.append(chunk)
    
    return chunks
    
