import os
import h5py
import numpy as np
import logging
import sys
from progressbar import ProgressBar

from .. import fitters
from ..lib import h5quick
from ..lib.misc import return_existing_files
from .combine_vels import combine_chunk_velocities, combine_chunk_velocities_dop
from .bary_vel_corr import bvc_wrapper

# Speed of light (m/s)
_c = 299792458.


class CombinedResults():
    """Container for timeseries fitting results
    
    This object class is a container for all the individual fit results from
    all observations of a star, and can be used as input in the final velocity 
    weighting to receive RVs. If that has already been performed, it can 
    additionally store the results from the weighting algorithm and RVs.
    
    :param filename: A string to the path from which to load an existing 
        :class:`CombinedResults` object. If None, then the object is
        initialized without data.
    :type filename: str, or None
    """
    
    def __init__(self, filename=None):
        
        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                                format='%(message)s')
        
        # Check whether one combined result should be loaded
        if isinstance(filename, str):
            try:
                self.filename = filename
                self.load_combined(self.filename)
            except Exception as e:
                logging.error('Problem loading combined results:', exc_info=True)
    
    
    def compute_bvcs(self, use_hip=True, bary_dict=None, precise=True, 
                     temp_vel=None, ref_vel=None, solar=False):
        """Compute barycentric velocities, using a wrapper function 
        for the barycorrpy package
        
        The keyword precise determines whether to use the multiplicative 
        algorithm, which makes use of the absolute measured Doppler-shift, or
        the less precise additive "predictive" algorithm. If True, make sure 
        that the algorithm has all the information needed to compute the 
        absolute Doppler-shifts (i.e. template velocity and/or velocity of 
        reference spectrum).
        
        :param use_hip: Whether to use the built-in hip catalogue to find 
            the star's coordinates; otherwise fall back to coordinates from the
            timeseries dictionary. Defaults to True.
        :type use_hip: bool
        :param use_hip: Whether to use the barycentric corrected JD (BJD) from
            barycorrpy as new date. Defaults to True.
        :type use_hip: bool
        :param bary_dict: A dictionary with stellar (and possible observatory) 
            information that should be used in the computation of barycentric 
            velocity corrections. If None, the info from the model results is 
            used.
        :type bary_dict: dict, or None
        :param precise: Whether or not to use the very precise multiplicative
            algorithm.
        :type precise: bool
        :param temp_vel: If precise=True, you may hand a template velocity 
            offset (in m/s) here to compute the absolute Doppler shifts. If 
            None is given, the algorithm will try and use the template velocity 
            of the results file(s).
        :type temp_vel: float, int, or None
        :param ref_vel: If precise=True, you may hand a velocity offset of a 
            reference spectrum (in m/s) here, which was used to compute the 
            template velocity offset. If None is given, the algorithm will 
            assume it to be 0.
        :type ref_vel: float, int, or None
        :param solar: If True, return the barycentric correction for the Sun as 
            target. Defaults to False. (But even then, if the 'star_name' in 
            the bvc_dict is 'Sun', the solar BC correction will be done.)
        :type solar: bool
        """
        
        bvc_dict = create_bvc_dict(self.info, self.timeseries, 
                                   bary_dict=bary_dict)
        
        # If not precise
        if precise == False:
            # Compute predictive bvcs
            bvcs, bjds = bvc_wrapper(bvc_dict, self.timeseries, 
                                     use_hip=use_hip, solar=solar)
            
            # Correct the RVs of the observations (if they exist already)
            if 'rv' in self.timeseries.keys():
                self.timeseries['rv_bc'] = self.timeseries['rv'] + bvcs
                
            # Update the timeseries dict
            self.timeseries['bary_vel_corr'] = bvcs
            self.timeseries['bary_date_corr'] = bjds
        
        # If precise
        else:
            # Compute the absolute Doppler shifts for each observation
            if not isinstance(temp_vel, (float, int)):
                temp_vel = self.info['temp_velocity']
            if not isinstance(ref_vel, (float, int)):
                ref_vel = 0.
            
            rv_absolute = np.array(self.timeseries['rv']) + temp_vel + ref_vel
            z_absolute  = rv_absolute / _c
            
            #print(rv_absolute)
            
            # Correct the RVs of the observations (this is done directly by
            # barycorrpy in this case)
            rv_corrected, bjds = bvc_wrapper(bvc_dict, self.timeseries,
                                             use_hip=use_hip, z_meas=z_absolute, 
                                             solar=solar)
            
            # Update the timeseries dict
            self.timeseries['rv_bc'] = rv_corrected
            self.timeseries['bary_vel_corr'] = rv_corrected - rv_absolute
            self.timeseries['bary_date_corr'] = bjds
        
        self.fill_timeseries_attributes()
        
        # And write the bvc_dict items into the info dictionary
        for key in ['star_ra', 'star_dec', 'star_pmra', 'star_pmdec', 'star_rv0', 'star_name']:
            self.info[key] = bvc_dict[key]
    
    
    def create_timeseries(self, weighting_pars=None, do_crx=True, crx_pars=None):
        """Create the timeseries data (weighted and unweighted RVs with 
        uncertainties, chunk-to-chunk scatter, RV precision measures, and
        optionally chromatic indices with uncertainties)
        
        :param weighting_pars: A dictionary of weighting parameter values 
            needed in the weighting algorithm. If None, a dictionary of 
            default values is used there.
        :type weighting_pars: dict, or None
        :param do_crx: Whether to also compute chromatic indices of the 
            observations. Defaults to True.
        :type do_crx: bool
        """
        
        velocities = self.params['velocity']
        #bvc = self.timeseries['bary_vel_corr']
        wavelengths = None
        if do_crx:
            wavelengths = self.params['wave_intercept']
        tseries, self.auxiliary, self.weighting_pars = combine_chunk_velocities(
                velocities, self.nr_orders, wavelengths=wavelengths, 
                weighting_pars=weighting_pars, crx_pars=crx_pars)
        
        self.timeseries.update(tseries)
        self.fill_timeseries_attributes()
    
    
    def create_timeseries_dop(self, weighting_pars=None, do_crx=True, crx_pars=None):
        """Create the timeseries data, using the algorithm as deployed in the
        dop code (weighted and unweighted RVs with 
        uncertainties, chunk-to-chunk scatter, RV precision measures, and
        optionally chromatic indices with uncertainties)
        
        :param weighting_pars: A dictionary of weighting parameter values 
            needed in the weighting algorithm. If None, a dictionary of 
            default values is used there.
        :type weighting_pars: dict, or None
        :param do_crx: Whether to also compute chromatic indices of the 
            observations. Defaults to True.
        :type do_crx: bool
        """
        
        velocities = self.params['velocity']
        bvc = self.timeseries['bary_vel_corr']
        wavelengths = None
        if do_crx:
            wavelengths = self.params['wave_intercept']
        tseries, self.auxiliary, self.weighting_pars = combine_chunk_velocities_dop(
                velocities, self.redchi2, self.medcnts, wavelengths=wavelengths, 
                weighting_pars=weighting_pars, crx_pars=crx_pars)
        
        self.timeseries.update(tseries)
        self.fill_timeseries_attributes()
        
    
    
    def fill_timeseries_attributes(self):
        """Create an object attribute for each entry in the the self.timeseries 
        dictionary (to make the results easier accessible).
        """
        for key, value in self.timeseries.items():
            setattr(self, key, value)
    
    
    def results_to_txt(self, filename, outkeys=None, delimiter='\t', 
                       header='', outformat=None, detailed=False, flux_chunk=None):
        """Write timeseries results to a txt-file
        
        :param filename: The output filepath.
        :type filename: str
        :param outkeys: Which of the self.timeseries items to write to file. If
            None, write the 'bary_date', 'rv' and 'rv_err' entries by default.
        :type outkeys: str, list, tuple, or None
        :param delimiter: The delimiter used in the txt-file. Defaults to '\t'.
        :type delimiter: str
        :param header: Potential header row to write before the data (e.g. the
            keys). If empty, no header row will be written.
        :type header: str
        :param outformat: The output format of each column. Make sure that this
            matches the data types (particularly for strings)!
        :type outformat: str, list, or None
        :param detailed: If True, ignore all other settings of keys etc. and 
            write a detailed results file with star info, filenames, dates 
            (uncorrected and corrected), RVs (uncorrected and corrected), chunk
            scatters, RV uncertainties, and flux measures.
        :type detailed: bool
        :param flux_chunk: A chunk index (or list of indices) which to use for
            the median flux estimate (if detailed=True). If None, then the 
            chunk closest to 5500 Angstrom is used (default).
        :type flux_chunk: list, tuple, int, or None
        """
        
        logging.info('Writing results to txt file: {}'.format(filename))
        
        if detailed:
            outkeys = ['orig_filename', 'bary_date', 'bary_date_corr',
                       'rv', 'rv_bc', 'c2c_scatter', 'rv_err', 'medcnts']
            outformat = ['%50s', '%20.8f', '%20.8f',
                         '%14.4f', '%14.4f', '%14.4f', '%14.4f', '%14.1f']
            delimiter = '\t'
            header = ' Instrument: {}\n'.format(self.info['instrument_name']) + \
                     ' Star = {}\n'.format(self.info['star_name']) + \
                     ' RA (ICRS,deg) = {}\n'.format(self.info['star_ra']) + \
                     ' DEC (ICRS,deg) = {}\n'.format(self.info['star_dec']) + \
                     ' pmRA (mas/yr) = {}\n'.format(self.info['star_pmra']) + \
                     ' pmDEC (mas/yr) = {}\n'.format(self.info['star_pmdec']) + \
                     ' parallax (mas) = {}\n'.format('?') + \
                     ' Vrad (km/s) = {}\n'.format(self.info['star_rv0']) + \
                     ' Vtemplate (km/s) = {}\n'.format(self.info['temp_velocity']) + \
                     '-' * 130 + '\n' + \
                     '    Filename          JD_MID           BJD_MID ' + \
                     '         RV_raw        RV_bc        c2c_scatter        RV_err          Flux\n' + \
                     '                      [days]           [days]  ' + \
                     '         [m/s]         [m/s]           [m/s]            [m/s]          [ADU]\n' + \
                     '-' * 130
            logging.debug('Detailed output.')
        else:
            if not isinstance(outkeys, (str,list,tuple)):
                outkeys = ['bary_date', 'rv', 'rv_err']
            elif isinstance(outkeys, str):
                outkeys = [outkeys]
            
            logging.debug('Keys:')
            logging.debug(outkeys)
        
        out_data = []
        data_types = []
        for key in outkeys:
            if key in self.timeseries.keys():
                out_data.append(self.timeseries[key])
                if isinstance(self.timeseries[key][0], str):
                    data_types += ['U100']
                else:
                    data_types += [type(self.timeseries[key][0])]
        
        if detailed:
            # Omit the path of the filenames (just the basename)
            for i in range(len(out_data[0])):
                out_data[0][i] = os.path.basename(out_data[0][i])
            # Include the median flux around wavelength 5500 Angstrom (roughly)
            # (or at the desired chunk indices)
            if not isinstance(flux_chunk, (int, list, tuple)):
                flux_chunk = np.argmin(np.abs(5500 - np.median(self.params['wave_intercept'], axis=0)))
            if isinstance(flux_chunk, int):
                out_data.append(self.medcnts[:,flux_chunk])
            else:
                out_data.append(np.median(self.medcnts[:,flux_chunk], axis=1))
            data_types += [type(out_data[-1][0])]
        
        out_array = np.zeros(len(out_data[0]), 
                             dtype=[('v{}'.format(i), data_types[i]) for i in range(len(data_types))])
        
        for i in range(len(data_types)):
            out_array['v{}'.format(i)] = out_data[i]
        
        np.savetxt(filename, out_array.T, delimiter=delimiter, header=header,
                   fmt=outformat)
    
    
    def load_individual_results(self, filenames, compact=False):
        """Load individual fit results
        
        :param filenames: The pathnames of the files to the individually saved
            results to load.
        :type filenames: list or tuple 
        :param compact: If True, use a compact version of the CombinedResults
            (only when loading individual results), where only the bare minimum 
            of parameters is loaded (to prevent memory crashs for very large 
            time series). Defaults to False.
        :type compact: bool
        """
        if compact:
            logging.info('Using compact version of CombinedResults.')
        # First of all, sort the input list
        filenames.sort()
        
        # First check whether all filenames supplied actually exist,
        # and only use the ones that do
        filenames, bad_files = return_existing_files(filenames)
        if len(bad_files) > 0:
            logging.info('Non-existing files:')
            for f in bad_files:
                logging.info(f)
        
        self.nr_files = len(filenames)        
        
        # Get param names and general info from first file. It can be either
        # 'h5py' or 'dill', so take care of that.
        filetype = fitters.results_io.filetype_from_ext(filenames[0])
        result = fitters.results_io.load_results(
                filenames[0], filetype=filetype)
                        
        # If it was a 'dill' file, transform the recovered object 
        # structure to a dictionary
        if not isinstance(result, dict):
            result = fitters.results_io.create_results_dict(result)
        
        # For the compact CombinedResults, only use absolutely necessary model
        # parameters (wavelength, velocity) and for chunks only their order nr
        if compact:
            self.param_names = [k for k in result['params'].keys() if 'wave' in k or
                                'velocity' in k]
            self.chunk_names = [k for k in result['chunks'].keys() if 'order' in k]
        else:
            self.param_names = [k for k in result['params'].keys()]
            self.chunk_names = [k for k in result['chunks'].keys()]
        
        # General info
        self.info = {
                'star_name': result['observation']['star_name'].decode(),
                'instrument_name': result['observation']['instrument_name'].decode(),
                'instrument_long': result['observation']['instrument_long'],
                'instrument_lat': result['observation']['instrument_lat'],
                'instrument_alt': result['observation']['instrument_alt']
                }
        
        # This was added only later, so make it backwards compatible
        if 'temp_velocity' in result['observation']:
            self.info['temp_velocity'] = result['observation']['temp_velocity']
        
        # Model info
        if 'model' in result.keys() and result['model'] != None:
            self.info['lsf_model'] = result['model']['lsf_model'].decode()
            self.info['stellar_template'] = result['model']['stellar_template'].decode()
            self.info['iodine_file'] = result['model']['iodine_file'].decode()
            self.info['osample_factor'] = result['model']['osample_factor']
            self.info['lsf_conv_width'] = result['model']['lsf_conv_width']
        
        self.nr_chunks = len(result['chunks'][self.chunk_names[0]])
        self.orders = np.unique(result['chunks']['order'])
        self.nr_orders = len(self.orders)
        
        # Allocate arrays
        self.timeseries = {
            'bary_date': np.zeros(self.nr_files),
            'bary_vel_corr': np.zeros(self.nr_files),
            'orig_filename': [''] * self.nr_files,
            'star_ra': np.zeros(self.nr_files),
            'star_dec': np.zeros(self.nr_files),
            'star_pmra': np.zeros(self.nr_files),
            'star_pmdec': np.zeros(self.nr_files)
        }
        
        self.params = {k: np.zeros((self.nr_files, self.nr_chunks)) for k in self.param_names}
        self.errors = {k: np.zeros((self.nr_files, self.nr_chunks)) for k in self.param_names}
        self.chunks = {k: np.zeros((self.nr_files, self.nr_chunks)) for k in self.chunk_names}
        self.redchi2 = np.zeros((self.nr_files, self.nr_chunks))
        self.residuals = np.zeros((self.nr_files, self.nr_chunks))
        self.medcnts = np.zeros((self.nr_files, self.nr_chunks))
        
        # Set up progressbar
        bar = ProgressBar(max_value=self.nr_files, redirect_stdout=True)
        bar.update(0)
        # Now load the results from all files and fill up the object properties,
        # again making sure about the file formats
        for i, file in enumerate(filenames):
            filetype = fitters.results_io.filetype_from_ext(file)
            result = fitters.results_io.load_results(file, filetype=filetype)
            # If it was a 'dill' file, transform the recovered object structure
            # to a dictionary
            if not isinstance(result, dict):
                result = fitters.results_io.create_results_dict(result)
            
            for k in self.timeseries.keys():
                self.timeseries[k][i] = result['observation'][k] if k in result['observation'] else np.nan
                if k == 'orig_filename':
                    self.timeseries[k][i] = self.timeseries[k][i].decode()
            for k in self.param_names:
                self.params[k][i] = result['params'][k]
                self.errors[k][i] = result['errors'][k]
            for k in self.chunk_names:
                self.chunks[k][i] = result['chunks'][k]
            self.redchi2[i] = result['redchi2']
            self.residuals[i] = result['residuals']
            self.medcnts[i] = result['medcnts']
            # Update progressbar
            bar.update(i+1)
        
        self.timeseries['res_filename'] = [os.path.abspath(f) for f in filenames]
        
        self.fill_timeseries_attributes()
        
        # Finish progressbar
        bar.finish()
        
        # Initiate an empty auxiliary and weighting pars attribute
        self.auxiliary = {}
        self.weighting_pars = {}
    
        
    def save_combined(self, filename):
        """Save the combined fit results to file
        
        :param filename: Save under this filename.
        :type filename: str
        """
        
        logging.info('Saving combined results to {}'.format(filename))
        
        # Make sure that the file extension matches the h5py format, and
        # correct if this is not the case
        match, new_filename = fitters.results_io.check_filename_format(
                filename, 'h5py', correct=True)
        
        with h5py.File(new_filename, 'w') as h:
            for k in ('res_filename', 'orig_filename'):
                self.timeseries[k] = [f.encode('utf8', 'replace') for f in self.timeseries[k]]
            h5quick.dict_to_group(self.timeseries, h, 'timeseries')
            h5quick.dict_to_group(self.auxiliary, h, 'auxiliary')
            h5quick.dict_to_group(self.params, h, 'params')
            h5quick.dict_to_group(self.errors, h, 'errors')
            h5quick.dict_to_group(self.chunks, h, 'chunks')
            for k in self.info:
                if isinstance(self.info, str):
                    self.info[k] = self.info[k].encode('utf8', 'replace')
            h5quick.dict_to_group(self.info, h, 'info')
            h5quick.dict_to_group(self.weighting_pars, h, 'weighting_pars')
            h['redchi2'] = self.redchi2
            h['residuals'] = self.residuals
            h['medcnts'] = self.medcnts
            #h['res_filenames'] = [f.encode('utf8', 'replace') for f in self.res_filenames]
    
    
    def load_combined(self, filename):
        """Load a combined fit results object from file
        
        :param filename: The pathname of the file.
        :type filename: str
        """
        
        logging.info('Loading combined results from {}'.format(filename))
        
        with h5py.File(filename, 'r') as h:
            self.timeseries = h5quick.h5data(h['timeseries'])
            for k in ('res_filename', 'orig_filename'):
                self.timeseries[k] = [f.decode() for f in self.timeseries[k]]
            self.auxiliary = h5quick.h5data(h['auxiliary'])
            self.params = h5quick.h5data(h['params'])
            self.errors = h5quick.h5data(h['errors'])
            self.chunks = h5quick.h5data(h['chunks'])
            self.info = h5quick.h5data(h['info'])
            for k in self.info:
                if isinstance(self.info[k], np.bytes_):
                    self.info[k] = self.info[k].decode()
            
            try:
                self.weighting_pars = h5quick.h5data(h['weighting_pars'])
            except:
                self.weighting_pars = {}
            self.redchi2 = h5quick.h5data(h['redchi2'])
            self.residuals = h5quick.h5data(h['residuals'])
            self.medcnts = h5quick.h5data(h['medcnts'])
            #self.res_filenames = [f.decode() for f in h5quick.h5data(h['res_filenames'])]
            
        self.nr_chunks = self.chunks['order'].shape[1]
        self.orders = np.unique(self.chunks['order'][0])
        self.nr_orders = len(self.orders)
        
        self.param_names = [k for k in self.params.keys()]
        self.chunk_names = [k for k in self.chunks.keys()]
        
        self.fill_timeseries_attributes()
        self.nr_files = len(self.res_filename)
        
        
    def remove_observations(self, res_names=None, obs_names=None):
        """Remove a number of individual results from the object, either
        by their individual result filenames or their original observation
        filenames.
        
        :param res_names: A list of individual result filenames to remove. If 
            None, supply 'obs_names' instead.
        :type res_names: list or tuple, or None
        :param obs_names: A list of original observation filenames to remove. If 
            None, supply 'res_names' instead.
        :type obs_names: list or tuple, or None
        """
        
        if isinstance(res_names, (list,tuple)):
            res_names.sort()
            inds = self._return_indices_of_filenames(res_names, self.timeseries['res_filename'])
        elif isinstance(obs_names, (list,tuple)):
            obs_names.sort()
            inds = self._return_indices_of_filenames(obs_names, self.timeseries['orig_filename'])
        else:
            raise KeyError('Either of "res_names" or "obs_names" must be list or tuple!')
        
        # Remove from timeseries
        for key in self.timeseries.keys():
            if isinstance(self.timeseries[key], np.ndarray):
                self.timeseries[key] = np.delete(self.timeseries[key], inds, axis=0)
            elif isinstance(self.timeseries[key], list):
                #self.timeseries[key] = [self.timeseries[key][i] for i in inds]
                for i in sorted(inds, reverse=True):
                    del self.timeseries[key][i]
        
        # Remove from auxiliary
        for key in self.auxiliary.keys():
            if isinstance(self.auxiliary[key], np.ndarray):
                self.auxiliary[key] = np.delete(self.auxiliary[key], inds, axis=0)
            elif isinstance(self.auxiliary[key], list):
                #self.auxiliary[key] = [self.auxiliary[key][i] for i in inds]
                for i in sorted(inds, reverse=True):
                    del self.auxiliary[key][i]
        
        # Remove from params and errors
        for key in self.params.keys():
            self.params[key] = np.delete(self.params[key], inds, axis=0)
            self.errors[key] = np.delete(self.errors[key], inds, axis=0)
        
        # Remove from chunks
        for key in self.chunks.keys():
            self.chunks[key] = np.delete(self.chunks[key], inds, axis=0)
        
        # Remove from redchi2, residuals and medcnts
        self.redchi2 = np.delete(self.redchi2, inds, axis=0)
        self.residuals = np.delete(self.residuals, inds, axis=0)
        self.medcnts = np.delete(self.medcnts, inds, axis=0)
        
        # Adapt the nr_files, and finally the timeseries attributes
        self.nr_files -= len(inds)
        self.fill_timeseries_attributes()        
    
    
    def _return_indices_of_filenames(self, filenames, all_names):
        """Return indices of filenames within all_names (if they are in there)
        
        :param filenames: A list of filenames to check for.
        :type filenames: list, tuple, np.ndarray
        :param all_names: The list of filenames to check in.
        :type all_names: list, tuple, np.ndarray
        
        :return: The indices of filenames within all_names.
        :rtype: list
        """
        inds = []
        for i, f in enumerate(all_names):
            if f in filenames:
                inds.append(i)
        
        return inds


def create_bvc_dict(info_dict, timeseries_dict, bary_dict=None):
    """Create a BVC info dict from modelling results and user input dictionary
    
    :param info_dict: The info dictionary from a :class:'CombinedResults'
        object.
    :type info_dict: dict
    :param timeseries_dict: The timeseries dictionary from a 
        :class:'CombinedResults' object.
    :type timeseries_dict: dict
    :param bary_dict: A user-supplied dictionary with stellar (and possible 
        observatory) information that should be used in the computation of 
        barycentric velocity corrections. If None, the info from the model 
        results is used.
    :type bary_dict: dict, or None
    
    :return: The BVC dictionary.
    :rtype: dict
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    # Setup the BVC info dictionary, first with info from the modelling results.
    out_dict = {
            'star_ra': None,
            'star_dec': None,
            'star_pmra': None,
            'star_pmdec': None,
            'star_rv0': None,
            'star_name': info_dict['star_name'],
            'instrument_lat': info_dict['instrument_lat'],
            'instrument_long': info_dict['instrument_long'],
            'instrument_alt': info_dict['instrument_alt']
            }
    
    # For stellar parameters find the first non-NaN entry and use that (just to make sure).
    for key in ('star_ra', 'star_dec', 'star_pmra', 'star_pmdec'):
        ind = np.where(np.isfinite(timeseries_dict[key]))
        if len(ind[0]) > 0:
            out_dict[key] = timeseries_dict[key][ind[0][0]]
        else:
            out_dict[key] = np.nan
    
    # Now fill the dictionary with the info from the user-supplied bary_dict
    # (if it exists)
    if isinstance(bary_dict, dict):
        logging.info('Using BVC information from user-supplied dictionary...')
        for key, value in bary_dict.items():
            if key in out_dict.keys():
                out_dict[key] = value
            else:
                logging.warning('Key {} in user-supplied dictionary is not useful!'.format(key))
                logging.warning('(Allowed are: {})'.format(out_dict.keys()))
    
    return out_dict
    
    
    
