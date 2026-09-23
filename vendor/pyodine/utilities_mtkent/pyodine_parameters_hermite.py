"""
    Set here all the important parameters for the I2 reduction pipeline, both
    modeling of individual observations as well as the template creation.
    
    Paul Heeren, 3/02/2021
"""

from pyodine import models

import logging
import os
import sys


utilities_dir_path = os.path.dirname(os.path.realpath(__file__))

###############################################################################
## Here we define the instrument-specific setup of the LSFs. These values can #
## then be used to change the class variables of the respective LSF models.   #
###############################################################################

# For the MultiGaussian model: 11 positions & sigmas for the central Gauss and
# the satellites
_multigauss_setup_dict = {
        'positions': [-3.8, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 3.8],
        #'positions': [-2.9, -2.5, -1.9, -1.4, -1.0, 0.0, 1.0, 1.4, 1.9, 2.5, 2.9],
        'sigmas':    [ 0.9,  0.9,  0.9,  0.9, 1.2, 0.9, 0.9, 0.9, 0.9]
        #'sigmas':    [ 0.9,  0.9,  0.9,  0.9,  0.9, 0.6, 0.9, 0.9, 0.9, 0.9, 0.9]
        }

# For the HermiteGaussian model: Set Hermite degrees 1, 2, and 9 to 0. (This
# is actually the default case, but here to explain how one can adapt it.)
_hermitegauss_setup_dict = {
        'weight_1': 0, 'weight_2': 0, 'weight_3': 1, 'weight_4': 1, 'weight_5': 1, 
        'weight_6': 1, 'weight_7': 1, 'weight_8': 1, 'weight_9': 0
        }


class Parameters:
    """The control commands for the main routine
    
    The exact details of the algorithm are defined entirely by the parameters
    in this class: Parameters for chunk creation, general model parameters,
    and details about how many runs are used in the modelling and which LSF
    models are employed (and more).
    
    Furthermore, in the class method :func:`self.constrain_parameters` you can
    specify and alter input parameter descriptions for the model, e.g. set
    bounds or fix parameters.
    """
    
    def __init__(self):
        
        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                                format='%(message)s')
        
        # General parameters:
        self.osample_obs = 6                    # Oversample factor for the observation modeling
        self.lsf_conv_width = 6.                # LSF is evaluated over this many pixels (times 2)
        self.number_cores = 4                   # Number of processor cores for multiprocessing
        
        self.log_config_file = os.path.join(utilities_dir_path, 'logging.json')   # The logging config file
        self.log_level = logging.INFO           # The logging level used for console and info file
        
        self.use_progressbar = False            # Use a progressbar during chunk modelling?
        
        # Tellurics:
        self.telluric_mask = None               # Telluric mask to use (carmenes, uves or hitran); 
                                                # (None: tellurics are not taken care of)
        self.tell_wave_range = (None,6500)      # Load tellurics only within this wavelength range
        self.tell_dispersion = 0.002            # Dispersion (i.e. wavelength grid) of telluric mask
        
        # Chunking: Which algorithm to use?
        # (currently supported: 'auto_wave_comoving')
        self.chunking_algorithm = 'auto_wave_comoving'
        # If the auto_wave_comoving algorithm is used, the chunks are shifted in wavelengths with
        # respect to the template chunks to account for the change in barycentric velocity. Supply a
        # different value to delta_v in order to define the shift yourself (e.g. 0 for solar observations).
        self.order_range = (None,None)          # Order range (min,max) to use in observation modeling;
                                                # (None,None) uses automatically the same as in the template
        # The chunk width is now determined by the template chunks
        #self.chunk_width = 91                   # Width of chunks in pixels in observation modeling
        self.chunk_padding = 10                 # Padding (left and right) of the chunks in pixels
        self.chunks_per_order = None            # Maximum number of chunks per order (optional)
        self.chunk_delta_v = None               # Velocity shift between template and observation 
                                                # (None: relative barycentric velocity)
        
        # Reference spectrum to use in normalizer and for the first velocity guess
        self.ref_spectrum = 'arcturus'          # Reference spectrum ('arcturus' or 'sun')
        self.velgues_order_range = (43,50)       # Orders used for velocity guess (should be outside I2 region)
        self.delta_v = 1000.                    # The velocity step size for the cross-correlation (in m/s)
        self.maxlag  = 500                      # The number of steps to each side in the cross-correlation
        
        # Normalize chunks in the beginning?
        self.normalize_chunks = False
        
        # Weighting of pixels:
        self.bad_pixel_mask = False             # Whether to run the bad pixel mask
        self.bad_pixel_cutoff = 0.22            # Cutoff parameter for the bad pixel mask
        self.correct_obs = False                # Whether to correct the observation in regions of weight = 0
        self.weight_type = 'flat'               # Type of weights (flat or inverse, as implemented in pyodine.components.Spectrum)
        self.rel_noise = 0.008                  # Only used if weight_type='inverse': The relative noise within a flatfield spectrum
        
        # I2 atlas:
        self.i2_to_use = 3                      # Index of I2 FTS to use (see archive/conf.py)
        self.wavelength_scale = 'air'           # Which wavelength scale to use ('air' or 'vacuum' - should always be the first)
        
        # If you want to create and save velocity analysis plots, put in the desired
        # run number here (these results will be plotted) - else put to None
        self.vel_analysis_plots = -1            # -1 corresponds to the last run
        
        # Now to the run info: For each modelling run, define a new entry in the following dictionary
        # with all the neccessary information needed
        # (except fitting parameters, those are defined further below in constrain_parameters())
        self.model_runs = {
                0:
                {# First define the LSF
                 'lsf_model': models.lsf.SingleGaussian,    # LSF model to use (this is absolutely neccessary)
                 # Then define the wavelength model
                 'wave_model': models.wave.LinearWaveModel,
                 # And define the continuum model
                 'cont_model': models.cont.LinearContinuumModel,
                 
                 # Before the chunks are modeled, you can smooth the wavelength guesses for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs):
                 'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Fitting keywords
                 'use_chauvenet_pixels': True,              # Chauvenet criterion for pixel outliers? (None: False)
                 
                 # Save the fit results from this run?
                 # You can also define the filetype:
                 #    - 'h5py': Saves the most important results to hdf5 (small filesize, harder to recover)
                 #    - 'dill': Saves the whole object structure to pickle (large filesize, easy to recover)
                 'save_result': True,                       # Save the result of this run (None: True)
                 'save_filetype': 'h5py',                   # Filetype to save in (None: 'h5py')
                 # After the chunks have been modeled, you can model the wavelength results for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for next run):
                 'wave_slope_deg': 5,                       # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'wave_intercept_deg': 5,                   # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Plotting keywords
                 'plot_success': True,                      # Create plot of fitting success (None: False)
                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (None: False)
                 'plot_chunks': [150, 250, 400],            # A list with indices of chunks that will be plotted and saved
                 'plot_lsf_pars': True,                     # Plot lsf parameter results (None: False)
                 
                 # Median parameter results
                 'save_median_pars': True,                  # Save median results to text file (None: False)
                 },
                
                1:
                {# First define the LSF
                 'lsf_model': models.lsf.HermiteGaussian,     # LSF model to use (this is absolutely neccessary)
                 'lsf_setup_dict': _hermitegauss_setup_dict,  # The instrument-specific LSF setup parameters
                 # Then define the wavelength model
                 'wave_model': models.wave.LinearWaveModel,
                 # And define the continuum model
                 'cont_model': models.cont.LinearContinuumModel,
                 
                 # Before the chunks are modeled, you can smooth the wavelength guesses for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs):
                 'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Fitting keywords
                 'use_chauvenet_pixels': True,              # Chauvenet criterion for pixel outliers? (None: False)
                 
                 # Save the fit results from this run?
                 # You can also define the filetype:
                 #    - 'h5py': Saves the most important results to hdf5 (small filesize, harder to recover)
                 #    - 'dill': Saves the whole object structure to pickle (large filesize, easy to recover)
                 'save_result': True,                       # Save the result of this run (None: True)
                 'save_filetype': 'h5py',                   # Filetype to save in (None: 'h5py')
                 # After the chunks have been modeled, you can model the wavelength results for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for next run):
                 'wave_slope_deg': 5,                       # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'wave_intercept_deg': 5,                   # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Plotting keywords
                 'plot_success': True,                      # Create plot of fitting success (None: False)
                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (None: False)
                 'plot_chunks': [150, 250, 400],            # A list with indices of chunks that will be plotted and saved
                 'plot_lsf_pars': True,                     # Plot lsf parameter results (None: False)
                 
                 # Median parameter results
                 'save_median_pars': True,                  # Save median results to text file (None: False)
                 }}
#        ,
#                
#                2:
#                {# First define the LSF
#                 'lsf_model': models.lsf.FixedLSF,          # LSF model to use (this is absolutely neccessary)
#                 # For fixed lsf consisting of smoothed lsf results from previous runs,
#                 # define the smoothing parameters here:
#                 'smooth_lsf_run': 1,                       # Smooth lsfs from this run (None: last run)
#                 'smooth_pixels': 160,                      # Pixels (in dispersion direction) to smooth over
#                 'smooth_orders': 3,                        # Orders (in cross-disp direction) to smooth over
#                 'order_separation': 15,                    # Avg. pixels between orders in raw spectrum
#                 'smooth_manual_redchi': False,             # If true, calculate smooth weights from manual redchi2
#                                                            # (otherwise: the lmfit redchi2)
#                 'smooth_osample': 0,                       # Oversampling to use in smoothing 
#                                                            # (None or 0: use the oversampling from the model)
#                
#                 # Then define the wavelength model
#                 'wave_model': models.wave.LinearWaveModel,
#                 # And define the continuum model
#                 'cont_model': models.cont.LinearContinuumModel,
#                 
#                 # Before the chunks are modeled, you can smooth the wavelength guesses for the chunks
#                 # over the orders with polynomials (in order to use the smoothed values as input for the run)
#                 # (probably only makes sense before first run, later use smoothed results from previous runs):
#                 'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (None: 0, no fitting)
#                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept (None: 0, no fitting)
#                 # Fitting keywords
#                 'use_chauvenet_pixels': True,              # Chauvenet criterion for pixel outliers? (None: False)
#                 
#                 # Save the fit results from this run?
#                 # You can also define the filetype:
#                 #    - 'h5py': Saves the most important results to hdf5 (small filesize, harder to recover)
#                 #    - 'dill': Saves the whole object structure to pickle (large filesize, easy to recover)
#                 'save_result': True,                       # Save the result of this run (None: True)
#                 'save_filetype': 'dill',                   # Filetype to save in (None: 'h5py')
#                 # After the chunks have been modeled, you can model the wavelength results for the chunks
#                 # over the orders with polynomials (in order to use the smoothed values as input for next run):
#                 'wave_slope_deg': 3,                       # Polynomial degree of dispersion fitting (None or 0: no fitting)
#                 'wave_intercept_deg': 3,                   # Same as above, for wavelength intercept (None or 0: no fitting)
#                 # Plotting keywords
#                 'plot_success': True,                      # Create plot of fitting success (None: False)
#                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (None: False)
#                 'plot_chunks': [150, 250, 400],            # A list with indices of chunks that will be plotted and saved
#                 'plot_lsf_pars': True,                     # Plot lsf parameter results (None: False)
#                 
#                 # Median parameter results
#                 'save_median_pars': True,                  # Save median results to text file (None: False)
#                 }
#                
#        }
        

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        """Constrain the lmfit_parameters for the models, however you wish!
        
        :param lmfit_params: A list of :class:`lmfit.Parameters` objects
            for each chunk.
        :type lmfit_params: list[:class:`lmfit.Parameters`]
        :param run_id: The current run_id, to allow different definitions
            from run to run.
        :type run_id: int
        :param run_results: Dictionary with important observation info and
            results from previous modelling runs.
        :type run_results: dict
        :param fitter: The fitter used in the modelling.
        :type fitter: :class:`LmfitWrapper`
        
        :return: The updated list of :class:`lmfit.Parameters` objects.
        :rtype: list[:class:`lmfit.Parameters`]
        """
        
        logging.info('')
        logging.info('Constraining parameters for RUN {}'.format(run_id))
        
        ###########################################################################
        # RUN 0
        # Mainly there for first wavelength solution to feed into the next runs.
        # Using a Single-Gaussian model.
        ###########################################################################
        if run_id == 0:
            # Loop over the chunks
            for i in range(len(lmfit_params)):
                # Set velocity to velocity guess (from the reference spectrum)
                lmfit_params[i]['velocity'].set(
                        value=run_results[run_id]['velocity_guess'])
                
                # SingleGaussian model - just constrain the lsf_fwhm
                lmfit_params[i]['lsf_fwhm'].set(
                        value=2.2, min=2.0, max=4.5)
                
                # Constrain the iodine to not become negative (just in case)
                lmfit_params[i]['iod_depth'].set(value=1., vary=False)
                
                lmfit_params[i]['tem_depth'].set(value=1., vary=False)
                
                # If the chunks were normalized beforehand:
                # Fix the continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
        
        ###########################################################################
        # RUN 1
        # A full model now, with the full LSF description (constrained?!).
        # Fit the median LSF result from the first run to find good starting
        # values for this run's LSF (and possibly define bounds).
        # Use (smoothed) wavelength results and continuum results from run 0 
        # as starting values.
        ###########################################################################
        elif run_id == 1:
            # Dictionary of median lsf parameters from previous run
            median_lsf_pars = run_results[0]['median_pars'].filter('lsf')  #{p[4:]: run_results[0]['median_pars'][p] for p in run_results[0]['median_pars'] if 'lsf' in p}
            # Fit the lsf from last run to get good starting parameters
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf_pars)
            
            logging.info('')
            logging.info('Fitted LSF parameters:')
            logging.info(lsf_fit_pars)
            
            # Loop over the chunks
            for i in range(len(lmfit_params)):
                # Set velocity to median velocity from last run
                lmfit_params[i]['velocity'].set(
                        value=run_results[0]['median_pars']['velocity']) #run_results[run_id]['velocity_guess'])
                # Set iodine and template depth to median results from last run
                lmfit_params[i]['iod_depth'].set(value=1., vary=False)
                lmfit_params[i]['tem_depth'].set(value=1., vary=False)
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[0]['wave_slope_fit'][i]) #,run_results[0]['results'][i].params['wave_slope']) #
                        #min=run_results[0]['wave_slope_fit'][i]*0.99,
                        #max=run_results[0]['wave_slope_fit'][i]*1.01)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['wave_intercept_fit'][i])#,run_results[0]['results'][i].params['wave_intercept']) #
                        #min=run_results[0]['wave_intercept_fit'][i]*0.99,
                        #max=run_results[0]['wave_intercept_fit'][i]*1.01)
                
                # Continuum parameters: use results from before
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[0]['results'][i].params['cont_intercept'])
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[0]['results'][i].params['cont_slope'])
                # If the chunks were normalized beforehand:
                # Fix the continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
                
                # Set the LSF parameters to the fitted values,
                # and constrain them (otherwise crazy things can happen)
                for p in lsf_fit_pars.keys():
                    if fitter.model.lsf_model.pars_dict[p]:
                        if abs(lsf_fit_pars[p]) >= 1e-13:
                            lmfit_params[i]['lsf_'+p].set(
                                    value=lsf_fit_pars[p],
                                    min=lsf_fit_pars[p]-abs(lsf_fit_pars[p])*2.,
                                    max=lsf_fit_pars[p]+abs(lsf_fit_pars[p])*2.)
                        else:
                            lmfit_params[i]['lsf_'+p].set(
                                    value=lsf_fit_pars[p],
                                    min=lsf_fit_pars[p]-2e-13,
                                    max=lsf_fit_pars[p]+2e-13)
                    else:
                        lmfit_params[i]['lsf_'+p].set(value=0., vary=False)
        
        return lmfit_params
    

class Parameters2:
    """The control commands for the main routine
    
    The exact details of the algorithm are defined entirely by the parameters
    in this class: Parameters for chunk creation, general model parameters,
    and details about how many runs are used in the modelling and which LSF
    models are employed (and more).
    
    Furthermore, in the class method :func:`self.constrain_parameters` you can
    specify and alter input parameter descriptions for the model, e.g. set
    bounds or fix parameters.
    """
    
    def __init__(self):
        
        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                                format='%(message)s')
        
        # General parameters:
        self.osample_obs = 4                    # Oversample factor for the observation modeling
        self.lsf_conv_width = 6.                # LSF is evaluated over this many pixels (times 2)
        self.number_cores = 4                   # Number of processor cores for multiprocessing
        
        self.log_config_file = os.path.join(utilities_dir_path, 'logging.json')   # The logging config file
        self.log_level = logging.INFO           # The logging level used for console and info file
        
        self.use_progressbar = False            # Use a progressbar during chunk modelling?
        
        # Tellurics:
        self.telluric_mask = None               # Telluric mask to use (carmenes, uves or hitran); 
                                                # (None: tellurics are not taken care of)
        self.tell_wave_range = (None,6500)      # Load tellurics only within this wavelength range
        self.tell_dispersion = 0.002            # Dispersion (i.e. wavelength grid) of telluric mask
        
        # Chunking: Which algorithm to use?
        # (currently supported: 'auto_wave_comoving')
        self.chunking_algorithm = 'auto_wave_comoving'
        # If the auto_wave_comoving algorithm is used, the chunks are shifted in wavelengths with
        # respect to the template chunks to account for the change in barycentric velocity. Supply a
        # different value to delta_v in order to define the shift yourself (e.g. 0 for solar observations).
        self.order_range = (None,None)          # Order range (min,max) to use in observation modeling;
                                                # (None,None) uses automatically the same as in the template
        # The chunk width is now determined by the template chunks
        #self.chunk_width = 91                   # Width of chunks in pixels in observation modeling
        self.chunk_padding = 10                 # Padding (left and right) of the chunks in pixels
        self.chunks_per_order = None            # Maximum number of chunks per order (optional)
        self.chunk_delta_v = None               # Velocity shift between template and observation 
                                                # (None: relative barycentric velocity)
        
        # Reference spectrum to use in normalizer and for the first velocity guess
        self.ref_spectrum = 'arcturus'          # Reference spectrum ('arcturus' or 'sun')
        self.velgues_order_range = (43,50)       # Orders used for velocity guess (should be outside I2 region)
        self.delta_v = 1000.                    # The velocity step size for the cross-correlation (in m/s)
        self.maxlag  = 500                      # The number of steps to each side in the cross-correlation
        
        # Normalize chunks in the beginning?
        self.normalize_chunks = False
        
        # Weighting of pixels:
        self.bad_pixel_mask = False             # Whether to run the bad pixel mask
        self.bad_pixel_cutoff = 0.22            # Cutoff parameter for the bad pixel mask
        self.correct_obs = False                # Whether to correct the observation in regions of weight = 0
        self.weight_type = 'flat'               # Type of weights (flat or inverse, as implemented in pyodine.components.Spectrum)
        self.rel_noise = 0.008                  # Only used if weight_type='inverse': The relative noise within a flatfield spectrum
        
        # I2 atlas:
        self.i2_to_use = 3                      # Index of I2 FTS to use (see archive/conf.py)
        self.wavelength_scale = 'air'           # Which wavelength scale to use ('air' or 'vacuum' - should always be the first)
        
        # If you want to create and save velocity analysis plots, put in the desired
        # run number here (these results will be plotted) - else put to None
        self.vel_analysis_plots = -1            # -1 corresponds to the last run
        
        # Now to the run info: For each modelling run, define a new entry in the following dictionary
        # with all the neccessary information needed
        # (except fitting parameters, those are defined further below in constrain_parameters())
        self.model_runs = {
                0:
                {# First define the LSF
                 'lsf_model': models.lsf.SingleGaussian,    # LSF model to use (this is absolutely neccessary)
                 # Then define the wavelength model
                 'wave_model': models.wave.LinearWaveModel,
                 # And define the continuum model
                 'cont_model': models.cont.LinearContinuumModel,
                 
                 # Before the chunks are modeled, you can smooth the wavelength guesses for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs):
                 'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Fitting keywords
                 'use_chauvenet_pixels': True,              # Chauvenet criterion for pixel outliers? (None: False)
                 
                 # Save the fit results from this run?
                 # You can also define the filetype:
                 #    - 'h5py': Saves the most important results to hdf5 (small filesize, harder to recover)
                 #    - 'dill': Saves the whole object structure to pickle (large filesize, easy to recover)
                 'save_result': True,                       # Save the result of this run (None: True)
                 'save_filetype': 'h5py',                   # Filetype to save in (None: 'h5py')
                 # After the chunks have been modeled, you can model the wavelength results for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for next run):
                 'wave_slope_deg': 5,                       # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'wave_intercept_deg': 5,                   # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Plotting keywords
                 'plot_success': True,                      # Create plot of fitting success (None: False)
                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (None: False)
                 'plot_chunks': [150, 250, 400],            # A list with indices of chunks that will be plotted and saved
                 'plot_lsf_pars': True,                     # Plot lsf parameter results (None: False)
                 
                 # Median parameter results
                 'save_median_pars': True,                  # Save median results to text file (None: False)
                 },
                
                1:
                {# First define the LSF
                 'lsf_model': models.lsf.HermiteGaussian,     # LSF model to use (this is absolutely neccessary)
                 'lsf_setup_dict': _hermitegauss_setup_dict,  # The instrument-specific LSF setup parameters
                 # Then define the wavelength model
                 'wave_model': models.wave.LinearWaveModel,
                 # And define the continuum model
                 'cont_model': models.cont.LinearContinuumModel,
                 
                 # Before the chunks are modeled, you can smooth the wavelength guesses for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs):
                 'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Fitting keywords
                 'use_chauvenet_pixels': True,              # Chauvenet criterion for pixel outliers? (None: False)
                 
                 # Save the fit results from this run?
                 # You can also define the filetype:
                 #    - 'h5py': Saves the most important results to hdf5 (small filesize, harder to recover)
                 #    - 'dill': Saves the whole object structure to pickle (large filesize, easy to recover)
                 'save_result': True,                       # Save the result of this run (None: True)
                 'save_filetype': 'h5py',                   # Filetype to save in (None: 'h5py')
                 # After the chunks have been modeled, you can model the wavelength results for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for next run):
                 'wave_slope_deg': 5,                       # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'wave_intercept_deg': 5,                   # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Plotting keywords
                 'plot_success': True,                      # Create plot of fitting success (None: False)
                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (None: False)
                 'plot_chunks': [150, 250, 400],            # A list with indices of chunks that will be plotted and saved
                 'plot_lsf_pars': True,                     # Plot lsf parameter results (None: False)
                 
                 # Median parameter results
                 'save_median_pars': True,                  # Save median results to text file (None: False)
                 }}
        

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        """Constrain the lmfit_parameters for the models, however you wish!
        
        :param lmfit_params: A list of :class:`lmfit.Parameters` objects
            for each chunk.
        :type lmfit_params: list[:class:`lmfit.Parameters`]
        :param run_id: The current run_id, to allow different definitions
            from run to run.
        :type run_id: int
        :param run_results: Dictionary with important observation info and
            results from previous modelling runs.
        :type run_results: dict
        :param fitter: The fitter used in the modelling.
        :type fitter: :class:`LmfitWrapper`
        
        :return: The updated list of :class:`lmfit.Parameters` objects.
        :rtype: list[:class:`lmfit.Parameters`]
        """
        
        logging.info('')
        logging.info('Constraining parameters for RUN {}'.format(run_id))
        
        ###########################################################################
        # RUN 0
        # Mainly there for first wavelength solution to feed into the next runs.
        # Using a Single-Gaussian model.
        ###########################################################################
        if run_id == 0:
            # Loop over the chunks
            for i in range(len(lmfit_params)):
                # Set velocity to velocity guess (from the reference spectrum)
                lmfit_params[i]['velocity'].set(
                        value=run_results[run_id]['velocity_guess'])
                
                # SingleGaussian model - just constrain the lsf_fwhm
                lmfit_params[i]['lsf_fwhm'].set(
                        value=2.2, min=2.0, max=4.5)
                
                # Constrain the iodine to not become negative (just in case)
                lmfit_params[i]['iod_depth'].set(value=1., min=0.9, max=1.1)
                
                lmfit_params[i]['tem_depth'].set(value=1., min=0.9, max=1.1)
                
                # If the chunks were normalized beforehand:
                # Fix the continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
        
        ###########################################################################
        # RUN 1
        # A full model now, with the full LSF description (constrained?!).
        # Fit the median LSF result from the first run to find good starting
        # values for this run's LSF (and possibly define bounds).
        # Use (smoothed) wavelength results and continuum results from run 0 
        # as starting values.
        ###########################################################################
        elif run_id == 1:
            # Dictionary of median lsf parameters from previous run
            median_lsf_pars = run_results[0]['median_pars'].filter('lsf')  #{p[4:]: run_results[0]['median_pars'][p] for p in run_results[0]['median_pars'] if 'lsf' in p}
            # Fit the lsf from last run to get good starting parameters
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf_pars)
            
            logging.info('')
            logging.info('Fitted LSF parameters:')
            logging.info(lsf_fit_pars)
            
            # Loop over the chunks
            for i in range(len(lmfit_params)):
                # Set velocity to median velocity from last run
                lmfit_params[i]['velocity'].set(
                        value=run_results[0]['median_pars']['velocity']) #run_results[run_id]['velocity_guess'])
                # Set iodine and template depth to median results from last run
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[0]['median_pars']['iod_depth'], min=0.9, max=1.1)
                lmfit_params[i]['tem_depth'].set(
                        value=run_results[0]['median_pars']['tem_depth'], min=0.9, max=1.1)
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[0]['wave_slope_fit'][i]) #,run_results[0]['results'][i].params['wave_slope']) #
                        #min=run_results[0]['wave_slope_fit'][i]*0.99,
                        #max=run_results[0]['wave_slope_fit'][i]*1.01)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['wave_intercept_fit'][i])#,run_results[0]['results'][i].params['wave_intercept']) #
                        #min=run_results[0]['wave_intercept_fit'][i]*0.99,
                        #max=run_results[0]['wave_intercept_fit'][i]*1.01)
                
                # Continuum parameters: use results from before
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[0]['results'][i].params['cont_intercept'])
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[0]['results'][i].params['cont_slope'])
                # If the chunks were normalized beforehand:
                # Fix the continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
                
                # Set the LSF parameters to the fitted values,
                # and constrain them (otherwise crazy things can happen)
                for p in lsf_fit_pars.keys():
                    lmfit_params[i]['lsf_'+p].set(
                        value=lsf_fit_pars[p],
                        min=lsf_fit_pars[p]-abs(lsf_fit_pars[p])*0.4,
                        max=lsf_fit_pars[p]+abs(lsf_fit_pars[p])*0.4)
        
        return lmfit_params

    
class Template_Parameters:
    """The control commands for the main template creation routine
    
    This is mostly the same as the :class:`Parameters` class, but with some
    extra parameters essential for the deconvolution and template generation.
    """
    
    def __init__(self):
        
        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                                format='%(message)s')
        
        # General parameters:
        self.osample_obs = 6                    # Oversample factor for the observation modeling
        self.lsf_conv_width = 6.                # LSF is evaluated over this many pixels (times 2)
        
        self.log_config_file = os.path.join(utilities_dir_path, 'logging.json')   # The logging config file
        self.log_level = logging.INFO           # The logging level used for console and info file
        
        self.use_progressbar = True             # Use a progressbar during chunk modelling?
        
        # Tellurics:
        self.telluric_mask = None               # Telluric mask to use (carmenes, uves or hitran); 
                                                # if None, tellurics are not taken care of
        self.tell_wave_range = (None,6500)      # Load tellurics only within this wavelength range
        self.tell_dispersion = 0.002            # Dispersion (i.e. wavelength grid) of telluric mask
        
        # Chunking: Which algorithm to use?
        # (currently supported: 'auto_equal_width' or 'wavelength_defined')
        self.chunking_algorithm = 'auto_equal_width'
        # If the auto_equal_width chunking algorithm is used, the chunks are defined by the user
        # through their width, padding, number of chunks per order, and pixel offset of the first chunk:
        self.temp_order_range = (22,43)         # Order range (min,max) to use in observation modeling;
                                                # (None,None) uses all orders
        self.chunk_width = 125                  # Width of chunks in pixels in observation modeling
        self.chunk_padding = 45                 # Padding (left and right) of the chunks in pixels
        self.chunks_per_order = 31              # Maximum number of chunks per order (optional)
        self.pix_offset0 = 100                  # The starting pixel of the first chunk within each order
                                                # (None: the chunks will be centered within orders)
        # Otherwise, if the wavelength_defined chunking algorithm is chosen, make sure you have
        # added a dictionary with start and end wavelengths for each chunk (see below constrain_parameters())
        self.wavelength_dict = self.chunk_wavelengths
        
        # Reference spectrum to use in normalizer and for the first velocity guess
        self.ref_spectrum = 'arcturus'          # Reference spectrum ('arcturus' or 'sun')
        self.velgues_order_range = (20,40)      # Orders used for velocity guess (should be outside I2 region)
        self.delta_v = 1000.                    # The velocity step size for the cross-correlation (in m/s)
        self.maxlag  = 500                      # The number of steps to each side in the cross-correlation
        
        # Normalize chunks in the beginning?
        self.normalize_chunks = False
        
        # Weighting of pixels:
        self.bad_pixel_mask = False             # Whether to run the bad pixel mask
        self.bad_pixel_cutoff = 0.22            # Cutoff parameter for the bad pixel mask
        self.correct_obs = False                # Whether to correct the observation in regions of weight = 0.
        self.weight_type = 'flat'               # Type of weights (flat or lick, as implemented in pyodine.components.Spectrum)
        
        # I2 atlas:
        self.i2_to_use = 3                      # Index of I2 FTS to use (see archive/conf.py)
        self.wavelength_scale = 'air'           # Which wavelength scale to use ('air' or 'vacuum' - should always be the first)
        
        # The parameters for the Jansson deconvolution algorithm.
        self.jansson_run_model = 1              # Model (LSF, wave, cont) from this run used in deconvolution
        self.chunk_weights_redchi = False       # Use fitting red.Chi2 as chunk weights for template? (otherwise analytic)
        # Deconvolution parameters
        self.deconvolution_pars = {
                'osample_temp': 10,             # Oversampling of template
                'jansson_niter': 1200,          # Max. number of iterations in Jansson deconvolution
                'jansson_zerolevel': 0.00,      # Spectrum zero-level in Jansson deconvolution
                'jansson_contlevel': 1.02,      # Spectrum continuum-level in Jansson deconvolution
                'jansson_conver': 0.2,          # Convergence parameter in Jansson deconvolution (careful with that!)
                'jansson_chi_change': 1e-6,     # Minimum change of red.Chi**2 in Jansson deconvolution after 
                                                # which iterations are stopped
                'lsf_conv_width': self.lsf_conv_width,  # The LSF is evaluated over this many pixels (times 2)
                }
        
        # If a smoothed LSF is used as input for the Jansson deconvolution, these
        # parameters specify it
        self.jansson_lsf_smoothing = {
                'do_smoothing': False,          # If False, do no LSF smoothing
                'smooth_lsf_run': 1,            # Smooth the LSFs from this run (None: last run)
                'smooth_pixels': 160,           # Pixels (in dispersion direction) to smooth over
                'smooth_orders': 3,             # Orders (in cross-disp direction) to smooth over
                'order_separation': 15,         # Avg. pixels between orders in raw spectrum
                'smooth_manual_redchi': False,  # If true, calculate smooth weights from manual redchi2
                                                # (otherwise: the lmfit redchi2)
                }
        
        
        # Now to the run info: For each modelling run, define a new entry in the following dictionary
        # with all the neccessary information needed
        # (except fitting parameters, those are defined further below in constrain_parameters())
        self.model_runs = {
                0:
                {# First define the LSF
                 'lsf_model': models.lsf.SingleGaussian,    # LSF model to use (this is absolutely neccessary)
                 # Then define the wavelength model
                 'wave_model': models.wave.LinearWaveModel,
                 # And define the continuum model
                 'cont_model': models.cont.LinearContinuumModel,
                 
                 # Before the chunks are modeled, you can smooth the wavelength guesses for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs):
                 'pre_wave_slope_deg': 3,                   # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'pre_wave_intercept_deg': 3,               # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Fitting keywords
                 'use_chauvenet_pixels': True,              # Chauvenet criterion for pixel outliers? (None: False)
                 
                 # Save the fit results from this run?
                 # You can also define the filetype:
                 #    - 'h5py': Saves the most important results to hdf5 (small filesize, harder to recover)
                 #    - 'dill': Saves the whole object structure to pickle (large filesize, easy to recover)
                 'save_result': True,                       # Save the result of this run (None: True)
                 'save_filetype': 'h5py',                   # Filetype to save in (None: 'h5py')
                 # After the chunks have been modeled, you can model the wavelength results for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for next run):
                 'wave_slope_deg': 3,                       # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'wave_intercept_deg': 3,                   # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Plotting keywords
                 'plot_success': True,                      # Create plot of fitting success (None: False)
                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (None: False)
                 'plot_chunks': [150, 250, 400],            # A list with indices of chunks that will be plotted and saved
                 'plot_lsf_pars': True,                     # Plot lsf parameter results (None: False)
                 
                 # Median parameter results
                 'save_median_pars': True,                  # Save median results to text file (None: False)
                 },
                
                1:
                {# First define the LSF
                 'lsf_model': models.lsf.MultiGaussian,     # LSF model to use (this is absolutely neccessary)
                 'lsf_setup_dict': _multigauss_setup_dict,  # The instrument-specific LSF setup parameters
                 # Then define the wavelength model
                 'wave_model': models.wave.LinearWaveModel,
                 # And define the continuum model
                 'cont_model': models.cont.LinearContinuumModel,
                 
                 # Before the chunks are modeled, you can smooth the wavelength guesses for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs):
                 'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Fitting keywords
                 'use_chauvenet_pixels': True,              # Chauvenet criterion for pixel outliers? (None: False)
                 
                 # Save the fit results from this run?
                 # You can also define the filetype:
                 #    - 'h5py': Saves the most important results to hdf5 (small filesize, harder to recover)
                 #    - 'dill': Saves the whole object structure to pickle (large filesize, easy to recover)
                 'save_result': True,                       # Save the result of this run (None: True)
                 'save_filetype': 'dill',                   # Filetype to save in (None: 'h5py')
                 # After the chunks have been modeled, you can model the wavelength results for the chunks
                 # over the orders with polynomials (in order to use the smoothed values as input for next run):
                 'wave_slope_deg': 3,                       # Polynomial degree of dispersion fitting (None or 0: no fitting)
                 'wave_intercept_deg': 3,                   # Same as above, for wavelength intercept (None or 0: no fitting)
                 # Plotting keywords
                 'plot_success': True,                      # Create plot of fitting success (None: False)
                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (None: False)
                 'plot_chunks': [150, 250, 400],            # A list with indices of chunks that will be plotted and saved
                 'plot_lsf_pars': True,                     # Plot lsf parameter results (None: False)
                 
                 # Median parameter results
                 'save_median_pars': True,                  # Save median results to text file (None: False)
                 }}
        #,
        #        
        #        2:
        #        {# First define the LSF
        #         'lsf_model': models.lsf.FixedLSF,          # LSF model to use (this is absolutely neccessary)
        #         # For fixed lsf consisting of smoothed lsf results from previous runs,
        #         # define the smoothing parameters here:
        #         'smooth_lsf_run': 1,                       # Smooth lsfs from this run (None: last run)
        #         'smooth_pixels': 160,                      # Pixels (in dispersion direction) to smooth over
        #         'smooth_orders': 3,                        # Orders (in cross-disp direction) to smooth over
        #         'order_separation': 15,                    # Avg. pixels between orders in raw spectrum
        #         'smooth_manual_redchi': False,             # If true, calculate smooth weights from manual redchi2
        #                                                    # (otherwise: the lmfit redchi2)
        #         'smooth_osample': 0,                       # Oversampling to use in smoothing 
        #                                                    # (None or 0: use the oversampling from the model)
        #        
        #         # Then define the wavelength model
        #         'wave_model': models.wave.LinearWaveModel,
        #         # And define the continuum model
        #         'cont_model': models.cont.LinearContinuumModel,
        #         
        #         # Before the chunks are modeled, you can smooth the wavelength guesses for the chunks
        #         # over the orders with polynomials (in order to use the smoothed values as input for the run)
        #         # (probably only makes sense before first run, later use smoothed results from previous runs):
        #         'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (None: 0, no fitting)
        #         'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept (None: 0, no fitting)
        #         # Fitting keywords
        #         'use_chauvenet_pixels': True,              # Chauvenet criterion for pixel outliers? (None: False)
        #         
        #         # Save the fit results from this run?
        #         # You can also define the filetype:
        #         #    - 'h5py': Saves the most important results to hdf5 (small filesize, harder to recover)
        #         #    - 'dill': Saves the whole object structure to pickle (large filesize, easy to recover)
        #         'save_result': True,                       # Save the result of this run (None: True)
        #         'save_filetype': 'dill',                   # Filetype to save in (None: 'h5py')
        #         # After the chunks have been modeled, you can model the wavelength results for the chunks
        #         # over the orders with polynomials (in order to use the smoothed values as input for next run):
        #         'wave_slope_deg': 3,                       # Polynomial degree of dispersion fitting (None or 0: no fitting)
        #         'wave_intercept_deg': 3,                   # Same as above, for wavelength intercept (None or 0: no fitting)
        #         # Plotting keywords
        #         'plot_success': True,                      # Create plot of fitting success (None: False)
        #         'plot_analysis': True,                     # Create analysis plots (residuals etc.) (None: False)
        #         'plot_chunks': [150, 250, 400],            # A list with indices of chunks that will be plotted and saved
        #         'plot_lsf_pars': True,                     # Plot lsf parameter results (None: False)
        #         
        #         # Median parameter results
        #         'save_median_pars': True,                  # Save median results to text file (None: False)
        #         }
        #        
        #}
        

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        """Constrain the lmfit_parameters for the models, however you wish!
        
        :param lmfit_params: A list of :class:`lmfit.Parameters` objects
            for each chunk.
        :type lmfit_params: list[:class:`lmfit.Parameters`]
        :param run_id: The current run_id, to allow different definitions
            from run to run.
        :type run_id: int
        :param run_results: Dictionary with important observation info and
            results from previous modelling runs.
        :type run_results: dict
        :param fitter: The fitter used in the modelling.
        :type fitter: :class:`LmfitWrapper`
        
        :return: The updated list of :class:`lmfit.Parameters` objects.
        :rtype: list[:class:`lmfit.Parameters`]
        """
        
        logging.info('')
        logging.info('Constraining parameters for RUN {}'.format(run_id))
        
        ###########################################################################
        # RUN 0
        # Mainly there for first wavelength solution to feed into the next runs.
        # Using a Single-Gaussian model.
        ###########################################################################
        if run_id == 0:
            # Loop over the chunks
            for i in range(len(lmfit_params)):
                # Fix velocity and template depth in template creation
                # (these are put to constant 0. and 1. then)
                lmfit_params[i]['velocity'].set(
                        vary=False)
                lmfit_params[i]['tem_depth'].set(
                        vary=False)
                
                # SingleGaussian model - just constrain the lsf_fwhm
                lmfit_params[i]['lsf_fwhm'].set(
                        value=3.0, min=2.0, max=4.5)
                
                # Constrain the iodine to not become negative (just in case)
                lmfit_params[i]['iod_depth'].set(min=0.1)
                
                # If the chunks were normalized beforehand:
                # Fix the continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
        
        ###########################################################################
        # RUN 1
        # A full model now, with the full LSF description (constrained?!).
        # Fit the median LSF result from the first run to find good starting
        # values for this run's LSF (and possibly define bounds).
        # Use (smoothed) wavelength results and continuum results from run 0 
        # as starting values.
        ###########################################################################
        elif run_id == 1:
            # Dictionary of median lsf parameters from previous run
            median_lsf_pars = run_results[0]['median_pars'].filter('lsf')  #{p[4:]: run_results[0]['median_pars'][p] for p in run_results[0]['median_pars'] if 'lsf' in p}
            # Fit the lsf from last run to get good starting parameters
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf_pars)
            
            logging.info('')
            logging.info('Fitted LSF parameters:')
            logging.info(lsf_fit_pars)
            
            # Loop over the chunks
            for i in range(len(lmfit_params)):
                # Fix velocity and template depth in template creation
                # (these are put to constant 0. and 1. then)
                lmfit_params[i]['velocity'].set(
                        vary=False)
                lmfit_params[i]['tem_depth'].set(
                        vary=False)
                
                # Set iodine depth to median results from last run
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[0]['median_pars']['iod_depth'])
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[0]['results'][i].params['wave_slope']) #run_results[0]['wave_slope_fit'][i])#), #,
                        #min=run_results[0]['wave_slope_fit'][i]*0.99,
                        #max=run_results[0]['wave_slope_fit'][i]*1.01)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['results'][i].params['wave_intercept']) #run_results[0]['wave_intercept_fit'][i])#,
                        #min=run_results[0]['wave_intercept_fit'][i]*0.99,
                        #max=run_results[0]['wave_intercept_fit'][i]*1.01)
                
                # Continuum parameters: use results from before
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[0]['results'][i].params['cont_intercept'])
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[0]['results'][i].params['cont_slope'])
                # If the chunks were normalized beforehand:
                # Fix the continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
                
                # Set the LSF parameters to the fitted values,
                # and constrain them (otherwise crazy things can happen)
                for p in lsf_fit_pars.keys():
                    if fitter.model.lsf_model.pars_dict[p]:
                        if abs(lsf_fit_pars[p]) >= 1e-13:
                            lmfit_params[i]['lsf_'+p].set(
                                    value=lsf_fit_pars[p],
                                    min=lsf_fit_pars[p]-abs(lsf_fit_pars[p])*2.,
                                    max=lsf_fit_pars[p]+abs(lsf_fit_pars[p])*2.)
                        else:
                            lmfit_params[i]['lsf_'+p].set(
                                    value=lsf_fit_pars[p],
                                    min=lsf_fit_pars[p]-2e-13,
                                    max=lsf_fit_pars[p]+2e-13)
                    else:
                        lmfit_params[i]['lsf_'+p].set(value=0., vary=False)
        """
        ###########################################################################
        # RUN 2
        # Final run. Use smoothed LSF results from run 1 and keep them fixed.
        # Only vary other parameters (use results from run 1 as starting values).
        ###########################################################################
        elif run_id == 2:
            # Loop over the chunks
            for i in range(len(lmfit_params)):
                # Fix velocity and template depth in template creation
                # (these are put to constant 0. and 1. then)
                lmfit_params[i]['velocity'].set(
                        vary=False)
                lmfit_params[i]['tem_depth'].set(
                        vary=False)
                
                # Set iodine depth to median results from last run
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[1]['median_pars']['iod_depth'])
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[1]['wave_slope_fit'][i]),
                        #min=run_results[1]['wave_slope_fit'][i]*0.99,
                        #max=run_results[1]['wave_slope_fit'][i]*1.01)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[1]['wave_intercept_fit'][i])
                        #min=run_results[1]['wave_intercept_fit'][i]*0.96,
                        #max=run_results[1]['wave_intercept_fit'][i]*1.04)
                        
                # Continuum parameters: use results from before
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[1]['results'][i].params['cont_intercept'])
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[1]['results'][i].params['cont_slope'])
                # If the chunks were normalized beforehand:
                # Fix the continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
                
                # For fixed, smoothed LSF: Don't vary order and pixel0 values!
                # (and better also not the amplitude)
                lmfit_params[i]['lsf_order'].vary = False
                lmfit_params[i]['lsf_pixel0'].vary = False
                lmfit_params[i]['lsf_amplitude'].vary = False
        """
        return lmfit_params

    # Start and end wavelengths of chunks if wave_defined chunking algorithm is used
    # in template creation
    chunk_wavelengths = {
            'start_wave': [],
            'end_wave': []
            }