#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 11:17:04 2022

@author: paul
"""

# Import packages
import pyodine
import pipe_lib

import os
import sys
import time
import numpy as np
from pathos.multiprocessing import Pool
import logging

import argparse
import importlib


def model_single_i2flat(utilities, Pars, obs_file, iod=None, plot_dir=None, 
                        res_names=None, error_log=None, info_log=None, 
                        quiet=False, live=False):
    """Model a flat-field spectrum with I2 cell in the light path
    
    This routine is essentially very similar to the template creation routine,
    as it also models a flat-field spectrum with I2 features (similar to the 
    hot-star spectra), the difference being that no deconvolved template is 
    created after the modelling.
    
    :param utilities: The utilities module for the instrument used in this 
        analysis.
    :type utilities: library
    :param Pars: The parameter input object to use. NOTE: We use the same type
        of class as in the template creation here, as it contains all 
        necessary parameters.
    :type Pars: :class:`Template_Parameters`
    :param obs_file: The pathname of the I2+flat observation to model.
    :type obs_file: str
    :param iod: The I2 template to use in the modelling. If None, it is loaded 
        as specified in the parameter input object.
    :type iod: :class:`IodineTemplate`, or None
    :param plot_dir: The directory name where to save plots and modelling 
        results. If the directory structure does not exist yet, it will be 
        created in the process. If None is given, no results/plots will be 
        saved (default).
    :type plot_dir: str or None
    :param res_names: The pathname under which to save the results file. If you 
        want to save results from multiple runs, you should supply a list with 
        pathnames for each run. If the directory structure does not exist yet, 
        it will be created in the process. If None is given, no results will be 
        saved (default).
    :type res_names: str or list or None
    :param error_log: A pathname of a log-file used for error messages. If 
        None, no errors are logged.
    :type error_log: str, or None
    :param info_log: A pathname of a log-file used for info messages. If 
        None, no info is logged.
    :type info_log: str, or None
    :param quiet: Whether or not to print info messages to terminal. Defaults 
        to False (messages are printed).
    :type quiet: bool
    :param live: If True, then the modelling is performed in live-mode, i.e.
        each modelled chunk is plotted and the best-fit parameters printed to
        terminal. Defaults to False.
    :type live: bool
    """
    
    # Check whether a logger is already setup. If no, setup a new one
    #if not logging.getLogger().hasHandlers():
    #print('Create logger')
    #print(logging.getLogger().hasHandlers())
    pyodine.lib.misc.setup_logging(
            config_file=Pars.log_config_file, level=Pars.log_level,
            error_log=error_log, info_log=info_log, quiet=quiet)
    
    try:
        
        # Start timer
        start_t = time.time()
        
        # Get git information
        branch_name = pyodine.lib.git_check.get_git_branch_name()
        branch_hash = pyodine.lib.git_check.get_git_revision_short_hash()
        
        logging.info('')
        logging.info('Branch: {}'.format(branch_name))
        logging.info('Hash: {}'.format(branch_hash))
        logging.info('---------------------------')
        logging.info('Modelling flat+I2 spectra...')
        logging.info('Working on: {}'.format(obs_file))
        
        ###########################################################################
        ## Set up the environment, and load all neccessary data and parameters
        ###########################################################################
        
        # Load observation
        obs = utilities.load_pyodine.ObservationWrapper(obs_file)
        
        # Load the iodine atlas from file
        iod = utilities.load_pyodine.IodineTemplate(Pars.i2_to_use)
        
        # Output directory for plots (setup the directory structure if non-existent)
        if isinstance(plot_dir, str):
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
        
        # Pathnames for the result files
        if isinstance(res_names, list) and isinstance(res_names[0], str):
            res_names = res_names
        elif isinstance(res_names, str):
            res_names = [res_names]
        else:
            res_names = None
        
        
        ###########################################################################
        ## Now prepare the modelling: Choose the orders, compute weights...
        ###########################################################################
        
        # Choose orders for modelling
        if isinstance(Pars.temp_order_range, (list, tuple)) and Pars.temp_order_range[0] is not None:
            orders = np.arange(Pars.temp_order_range[0], Pars.temp_order_range[1]+1, dtype='int')
        else:
            orders = np.array(obs.orders) 
        logging.info('')
        logging.info('Orders to use: {}'.format(orders))
        
        # Compute weights array for the combined flat+I2 spectrum
        weight = obs.compute_weight(weight_type=Pars.weight_type)
        
        # If the flat+I2 spectra should be normalized prior to fitting, this is
        # done here
        if Pars.normalize_chunks is True:
            for o in obs.orders:
                obs._flux[o] = (obs[o].flux / obs[o].cont)
        
        # Now create the chunks, using the algorithm (and corresponding parameters) 
        # as defined in the parameter input file
        if Pars.chunking_algorithm == 'auto_equal_width':
            obs_chunks = pyodine.chunks.auto_equal_width(
                    obs, width=Pars.chunk_width, orders=orders, padding=Pars.chunk_padding, 
                    chunks_per_order=Pars.chunks_per_order, pix_offset0=Pars.pix_offset0
                    )
        elif Pars.chunking_algorithm == 'wavelength_defined':
            obs_chunks = pyodine.chunks.wavelength_defined(
                    obs, Pars.wavelength_dict, Pars.chunk_padding
                    )
        else:
            raise KeyError('Algorithm {} not known! (Must be one of auto_equal_width, wavelength_defined)'.format(
                    Pars.chunking_algorithm))
        
        nr_chunks_total  = len(obs_chunks)
        nr_chunks_order0 = len(obs_chunks.get_order(obs_chunks.orders[0]))
        nr_orders_chunks = len(obs_chunks.orders)
        
        logging.info('')
        logging.info('Total number of chunks: {}'.format(nr_chunks_total))
        logging.info('Nr. chunks in order 0: {}'.format(nr_chunks_order0))
        logging.info('First and last covered pixel of chunks in order 0: {}, {}'.format(
                obs_chunks[0].abspix[0], obs_chunks[nr_chunks_order0-1].abspix[-1]))
        logging.info('Orders: {} - {} ({} in total)'.format(
                obs_chunks[0].order, obs_chunks[-1].order, nr_orders_chunks))
        
        # Produce the chunk weight array
        chunk_weight = []
        for chunk in obs_chunks:
            chunk_weight.append(weight[chunk.order, chunk.abspix[0]:chunk.abspix[-1]+1])
        #chunk_weight = np.array(chunk_weight)
        
        ###########################################################################
        ## Build and fit the model in as many runs as specified in the 
        ## input parameter file.
        ## Other things specified there:
        ## - which parameters to vary, and possibly limits and starting values
        ## - should output plots be created
        ## - smoothing of wavelength parameters
        ###########################################################################
        
        # Set up the dictionary holding all results from the runs
        run_results = {}
        for run_id in range(len(Pars.model_runs)):
            run_results[run_id] = {
                    'results': [],
                    'fitting_failed': [],
                    'red_chi_sq': np.zeros((nr_chunks_total)),
                    'chunk_w': [], #np.ones((chunk_weight.shape)),
                    'median_pars': pyodine.models.base.ParameterSet(),
                    'std_pars': pyodine.models.base.ParameterSet()
                    }
        
        # Run for run now...
        for run_id, run_dict in Pars.model_runs.items():
            logging.info('')
            logging.info('----------------------')
            logging.info('RUN {}'.format(run_id))
            logging.info('----------------------')
            
            # Build the desired LSF, wavelength and continuum models
            lsf_model = run_dict['lsf_model']
            # For the LSF model: Potentially adapt it to the instrument
            if 'lsf_setup_dict' in run_dict.keys() and isinstance(run_dict['lsf_setup_dict'], dict):
                lsf_model.adapt_LSF(run_dict['lsf_setup_dict'])
            
            wave_model = run_dict['wave_model'] #pyodine.models.wave.LinearWaveModel
            cont_model = run_dict['cont_model'] #pyodine.models.cont.LinearContinuumModel
            
            # If the LSF model is a fixed LSF, try and smooth LSF results from
            # an earlier run
            if lsf_model.name() == 'FixedLSF':
                if 'smooth_lsf_run' not in run_dict.keys():
                    smooth_lsf_run = run_id - 1
                else:
                    smooth_lsf_run = run_dict['smooth_lsf_run']
                # Check if the run LSFs to smooth over actually exist
                if smooth_lsf_run in run_results.keys() and len(run_results[smooth_lsf_run]['results']) == nr_chunks_total:
                    # Smooth the desired LSFs over adjacent chunks and orders
                    logging.info('')
                    logging.info('Smoothing LSF...')
                    manual_redchi2 = None
                    if 'smooth_manual_redchi' in run_dict.keys() and run_dict['smooth_manual_redchi']:
                        manual_redchi2 = run_results[smooth_lsf_run]['red_chi_sq']
                    smooth_osample = None
                    if 'smooth_osample' in run_dict.keys() and run_dict['smooth_osample'] > 0:
                        smooth_osample = int(run_dict['smooth_osample'])
                    
                    lsf_smoothed = pyodine.lib.misc.smooth_lsf(
                            obs_chunks, run_dict['smooth_pixels'], run_dict['smooth_orders'], 
                            run_dict['order_separation'], run_results[smooth_lsf_run]['results'],
                            redchi2=manual_redchi2, osample=smooth_osample)
                    logging.info('LSFs with nans: ', len(np.unique(np.argwhere(np.isnan(lsf_smoothed))[:,0])))
                    
                    LSFarr = pyodine.models.lsf.LSF_Array(lsf_smoothed, np.array([ch.order for ch in obs_chunks]),
                                                          np.array([ch.abspix[0] for ch in obs_chunks]))
                    
                    # Now create the full model, including the fixed LSF
                    model = pyodine.models.spectrum.SimpleModel(
                            lsf_model, wave_model, cont_model, iod, lsf_array=LSFarr,
                            osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
                        
                else:
                    raise KeyError('smooth_lsf_run not in run_results or no results for that run yet!')
            
            else:
                lsf_smoothed = None
                # If no fixed LSF model, set up a "normal" full model
                model = pyodine.models.spectrum.SimpleModel(
                        lsf_model, wave_model, cont_model, iod,
                        osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
            
            # Set up the fitter
            fitter = pyodine.fitters.lmfit_wrapper.LmfitWrapper(model)
            
            # Loop through all chunks to pre-guess the parameters
            starting_pars = []
            for i, chunk in enumerate(obs_chunks):
                starting_pars.append(model.guess_params(chunk))
            
            # Fit wavelength guess within each order if desired (smoother input)
            # (this really only makes sense for the very first run, otherwise
            # use smoothed wavelength results from previous runs)
            if 'pre_wave_slope_deg' in run_dict.keys() and run_dict['pre_wave_slope_deg'] > 0:
                poly_pars = pyodine.lib.misc.smooth_parameters_over_orders(
                        starting_pars, 'wave_slope', obs_chunks, deg=run_dict['pre_wave_slope_deg'])
                # Write the smoothed dispersion into starting_pars
                for i in range(nr_chunks_total):
                    starting_pars[i]['wave_slope'] = poly_pars[i]
                    
            if 'pre_wave_intercept_deg' in run_dict.keys() and run_dict['pre_wave_intercept_deg'] > 0:
                poly_pars = pyodine.lib.misc.smooth_parameters_over_orders(
                        starting_pars, 'wave_intercept', obs_chunks, deg=run_dict['pre_wave_intercept_deg'])
                # Write the smoothed intercepts into starting_pars
                for i in range(nr_chunks_total):
                    starting_pars[i]['wave_intercept'] = poly_pars[i]
            
            # Convert starting parameters to lmfit-parameters 
            lmfit_params = []
            for i in range(nr_chunks_total):
                lmfit_params.append(fitter.convert_params(starting_pars[i], to_lmfit=True))
            
            # Constrain the parameters, using the definitions as supplied in
            # the parameter input file
            lmfit_params = Pars.constrain_parameters(
                    lmfit_params, run_id, run_results, fitter
                    )
            
            ###########################################################################
            ## Finally loop over the chunks and model each of them.
            ###########################################################################
            
            use_chauvenet = False
            if 'use_chauvenet_pixels' in run_dict.keys() and run_dict['use_chauvenet_pixels']:
                use_chauvenet = True
            
            modelling_return = pipe_lib.model_all_chunks(
                    obs_chunks, chunk_weight, fitter, lmfit_params, 
                    use_chauvenet=use_chauvenet, compute_redchi2=True, 
                    use_progressbar=Pars.use_progressbar, live=live)
            
            (run_results[run_id]['results'], run_results[run_id]['chunk_w'], 
             run_results[run_id]['fitting_failed'], chauvenet_outliers, 
             run_results[run_id]['red_chi_sq']) = modelling_return
            
            
            ###########################################################################
            ## Now we do some run diagnostics and analysis of the model results.
            ## Plotting and results saving is done only if a plot_dir was supplied.
            ###########################################################################
            
            # For which chunks could uncertainties not be estimated?! (Why?)
            uncertainties_failed = [i for i in range(nr_chunks_total) if \
                                     'uncertainties could not be estimated' in run_results[run_id]['results'][i].report]
            # For which chunks is the red. Chi2 nan?!
            nan_rchi_fit = [i for i in range(nr_chunks_total) if np.isnan(run_results[run_id]['results'][i].redchi)]
            
            logging.info('')
            logging.info('Number of chunks with no uncertainties: {}'.format(len(uncertainties_failed)))
            logging.info('Number of chunks with outliers: {}'.format(len(chauvenet_outliers)))
            logging.info('Number of chunks with nan fitted red. Chi2: {}'.format(len(nan_rchi_fit)))
            
            # Save the fitting result to file
            # Check in the list for the correct results name
            if isinstance(res_names, list) and 'save_result' in run_dict.keys() and run_dict['save_result'] is True:
                if (len(res_names) - 1) >= run_id:
                    res_save_name = res_names[run_id]
                else:
                    res_save_name = res_names[0]
                # Create the directory structure if it does not exist yet
                res_save_name_dir = os.path.dirname(res_save_name)
                if not os.path.exists(res_save_name_dir) and res_save_name_dir != '':
                    os.makedirs(res_save_name_dir)
                # Save it under the correct file type
                if 'save_filetype' in run_dict.keys() and run_dict['save_filetype'] == 'dill':
                    pyodine.fitters.results_io.save_results(
                            res_save_name, run_results[run_id]['results'], filetype='dill')
                else:
                    pyodine.fitters.results_io.save_results(
                            res_save_name, run_results[run_id]['results'], filetype='h5py')
            
            
            wave_slope_fit = None
            wave_intercept_fit = None
            
            # Fit best-fit wavelength slopes within each order with a polynomial
            if 'wave_slope_deg' in run_dict.keys() and run_dict['wave_slope_deg'] > 0:
                wave_slope_fit = pyodine.lib.misc.smooth_fitresult_over_orders(
                        run_results[run_id]['results'], 'wave_slope', deg=run_dict['wave_slope_deg'])
                run_results[run_id]['wave_slope_fit'] = wave_slope_fit
            
            # Fit best-fit wavelength intercepts within each order with a polynomial
            if 'wave_intercept_deg' in run_dict.keys() and run_dict['wave_intercept_deg'] > 0:
                wave_intercept_fit = pyodine.lib.misc.smooth_fitresult_over_orders(
                        run_results[run_id]['results'], 'wave_intercept', deg=run_dict['wave_intercept_deg'])
                run_results[run_id]['wave_intercept_fit'] = wave_intercept_fit
            
            
            # Calculate median parameter results from this run
            param_names = [k for k in run_results[run_id]['results'][0].params.keys()]
            params = pyodine.models.base.ParameterSet() #{k: np.zeros(nr_chunks_total) for k in param_names}
            errors = pyodine.models.base.ParameterSet() #{k: np.zeros(nr_chunks_total) for k in param_names}
            for p in param_names:
                params[p] = np.array([r.params[p] for r in run_results[run_id]['results']])
                errors[p] = np.array([r.errors[p] for r in run_results[run_id]['results']])
                run_results[run_id]['median_pars'][p] = np.nanmedian(params[p])
                run_results[run_id]['std_pars'][p] = np.nanstd(params[p])
            
            # Write them to file
            if plot_dir and 'save_median_pars' in run_dict.keys() and run_dict['save_median_pars']:
                # Print median parameters into file
                pars_file = os.path.join(plot_dir, 'r{}_median_pars.txt'.format(run_id))
                
                with open(pars_file, 'a') as f:
                    for p in param_names:
                        f.write('\t' + p + ':\t' + str(run_results[run_id]['median_pars'][p]) + \
                                '\t +/- \t' + str(run_results[run_id]['std_pars'][p]) + '\n')
                    f.write('\n')
            
            
            # More plots
            if plot_dir and 'plot_analysis' in run_dict.keys() and run_dict['plot_analysis'] is True:
                # Plot residuals scatter and histogram, red. Chi**2, possibly evaluated chunks,
                # the wavelength slopes and intercepts along with their fits (if accessible),
                # the LSF parameters (if desired), the fitting success (if desired),
                # and some exemplary smoothed LSFs (if accessible)
                plot_chunks = None
                plot_lsf_pars = False
                nr_chunks_order, nr_orders = None, None
                if 'plot_chunks' in run_dict.keys() and isinstance(run_dict['plot_chunks'], (list,tuple)):
                    plot_chunks = run_dict['plot_chunks']
                if 'plot_lsf_pars' in run_dict.keys() and run_dict['plot_lsf_pars']:
                    plot_lsf_pars = True
                if 'plot_success' not in run_dict.keys() or not run_dict['plot_success']:
                    uncertainties_failed = None
                    nan_rchi_fit = None
                    chauvenet_outliers = None
                # If all orders are split into equal number of chunks
                if Pars.chunking_algorithm == 'auto_equal_width':
                    nr_chunks_order = nr_chunks_order0
                    nr_orders = nr_orders_chunks
                
                logging.info('')
                logging.info('Creating analysis plots...')
                
                pipe_lib.create_analysis_plots(
                        run_results[run_id]['results'], plot_dir, run_id=run_id, 
                        red_chi_sq=run_results[run_id]['red_chi_sq'], 
                        nr_chunks_order=nr_chunks_order, nr_orders=nr_orders, 
                        chunk_weight=run_results[run_id]['chunk_w'], plot_chunks=plot_chunks, 
                        chunks=obs_chunks, 
                        wave_intercept_fit=wave_intercept_fit, wave_slope_fit=wave_slope_fit,
                        plot_lsf_pars=plot_lsf_pars,
                        uncertainties_failed=uncertainties_failed,
                        nan_rchi_fit=nan_rchi_fit, chauvenet_outliers=chauvenet_outliers,
                        lsf_array=lsf_smoothed)
            
            ###########################################################################
            ## Run finished, proceeding to next run (unless all through)
            ###########################################################################
        
        ###########################################################################
        ## Now it's done, exit
        ###########################################################################
        
        modelling_time = time.time() - start_t
        logging.info('')
        logging.info('----------------------')
        logging.info('Done, full working time: {}'.format(modelling_time))
    
    except Exception as e:
        logging.error('Something went wrong!', exc_info=True)


def model_multi_i2flat(utilities, Pars, obs_files, plot_dirs=None, 
                       res_files=None, error_files=None, info_files=None, 
                       quiet=False, nr_cores=8):
    """Model multiple flat-field spectra with I2 cell in the light path at the 
    same time
    
    This function can parallelize the modelling of multiple flat+I2 spectra,
    taking advantage of Python's :class:`pathos.multiprocessing.Pool` capabilities.
    The number of parallel processes is defined in the parameter input object.
    
    :param utilities: The utilities module for the instrument used in this 
        analysis.
    :type utilities: library
    :param Pars: The parameter input object for the used instrument. NOTE: We 
        use the same type of class as in the template creation here, as it 
        contains all necessary parameters.
    :type Pars: :class:`Template_Parameters`
    :param obs_files: A pathname to a text-file with pathnames of the flat+I2 
        spectra for the modelling, or, alternatively, a list with the 
        pathnames.
    :type obs_files: str or list
    :param plot_dirs: A pathname to a text-file with directory names for each 
        spectrum where to save plots, or, alternatively, a list with the 
        directory names. If the directory structure does not exist yet, it will 
        be created in the process. If None is given, no results/plots will be 
        saved (default).
    :type plot_dirs: str, list, or None
    :param res_files: A pathname to a text-file with pathnames under which to 
        save the results file(s), or, alternatively, a list with the 
        pathname(s). If you want to save results from multiple runs, you should 
        supply pathnames for each run. If the directory structure does not 
        exist yet, it will be created in the process. If None is given, no 
        results will be saved (default).
    :type res_files: str, list, or None
    :param error_files: A pathname to a text-file with pathnames to log-files 
        used for error messages, or, alternatively, a list with the 
        pathname(s). If None, no errors are logged (default).
    :type error_files: str, list, or None
    :param info_files: A pathname to a text-file with pathnames to log-files 
        used for info messages, or, alternatively, a list with the pathname(s). 
        If None, no info is logged (default).
    :type info_files: str, or None
    :param quiet: Whether or not to print info messages to terminal. Defaults 
        to False (messages are printed).
    :type quiet: bool
    :param nr_cores: The number of cores to use in the parallel modelling. 
        Defaults to 8.
    :type nr_cores: int
    """
    
    ###########################################################################
    ## Some modules and data are loaded here already and then passed to the
    ## individual parallel modelling sessions, so that they do not need to
    ## be loaded in each one individually:
    ## - the I2 template spectrum
    ## - the plot directories and results files
    ## - the log files
    ###########################################################################
    
    # Start timer
    fulltime_start = time.time()
    
    # Load the pathnames of the observations
    if isinstance(obs_files, list):
        obs_names = obs_files
    elif isinstance(obs_files, str):
        with open(obs_files, 'r') as f:
            obs_names = [l.strip() for l in f.readlines()]
    
    # Load the iodine atlas from file
    iod = utilities.load_pyodine.IodineTemplate(Pars.i2_to_use)
    
    # Load the plot directory names for each observation
    if isinstance(plot_dirs, list) and isinstance(plot_dirs[0], str):
        plot_dir_names = plot_dirs
    elif isinstance(plot_dirs, str):
        with open(plot_dirs, 'r') as f:
            plot_dir_names = [l.strip() for l in f.readlines()]
    else:
        plot_dir_names = [None] * len(obs_names)
    
    # Pathnames for the result files
    # Several result names for each run of each observation can be supplied
    # either as a list of lists of names, or a file with several pathnames
    # in each line (for each observation)
    if isinstance(res_files, list):
        res_names = res_files
    elif isinstance(res_files, str):
        res_names = []
        with open(res_files, 'r') as f:
            for l in f.readlines():
                names = l.split() #.strip()
                res_names.append(names)
    else:
        res_names = [None] * len(obs_names)
    
    # Pathnames for the error log files
    if isinstance(error_files, list):
        error_logs = error_files
    elif isinstance(error_files, str):
        error_logs = []
        with open(error_files, 'r') as f:
            for l in f.readlines():
                error_logs.append(l.strip())
    else:
        error_logs = [None] * len(obs_names)
    
    # Pathnames for the info log files
    if isinstance(info_files, list):
        info_logs = info_files
    elif isinstance(info_files, str):
        info_logs = []
        with open(info_files, 'r') as f:
            for l in f.readlines():
                info_logs.append(l.strip())
    else:
        info_logs = [None] * len(obs_names)
    
    
    ###########################################################################
    ## Now parallelize the modelling of the spectra, by initializing the
    ## Pool workers and distribute the files.
    ## Changed this (similar to https://www.py4u.net/discuss/237878)
    ###########################################################################
    
    # Prepare the input arguments list for all the jobs (corresponding
    # to the arguments of the function model_single_i2flat)
    input_arguments = [
            (utilities, Pars, obs_name,) for obs_name in obs_names]
    # Prepare the keyword arguments list for all the jobs (corresponding
    # to the keywords of the function model_single_i2flat)
    input_keywords  = [
            {'iod': iod, 'plot_dir': plot_dir_name, 'res_names': res_name,
             'error_log': error_log, 'info_log': info_log, 'quiet': quiet
             } for plot_dir_name, res_name, error_log, info_log in zip(
             plot_dir_names, res_names, error_logs, info_logs)]
    
    # Setup the Pool object, distribute the arguments and start the jobs
    with Pool(nr_cores) as p:
        jobs = [
            p.apply_async(model_single_i2flat,
                          args=input_arg, kwds=input_kwd
                          ) for input_arg, input_kwd in zip(input_arguments, input_keywords)]
        
        for job in jobs:
            # Wait for the modelling in all workers to finish
            job.wait()
            #job.get()
            # regularly you'd use `job.get()`, but it would `raise` the exception,
            # which is not suitable for this example, so we dig in deeper and just use
            # the `._value` it'd return or raise:
            #print(params, type(job._value), job._value)
    
    # And done
    
    full_modelling_time = time.time() - fulltime_start
    print('\nDone, full working time: ', full_modelling_time)


if __name__ == '__main__':
    
    # Set up the parser for input arguments
    parser = argparse.ArgumentParser(
            description='Model flat-field+I2 spectra')
    
    # Required input arguments:
    # utilities_dir, ostar_files, temp_files, temp_outname, (plot_dir=None, par_file=None)
    parser.add_argument('utilities_dir', type=str, help='The pathname to the utilities directory for this instrument.')
    parser.add_argument('obs_files', type=str, help='A pathname to a text-file with pathnames of flat+I2 spectra for the modelling.')
    parser.add_argument('--plot_dir', type=str, help='A pathname to a text-file with directory names where to save plots.')
    parser.add_argument('--res_files', type=str, help='A pathname to a text-file with pathnames under which to save modelling results.')
    parser.add_argument('--par_file', type=str, help='The pathname of the parameter input file to use.')
    parser.add_argument('--error_files', type=str, help='The pathname to a text-file with pathnames of error log files.')
    parser.add_argument('--info_files', type=str, help='The pathname to a text-file with pathnames of info log files.')
    parser.add_argument('--nr_cores', type=int, help='The number of cores to use in the parallel modelling.')
    parser.add_argument('-q', '--quiet', action='store_true', dest='quiet', help='Do not print messages to the console.')
    
    # Parse the input arguments
    args = parser.parse_args()
    
    utilities_dir = args.utilities_dir
    obs_files     = args.obs_files
    plot_dirs     = args.plot_dirs
    res_files     = args.res_files
    par_file      = args.par_file
    error_files   = args.error_files
    info_files    = args.info_files
    quiet         = args.quiet
    nr_cores      = args.nr_cores
    
    # Import and load the reduction parameters
    if par_file == None:
        module = utilities_dir + '.pyodine_parameters'
        pyodine_parameters = importlib.import_module(module)
        Pars = pyodine_parameters.Template_Parameters()
    else:
        par_file = os.path.splitext(par_file)[0].replace('/', '.')
        pyodine_parameters = importlib.import_module(par_file)
        Pars = pyodine_parameters.Template_Parameters()
    
    # Import and load the utilities
    sys.path.append(os.path.abspath(utilities_dir))
    utilities_dir = utilities_dir.strip('/').split('/')[-1]
    utilities = importlib.import_module(utilities_dir)
    
    # And run the modelling routine
    model_multi_i2flat(utilities, Pars, obs_files, plot_dirs=plot_dirs, 
                        res_files=res_files, error_files=error_files, 
                        info_files=info_files, quiet=quiet, nr_cores=nr_cores)
