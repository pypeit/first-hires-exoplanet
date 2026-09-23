"""
    Set here all the important parameters for the timeseries combination
    of individual chunk results.
    
    Paul Heeren, 8/03/2021
"""

import logging
import os


utilities_dir_path = os.path.dirname(os.path.realpath(__file__))

class Timeseries_Parameters:
    """The control parameters for the chunk velocity combination
    
    In this class you can define details for the chunk velocity combination: 
    Should barycentric velocities be computed, and if so how? Which weighting
    algorithm should be used, with which parameters, and what outputs should be
    generated?
    """
    
    def __init__(self):    
        
        # Logging options
        self.log_config_file = os.path.join(utilities_dir_path, 'logging.json')   # The logging config file
        self.log_level = logging.INFO           # The logging level used for console and info file
        
        # If you hand a list of filenames to reject, are these the names of the
        # individual modelling results ('res_files') or of the original 
        # observations ('obs_files')?
        self.reject_type = 'res_files'
        
        # Should barycentric velocities be computed here using barycorrpy? 
        # And if so, in the precise or predictive way?
        self.compute_bvc = 'precise' # 'precise', 'predictive', or 'no'
        # Whether the code should aim to use the HIP-Nr. from the star_name
        # (otherwise coordinates)
        self.use_hip_for_bvc = True
        # If the observation spectra are of the Sun, you might like to set this
        # to True so that the solar BV correction of barycorrpy is used!
        self.solar_bvc = False
        
        # Which weighting algorithm should be used? Either 'song' or 'lick'.
        # The latter is not tested very well yet.
        self.weighting_algorithm = 'song'
        
        # This dictionary defines the parameters used in the weighting
        # algorithm. For song weighting algorithm:
        # - 'good_chunks' and 'good_orders' define which chunks to use in the
        #   computation of observation means (to offset-correct all chunks)
        # - 'sig_limit_low' and 'sig_limit_up' define the lower and upper limit 
        #   in m/s that chunk timeseries are allowed to have - outliers are 
        #   corrected to 'sig_correct';
        # - 'reweight_alpha', 'reweight_beta' and 'reweight_sigma' are
        #   used in the reweight function to create the chunk weights;
        # - 'weight_correct' is the value that weights of 0 or NaN are 
        #   corrected to.
        # For lick weighting algorithm:
        # - 'percentile': Best percentage of chunks to use in the computation
        #   of weights.
        # - 'maxchi', 'min_counts': Max. red. Chi^2 and minimum counts that 
        #   chunks should have - all others are set to 0 in weights.
        # - 'default_sigma': Default sigma for chunks with fewer than 4 non-zero
        #   chunk weights.
        # - 'useage_percentile': Best percentage of chunks within each observation 
        #   to use in the velocity combination.
        self.weighting_pars_song = {
                'good_chunks': list(range(3,15)), #(150, 350) # chunk indices within orders
                'good_orders': list(range(6,14)),
                'sig_limit_low': 4., 
                'sig_limit_up': 1000.,
                'sig_correct': 1000.,
                'reweight_alpha': 1.8,
                'reweight_beta': 8.0,
                'reweight_sigma': 2.0,
                'weight_correct': 0.01,
                }
        self.weighting_pars_lick = {
                'percentile': 0.997,
                'maxchi': 100000000.,
                'min_counts': 1000.,
                'default_sigma': 1000.,
                'useage_percentile': 0.997,
                }
        
        # Do chromatic index computation?
        self.do_crx = True
        # This dictionary defines the parameters used in the crx modelling:
        # - crx_sigma: If this is not 0, perform sigma-clipping of chunk 
        #   velocities within each observation. Defaults to 0., so no 
        #   sigma-clipping.
        # - crx_iterative: If True, then velocity outliers from the CRX model 
        #   are sigma-clipped iteratively. Defaults to False.
        # - max_iters: If iterative=True, this gives the maximum number of 
        #   iterations to perform in the CRX modelling. Defaults to 10.
        self.crx_pars = {
                'crx_sigma': 0.,
                'crx_iterative': False,
                'crx_max_iters': 10
                }
        
        # For writing timeseries results to a text-file:
        self.txt_outkeys = ['bary_date', 'rv_bc', 'rv_err']     # Write these results
        self.txt_delimiter = '\t'                               # Delimiter to use
        self.txt_header = ''                                    # Header line
        self.txt_outformat = ['%10.5f', '%6.4f', '%3.4f']       # Output format (make sure
                                                                # this matches the keys!)
        self.txt_detailed = True                                # If True, write a detailed
                                                                # output with more metrics
        self.txt_flux_chunk = [251, 252, 253]                   # Chunk indices for estimate
                                                                # of median flux
        
        # Save the final results to file?
        self.save_comb_res = True
        
        # Create and save analysis plots?
        self.plot_analysis = True
        
        # Print observation names of RV outliers? (Only if plot_analysis is True!)
        self.print_outliers = True
        
