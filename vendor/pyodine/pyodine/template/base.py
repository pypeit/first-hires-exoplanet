import os.path
import h5py
import astropy.time
import numpy as np
import logging
import sys

from ..components import MultiOrderSpectrum, Spectrum, Observation, NoDataError, TemplateChunk


class StellarTemplate(MultiOrderSpectrum):
    """A deconvolved stellar template, with chunks stitched together
    
    A :class:`StellarTemplate` is always represented as a 
    :class:`MultiOrderSpectrum`, even if there is only one order.
    To be subclassed together with either :class:`Spectrum` or 
    :class:`MultiOrderSpectrum`.
    
    NOTE: THIS HAS NOT BEEN USED OR TESTED IN THE LATER STAGES OF DEVELOPMENT
    OF THE SOFTWARE!
    
    :param observation: Hand a :class:`Observation` object upon template 
        creation. Otherwise, if a pathname is given, the 
        :class:`StellarTemplate` is loaded from there.
    :type observation: :class:`Observation`, or str
    :param velocity_offset: Upon template creation, hand the velocity-offset 
        between the observation and reference spectrum. Leave as None when 
        loading from file (default).
    :type velocity_offset: float, or None
    :param bary_vel_corr: Upon template creation, hand the barycentric velocity 
        of the observation. Leave as None when loading from file (default).
    :type bary_vel_corr: float, or None
    :param osample: Upon template creation, hand the oversampling factor used 
        in the template creation. Leave as None when loading from file 
        (default).
    :type osample: int, or None
    
    """

    def __init__(self, observation, velocity_offset=None, bary_vel_corr=None, 
                 osample=None):
        
        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                                format='%(message)s')
        
        self._data = {}
        # Load existing template, if first argument is a string
        if type(observation) is str:
            logging.info('Loading template from file: {}'.format(observation))
            
            try:
                with h5py.File(observation, 'r') as h:
                    try:
                        self.orig_filename = h['/orig_filename'][()]
                    except:
                        self.orig_filename = os.path.abspath(observation)
                    self.starname = h['/starname'][()]
                    self.time_start = astropy.time.Time(
                        h['/time_start'][()],
                        format='isot',
                        scale='utc'
                    )
                    self.velocity_offset = h['/velocity_offset'][()]
                    self.bary_vel_corr = h['/bary_vel_corr'][()]
                    self.osample = h['/osample'][()]
                    orders = h['/orders']
                    for i in orders:
                        wave = h['/data/{}/wave'.format(i)][()]
                        flux = h['/data/{}/flux'.format(i)][()]
                        self._data[i] = Spectrum(flux, wave)
            except Exception as e:
                raise(e)
        # Initialize template from observation
        elif isinstance(observation, Observation) and velocity_offset is not None:
            # FIXME: Handle different velocity offsets (relative to Arcturus atlas etc.)
            self.orig_filename = None
            self.velocity_offset = velocity_offset
            self.starname = observation.star.name
            self.time_start = observation.time_start  # Start time of the observation <astropy.time.Time>
            self.bary_vel_corr = bary_vel_corr
            if self.bary_vel_corr is None:
                self.bary_vel_corr = observation[0].bary_vel_corr
            self.osample = osample
        else:
            raise KeyError("""Invalid input arguments. Usage:
                1) StellarTemplate('filename.h5')
                2) StellarTemplate(template_observation, velocity_offset)""")

    @property
    def orders(self):
        """The order numbers of the template
        """
        return list(self._data.keys())

    def save(self, filename):
        """Save as HDF5 file (.h5)
        
        :param filename: The pathname of the directory where to save the
            file, or the filename itself.
        :type filename: str
        """
        
        logging.info('Saving deconvolved template to {}'.format(filename))
        
        if os.path.isdir(filename):
            date = self.time_start.datetime.strftime("%Y%m%d")
            filename = os.path.join(filename, self.starname + '_' + date + '.h5')
            logging.warning('Supplied filename was a directory.')
            logging.warning('New filename: {}'.format(filename))
            
        with h5py.File(filename, 'w') as h:
            h.create_dataset('/orig_filename', data=os.path.abspath(filename))
            h.create_dataset('/starname', data=self.starname)
            h.create_dataset('/time_start', data=self.time_start.isot)
            h.create_dataset('/velocity_offset', data=self.velocity_offset)
            h.create_dataset('/bary_vel_corr', data=self.bary_vel_corr)
            h.create_dataset('/osample', data=self.osample)
            h.create_dataset('/orders', data=self.orders)
            # TODO: Save timestamp
            for i in self.orders:
                grp = '/data/{}/'.format(i)
                h.create_dataset(grp + 'flux', data=self[i].flux)
                h.create_dataset(grp + 'wave', data=self[i].wave)

    def __len__(self):
        return len(self._data)

    def __getitem__(self, order):
        try:
            return self._data[order]  # TODO: Deepcopy this?
        except:
            logging.error(NoDataError('No data for order {}!'.format(order)))

    def __setitem__(self, order, spectrum):
        self._data[order] = spectrum  # TODO: Deepcopy this?

    def __str__(self):
        return '<StellarTemplate of {} ({} orders)>'.format(
            self.starname,
            len(self)
        )


class StellarTemplate_Chunked:
    """A deconvolved stellar template, with chunks remaining separated
    
    The :class:`StellarTemplate_Chunked` is represented as a list of individual
    :class:`TemplateChunks`, with additional parameters and methods.
    Used for I2 reduction similar as in Lick dop code, where observation
    chunks are defined over the same wavelengths as the template chunks
    (shifted by relative barycentric velocity).
    
    :param observation: Hand a :class:`Observation` object upon template 
        creation. Otherwise, if a pathname is given, the 
        :class:`StellarTemplate_Chunked` is loaded from there.
    :type observation: :class:`Observation` or str
    :param velocity_offset: Upon template creation, hand the velocity-offset 
        between the observation and reference spectrum. Leave as None when 
        loading from file (default).
    :type velocity_offset: float, or None
    :param bary_vel_corr: Upon template creation, hand the barycentric velocity 
        of the observation. Leave as None when loading from file (default).
    :type bary_vel_corr: float, or None
    :param osample: Upon template creation, hand the oversampling factor used 
        in the template creation. Leave as None when loading from file 
        (default).
    :type osample: int, or None
    
    """

    def __init__(self, observation, velocity_offset=None, bary_vel_corr=None, 
                 osample=None):
        
        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                                format='%(message)s')
        
        self.chunks = []
        
        # Load existing template, if first argument is a string
        if type(observation) is str:
            logging.info('Loading template from file: {}'.format(observation))
            
            try:
                with h5py.File(observation, 'r') as h:
                    try:
                        self.orig_filename = h['/orig_filename'][()].decode()
                    except:
                        self.orig_filename = os.path.abspath(observation)
                    self.starname = h['/starname'][()]
                    self.time_start = astropy.time.Time(
                        h['/time_start'][()],
                        format='isot',
                        scale='utc'
                    )
                    self.velocity_offset = h['/velocity_offset'][()]
                    self.bary_vel_corr = h['/bary_vel_corr'][()]
                    self.osample = h['/osample'][()]
                    orders = h['/orders'][()]
                    pix0 = h['/pix0'][()]
                    weight = h['/weight'][()]
                    w0 = h['/w0'][()]
                    w1 = h['/w1'][()]
                    for i in range(len(w1)):
                        flux = h['/chunks/{}/flux'.format(i)][()]
                        wave = h['/chunks/{}/wave'.format(i)][()]
                        pixel = h['/chunks/{}/pixel'.format(i)][()]
                        self.chunks.append(TemplateChunk(flux, wave, pixel, w0[i], w1[i], 
                                            orders[i], pix0[i], weight[i]))
            except Exception as e:
                raise e
        
        # Initialize template from observation
        elif isinstance(observation, Observation) and velocity_offset is not None:
            # FIXME: Handle different velocity offsets (relative to Arcturus atlas etc.)
            self.orig_filename = None
            self.velocity_offset = velocity_offset
            self.starname = observation.star.name
            self.time_start = observation.time_start  # Start time of the observation <astropy.time.Time>
            self.bary_vel_corr = bary_vel_corr
            if self.bary_vel_corr is None:
                self.bary_vel_corr = observation[0].bary_vel_corr
            self.osample = osample
            
        else:
            raise ValueError("""Invalid input arguments. Usage:
                1) StellarTemplate_Chunked('filename.h5')
                2) StellarTemplate_Chunked(template_observation, velocity_offset, bary_vel_corr, osample)""")
    
    def append(self, templatechunk):
        """Append chunks to the list
        
        :param templatechunk: The chunk to append.
        :type: :class:`TemplateChunk`
        
        """
        self.chunks.append(templatechunk)
    
    @property
    def orders(self):
        """A ndarray containing the order number of each :class:`TemplateChunk`
        
        :return: The order numbers.
        :rtype: ndarray[nr_chunks]
        """
        return np.array([chunk.order for chunk in self.chunks], dtype='int')
    
    @property
    def orders_unique(self):
        """A ndarray containing the orders covered by the template
        
        :return: The unique order numbers.
        :rtype: ndarray[nr_orders]
        """
        return np.sort(np.unique(self.orders))
    
    @property
    def pix0(self):
        """A ndarray containing the first pixel index of each
        :class:`TemplateChunk`
        
        :return: The pixel 0 values.
        :rtype: ndarray[nr_chunks]
        """
        return np.array([chunk.pix0 for chunk in self.chunks])
    
    @property
    def weight(self):
        """A ndarray containing the weight of each :class:`TemplateChunk`
        
        :return: The weights of the chunks.
        :rtype: ndarray[nr_chunks]
        """
        return np.array([chunk.weight for chunk in self.chunks])
    
    @property
    def w0(self):
        """A ndarray containing the wavelength intercepts of each
        :class:`TemplateChunk`
        
        :return: The wavelength intercepts of the chunks.
        :rtype: ndarray[nr_chunks]
        """
        return np.array([chunk.w0 for chunk in self.chunks])
    
    @property
    def w1(self):
        """A ndarray containing the wavelength slopes of each 
        :class:`TemplateChunk`
        
        :return: The wavelength slopes of the chunks.
        :rtype: ndarray[nr_chunks]
        """
        return np.array([chunk.w1 for chunk in self.chunks])
    
    def get_order_indices(self, order) -> list:
        """Return a list of indices of the chunks within a given order
        
        :param order: The order of interest.
        :type order: int
        
        :return: The chunk indices within the order.
        :rtype: list
        """
        return [i for i in range(len(self)) if self.chunks[i].order == order]
    
    def check_wavelength_range(self, wave_start, wave_stop):
        """Find the chunk with the best coverage of the wavelength range 
        defined by wave_start and wave_stop.
        If no such chunk is found, raise a NoDataError. Otherwise return the
        chunk index and coverage as a tuple.
        
        :param wave_start: Starting wavelength.
        :type wave_start: float
        :param wave_stop: Stopping wavelength.
        :type wave_stop: float
        
        :return: The index of the chunk which best covers the wavelength
            range.
        :rtype: int
        :return: A value between 0.0 and 1.0, telling how big a fraction of 
            the wavelength range [wave_start:wave_stop] is covered by data.
        :rtype: float
        
        """
        selected_chunk = None
        best_coverage = 0.
        # Loop through chunks and look for requested wavelengths
        for i in range(len(self)):
            try:
                coverage = self.chunks[i].check_wavelength_range(wave_start, wave_stop)
                # Select chunk if coverage better than previous chunk
                if coverage > best_coverage:
                    best_coverage = coverage
                    selected_chunk = i
            except NoDataError:
                pass
        if selected_chunk is None:
            raise NoDataError(
                'Could not find wavelength range {}-{} Å'.format(
                    wave_start, wave_stop
                ))
        return selected_chunk, best_coverage
    
    def get_wavelength_range(self, wave_start, wave_stop, require=None):
        """Loop through template chunks and search for a chunk with these
        wavelengths. Return a wavelength inverval of that chunk, defined by 
        wave_start and wave_stop. If multiple results, return the one with best 
        coverage.
        
        :param wave_start: Starting wavelength.
        :type wave_start: float
        :param wave_stop: Stopping wavelength.
        :type wave_stop: float
        :param require: If require='full', ensure that the full interval is 
            covered (one pixel outside in each end).
        :type require: str, or None
        
        :return: The chunk which best covers the selected wavelength range.
        :rtype: :class:`TemplateChunk`
        
        """
        selected_chunk, best_coverage = \
            self.check_wavelength_range(wave_start, wave_stop)
        return self.chunks[selected_chunk].get_wavelength_range(
            wave_start, wave_stop, require=require)
        

    def save(self, filename):
        """Save as HDF5 file (.h5)
        
        :param filename: The pathname of the directory where to save the
            file, or the filename itself.
        :type filename: str
        
        """
        
        logging.info('Saving deconvolved template to {}'.format(filename))
        
        if os.path.isdir(filename):
            date = self.time_start.datetime.strftime("%Y%m%d")
            filename = os.path.join(filename, self.starname + '_' + date + '.h5')
            logging.warning('Supplied filename was a directory.')
            logging.warning('New filename: {}'.format(filename))
            
        with h5py.File(filename, 'w') as h:
            h.create_dataset('/orig_filename', data=os.path.abspath(filename))
            h.create_dataset('/starname', data=self.starname)
            h.create_dataset('/time_start', data=self.time_start.isot)
            h.create_dataset('/velocity_offset', data=self.velocity_offset)
            h.create_dataset('/bary_vel_corr', data=self.bary_vel_corr)
            h.create_dataset('/osample', data=self.osample)
            h.create_dataset('/orders', data=self.orders)
            h.create_dataset('/pix0', data=self.pix0)
            h.create_dataset('/weight', data=self.weight)
            h.create_dataset('/w0', data=self.w0)
            h.create_dataset('/w1', data=self.w1)
            
            for i in range(len(self)):
                grp = '/chunks/{}/'.format(i)
                h.create_dataset(grp + 'flux', data=self.chunks[i].flux)
                h.create_dataset(grp + 'wave', data=self.chunks[i].wave)
                h.create_dataset(grp + 'pixel', data=self.chunks[i].pixel)
    
    def __len__(self):
        return len(self.chunks)
    
    def __getitem__(self, index):
        try:
            return self.chunks[index]
        except:
            raise NoDataError('No data for chunk {}!'.format(index))

    def __setitem__(self, index, templatechunk):
        self.chunks[index] = templatechunk
    
    def __str__(self):
        return '<StellarTemplate_Chunked of {} ({} chunks, {} orders)>'.format(
            self.starname,
            len(self),
            len(self.orders_unique)
        )
