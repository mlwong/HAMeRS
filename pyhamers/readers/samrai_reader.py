"""
Module for reading and handling samrai data.
"""

import copy
import h5py
import numpy
import re

from pyhamers.upsampling import Lagrange_upsampler

from pyhamers.readers.base_reader import BaseReader

class SamraiDataReader(BaseReader):
    """
    Class to read samrai data.
    """

    def __init__(self, data_directory_path, periodic_dimensions = (False, False, False), \
                 upsampling_method = 'constant', processor_zero_padding_length = 5, data_order = 'F'):
        """
        Constructor of the class.
        The current time step of the class is set to the first time step in dump file.
        """

        self._data_directory_path = data_directory_path

        # Get the full paths to data at different time steps.

        full_dumps_path = data_directory_path + '/' + 'dumps.visit'
        dumps_f = open(full_dumps_path)

        self._full_viz_folder_paths = {}
        self._steps = []

        for line in dumps_f.readlines():
            sub_folder = line
            sub_folder = re.sub('/summary.samrai\n', '', sub_folder)
            step = int(re.sub('visit_dump.', '', sub_folder))
            self._steps.append(step)
            self._full_viz_folder_paths[step] = data_directory_path + '/' + sub_folder + '/'

        # Set the data order.

        if data_order == 'C':
            self._data_order = 'C'
        elif data_order == 'F':
            self._data_order = 'F'
        else:
            raise RuntimeError("Unknown data order '" + data_order + "'!")

        # Step current time step to be the first time step in dump file and read the summary file
        # at that time step.

        self._step = self._steps[0]

        self._basic_info = {}
        self._readSummary(self._step)

        # Set the periodic dimensions.

        dim = self._basic_info['dim']

        if dim == 1:
            if len(periodic_dimensions) < 1:
                raise RuntimeError('Dimension of periodic_dimensions is not correct!')

            self._periodic_dimensions = (periodic_dimensions[0],)

        elif dim == 2:
            if len(periodic_dimensions) < 2:
                raise RuntimeError('Dimension of periodic_dimensions is not correct!')

            self._periodic_dimensions = (periodic_dimensions[0], periodic_dimensions[1])

        elif dim == 3:
            if len(periodic_dimensions) < 3:
                raise RuntimeError('Dimension of periodic_dimensions is not correct!')

            self._periodic_dimensions = (periodic_dimensions[0], periodic_dimensions[1], periodic_dimensions[2])

        # Set up the upsampling classes.

        self._upsampler = []
        if upsampling_method == 'constant':
            self._upsampler = Lagrange_upsampler.LagrangeUpsampler('constant', data_order=self._data_order)
        elif upsampling_method == 'second_order_Lagrange':
            self._upsampler = Lagrange_upsampler.LagrangeUpsampler('second_order', data_order=self._data_order)
        elif upsampling_method == 'fourth_order_Lagrange':
            self._upsampler = Lagrange_upsampler.LagrangeUpsampler('fourth_order', data_order=self._data_order)
        elif upsampling_method == 'sixth_order_Lagrange':
            self._upsampler = Lagrange_upsampler.LagrangeUpsampler('sixth_order', data_order=self._data_order)
        else:
            raise RuntimeError("Unknown method '" + upsampling_method + "' for upsampling!")

        self._upsampler_constant = Lagrange_upsampler.LagrangeUpsampler('constant', data_order=self._data_order)

        # Set the length of zero-padding for the processor number.

        if not isinstance(processor_zero_padding_length, int):
            raise RunTimeError("Length of zero-padding for processor number is not integer!")

        if processor_zero_padding_length < 1:
            raise RunTimeError("Length of zero-padding for processor number is smaller than 1!")

        self._processor_zero_padding_length = processor_zero_padding_length

        # Initialize subdomain.

        self._domain_size = self.getRefinedDomainSize()

        self._lo_subdomain = ()
        self._hi_subdomain = ()

        if dim == 1:
            self._lo_subdomain = (0,)
            self._hi_subdomain = (self._domain_size[0] - 1)

        elif dim == 2:
            self._lo_subdomain = (0, 0)
            self._hi_subdomain = (self._domain_size[0] - 1, self._domain_size[1] - 1)

        elif dim == 3:
            self._lo_subdomain = (0, 0, 0)
            self._hi_subdomain = (self._domain_size[0] - 1, self._domain_size[1] - 1, self._domain_size[2] - 1)

        # Initialize other containers.

        self._data_loaded = False
        self._data = {}


    @property
    def dimension(self):
        """
        Return the dimension of the domain.
        """

        return self._basic_info['dim']


    @property
    def periodic_dimensions(self):
        """
        Return a tuple indicating if data is periodic in each dimension.
        """

        return self._periodic_dimensions


    @property
    def time(self):
        """
        Return the simulation time at current time step.
        """

        return self._basic_info['t']


    @property
    def data_order(self):
        """
        Return the data order.
        """

        return self._data_order


    @property
    def steps(self):
        """
        Return all of the steps.
        """

        return self._steps


    def getDataDirectoryPath(self):
        """
        Get the stored absolute path of the data directory.
        """

        return self._data_directory_path


    def setStep(self, step):
        """
        Update the metadata from the summary file in the data directory at new time step and
        change the current time step to the new time step.
        """

        assert (step in self._steps), "Step to read in is not available in the dataset."

        self._step = step
        self._readSummary(step)


    def getStep(self):
        """
        Return the time step that is currently set.
        """

        return self._step


    step = property(getStep, setStep)


    def _readSummary(self, step):
        """
        Get the basic information, patch extents and path map from the summary file.
        """

        # Open the summary file.

        summary_file_path = self._full_viz_folder_paths[step] + '/' + 'summary.samrai'
        f_summary = h5py.File(summary_file_path, 'r')

        # Clear metadata in self._basic_info.

        self._basic_info = {}

        # Get the basic information.

        basic_info = f_summary['BASIC_INFO']

        # Get the number of file clusters.
        self._basic_info['num_file_clusters'] = basic_info['number_file_clusters'][()][0]

        # Get the time and dimension of the data.

        self._basic_info['t'] = basic_info['time'][()][0]
        self._basic_info['t_dump'] = basic_info['time_of_dump'][()][0]
        self._basic_info['n'] = basic_info['time_step_number'][()][0]
        self._basic_info['dim'] = basic_info['number_dimensions_of_problem'][()][0]

        # Get and check the grid type.

        self._basic_info['grid_type'] = basic_info['grid_type'][()][0].decode()

        if numpy.char.strip(self._basic_info['grid_type']) != 'CARTESIAN':
            raise RuntimeError("Grid type other than 'CARTESIAN' not supported!")

        # Get the number of levels and number of patches at different levels.

        self._basic_info['num_levels'] = basic_info['number_levels'][()][0]
        self._basic_info['num_patches'] = basic_info['number_patches_at_level'][()]
        self._basic_info['num_global_patches'] = basic_info['number_global_patches'][()][0]

        # Get the ratios to coarser levels at different levels.

        self._basic_info['ratios_to_coarser_levels'] = basic_info['ratios_to_coarser_levels'][()]

        # Get the variable names, number of variables and number of components in each variable.

        self._basic_info['var_names'] = basic_info['var_names'][()]
        var_names_decoded = []
        for i in range(len(self._basic_info['var_names'])):
            var_names_decoded.append(self._basic_info['var_names'][i].decode())
        self._basic_info['var_names'] = numpy.array(var_names_decoded)
        self._basic_info['num_variables'] = basic_info['number_visit_variables'][()][0]
        self._basic_info['num_var_components'] = basic_info['var_number_components'][()]

        # Get the geometry and check the dimension.

        self._basic_info['x_lo'] = basic_info['XLO'][()]
        self._basic_info['dx'] = basic_info['dx'][()]

        extents = f_summary['extents']

        # Get the patch extents.

        self._patch_extents = extents['patch_extents'][()]

        # Geth the patch map.

        self._patch_map =  extents['patch_map'][()]

        # Set the flag for loading summary file to be true.

        self._summary_loaded = True

        # Close the summary file.

        f_summary.close()


    def getBasicInfo(self):
        """
        Return the loaded basic information.
        """

        if not self._summary_loaded:
            raise RuntimeError('The summary file is not read yet!')

        return self._basic_info


    def getPatchExtents(self):
        """
        Return the loaded patch extents.
        """

        if not self._summary_loaded:
            raise RuntimeError('The summary file is not read yet!')

        return self._patch_extents


    def getPatchMap(self):
        """
        Return the loaded patch map.
        """

        if not self._summary_loaded:
            raise RuntimeError('The summary file is not read yet!')

        return self._patch_map


    def getData(self, var_name):
        """
        Return the loaded data.
        """

        if not self._data_loaded:
            raise RuntimeError('No data is read yet!')

        return self._data[var_name]


    def clearData(self):
        """
        Clear any loaded data.
        """

        self._data.clear()
        self._data_loaded = False


    def getDomainSizeAtOneLevel(self, \
            level_num):
        """
        Get the domain size at one particular level.
        """

        # Read the summary file if it is not yet read.

        if not self._summary_loaded:
            self.readSummary()

        dim = self._basic_info['dim']
        num_levels = self._basic_info['num_levels']
        num_patches = self._basic_info['num_patches']

        # Check whether the required level is valid.

        if level_num >= num_levels:
            raise RuntimeError('Level number is greater than number of levels!')
        elif level_num < 0:
            raise RuntimeError('Level number is negative!')

        # Get the domain shape at the level.

        num_patches_level = num_patches[level_num]

        patch_level_start_idx = 0
        for level_idx in range(0, level_num):
            patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]

        lo_level = self._patch_extents[patch_level_start_idx][0]
        hi_level = self._patch_extents[patch_level_start_idx][1]

        for global_patch_idx in range(patch_level_start_idx + 1, patch_level_start_idx + num_patches_level):
            lo_level = numpy.minimum(lo_level, self._patch_extents[global_patch_idx][0])
            hi_level = numpy.maximum(hi_level, self._patch_extents[global_patch_idx][1])

        domain_shape = hi_level[0:dim] - lo_level[0:dim] + numpy.ones(dim, dtype = numpy.int)

        return tuple(domain_shape)


    def getRefinedDomainSize(self):
        """
        Get the full domain size refined to the highest level.
        """

        # Read the summary file if it is not yet read.

        if not self._summary_loaded:
            self.readSummary()

        dim = self._basic_info['dim']
        num_levels = self._basic_info['num_levels']
        num_patches = self._basic_info['num_patches']

        # Get the ratio from the coarest level to the finest level.
        ratios_to_coarser_levels = self._basic_info['ratios_to_coarser_levels']
        ratio_of_coarest_to_finest = numpy.ones(dim, dtype = ratios_to_coarser_levels.dtype)

        for level_idx in range(1, num_levels):
            ratio_of_coarest_to_finest = numpy.multiply(ratios_to_coarser_levels[level_idx][0:dim], ratio_of_coarest_to_finest)

        # Get the refined domain shape at the root level.

        num_patches_root_level = num_patches[0]

        lo_root_level = self._patch_extents[0][0]
        hi_root_level = self._patch_extents[0][1]

        for patch_idx in range(1, num_patches_root_level):
            lo_root_level = numpy.minimum(lo_root_level, self._patch_extents[patch_idx][0])
            hi_root_level = numpy.maximum(hi_root_level, self._patch_extents[patch_idx][1])

        domain_shape = hi_root_level[0:dim] - lo_root_level[0:dim] + numpy.ones(dim, dtype = lo_root_level.dtype)
        domain_shape = numpy.multiply(domain_shape, ratio_of_coarest_to_finest)

        return tuple(domain_shape)


    @property
    def domain_size(self):
        """
        Return a tuple containing the full domain size of this dataset.
        """

        return self._domain_size


    def setSubDomain(self, lo_and_hi):
        """
        Set the sub-domain for reading coordinates and data in a subdomain.
        """

        try:
            lo, hi = lo_and_hi
        except ValueError:
            raise ValueError("Pass an iterable with two items!")

        dim = self._basic_info['dim']

        for i in range(dim):
            if lo[i] < 0 or lo[i] > self._domain_size[i]:
                raise ValueError('Invalid indices in sub-domain. Cannot be < 0 or > domain size!')
            if hi[i] < 0 or hi[i] > self._domain_size[i]:
                raise ValueError('Invalid indices in sub-domain. Cannot be < 0 or > domain size!')
            if hi[i] < lo[i]:
                raise ValueError('Invalid indices in sub-domain. Upper bound cannot be smaller than lower bound!')

        if self._data_loaded == True:
            self.clearData()

        if dim == 1:
            if len(lo) < 1 or len(hi) < 1:
                raise RuntimeError('Dimension of lo or hi is not correct!')

            if lo[0] > hi[0]:
                raise RuntimeError('lo is greater than hi!')

            self._lo_subdomain = (lo[0],)
            self._hi_subdomain = (hi[0],)

        elif dim == 2:
            if len(lo) < 2 or len(hi) < 2:
                raise RuntimeError('Dimension of lo or hi is not correct!')

            if lo[0] > hi[0] or lo[1] > hi[1]:
                raise RuntimeError('lo is greater than hi!')

            self._lo_subdomain = (lo[0], lo[1])
            self._hi_subdomain = (hi[0], hi[1])

        elif dim == 3:
            if len(lo) < 3 or len(hi) < 3:
                raise RuntimeError('Dimension of lo or hi is not correct!')

            if lo[0] > hi[0] or lo[1] > hi[1] or lo[2] > hi[2]:
                raise RuntimeError('lo is greater than hi!')

            self._lo_subdomain = (lo[0], lo[1], lo[2])
            self._hi_subdomain = (hi[0], hi[1], hi[2])


    def getSubDomain(self):
        """
        Return two tuples containing the sub-domain used in this reader
        as a lower bound (lo) and upper bound (hi).
        """

        return self._lo_subdomain, self._hi_subdomain


    sub_domain = property(getSubDomain, setSubDomain)


    def getCoordinatesAtOneLevel(self, level_num):
        """
        Get coordinates at one particular level.
        """

        # Get the dimension of the problem, number of levels and number of patches.

        dim = self._basic_info['dim']
        num_levels = self._basic_info['num_levels']
        num_patches = self._basic_info['num_patches']

        # Get the domain shape.

        domain_shape = self.getDomainSizeAtOneLevel(level_num)

        # Get the number of patches at this level and grid spacings at different levels.

        num_patches_level = num_patches[level_num]
        dx = self._basic_info['dx']

        # Get the coordinates of the corners at this level.

        patch_level_start_idx = 0
        for level_idx in range(0, level_num):
            patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]

        x_lo_level = self._patch_extents[patch_level_start_idx][2]
        x_hi_level = self._patch_extents[patch_level_start_idx][3]

        for global_patch_idx in range(patch_level_start_idx + 1, patch_level_start_idx + num_patches_level):
            x_lo_level = numpy.minimum(x_lo_level, self._patch_extents[global_patch_idx][2])
            x_hi_level = numpy.maximum(x_hi_level, self._patch_extents[global_patch_idx][3])

        # Compute the coordinates at this level.

        x_coords = []
        y_coords = []
        z_coords = []

        if dim == 1:
            x_coords = numpy.linspace(x_lo_level[0] + 0.5*dx[level_num][0], \
                x_hi_level[0] - 0.5*dx[level_num][0], \
                num = domain_shape[0])

            return x_coords

        elif dim == 2:
            x_coords = numpy.linspace(x_lo_level[0] + 0.5*dx[level_num][0], \
                x_hi_level[0] - 0.5*dx[level_num][0], \
                num = domain_shape[0])
            y_coords = numpy.linspace(x_lo_level[1] + 0.5*dx[level_num][1], \
                x_hi_level[1] - 0.5*dx[level_num][1], \
                num = domain_shape[1])

            return x_coords, y_coords

        elif dim == 3:
            x_coords = numpy.linspace(x_lo_level[0] + 0.5*dx[level_num][0], \
                x_hi_level[0] - 0.5*dx[level_num][0], \
                num = domain_shape[0])
            y_coords = numpy.linspace(x_lo_level[1] + 0.5*dx[level_num][1], \
                x_hi_level[1] - 0.5*dx[level_num][1], \
                num = domain_shape[1])
            z_coords = numpy.linspace(x_lo_level[2] + 0.5*dx[level_num][2], \
                x_hi_level[2] - 0.5*dx[level_num][2], \
                num = domain_shape[2])

            return x_coords, y_coords, z_coords


    def readDataAtOneLevel(self, \
            var_names,
            level_num):
        """
        Read data at one particular level.
        """

        # Get the number of file clusters.

        num_file_clusters = self._basic_info['num_file_clusters']

        # Get the dimension of the problem, number of levels and number of patches.

        dim = self._basic_info['dim']
        num_levels = self._basic_info['num_levels']
        num_patches = self._basic_info['num_patches']

        # Check whether the required level is valid.

        if level_num >= num_levels:
            raise RuntimeError('Level number is greater than number of levels!')
        elif level_num < 0:
            raise RuntimeError('Level number is negative!')

        # Get the variable names.

        var_num_components = {}
        var_component_names = {}

        for var_name in var_names:
            var_idx = numpy.where(self._basic_info['var_names'] == var_name)[0][0]
            var_num_components[var_name] = self._basic_info['num_var_components'][var_idx]

            var_component_names[var_name] = [None]*var_num_components[var_name]
            if var_num_components[var_name] == 1:
                var_component_names[var_name] = [var_name]
            else:
                for component_idx in range(0, var_num_components[var_name]):
                    var_component_names[var_name][component_idx] = var_name + '.' + str(component_idx).zfill(2)

        # Get the domain shape.

        domain_shape = self.getDomainSizeAtOneLevel(level_num)

        # Get the index of the lower corner at this level.

        num_patches_level = num_patches[level_num]

        patch_level_start_idx = 0
        for level_idx in range(0, level_num):
            patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]

        lo_level = self._patch_extents[patch_level_start_idx][0]

        for global_patch_idx in range(patch_level_start_idx + 1, patch_level_start_idx + num_patches_level):
            lo_level = numpy.minimum(lo_level, self._patch_extents[global_patch_idx][0])

        # Initialize container to store the data. The elements in the container are initialized as NAN values.

        for var_name in var_names:
            if self._data_order == 'C':
                data_shape = numpy.insert(domain_shape, 0, var_num_components[var_name])
            else:
                data_shape = numpy.append(domain_shape, var_num_components[var_name])

            self._data[var_name] = numpy.empty(data_shape, dtype = numpy.float64, order = self._data_order)
            self._data[var_name][:] = numpy.NAN

        # Get the data from all patches at the specified level.

        if dim == 1:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(self._processor_zero_padding_length) + '.samrai'
                full_path = self._full_viz_folder_paths[self._step] + '/' + file_name
                f_input = h5py.File(full_path, 'r')

                file_cluster = f_input['processor.' + str(process_idx).zfill(self._processor_zero_padding_length)]
                file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]

                for var_name in var_names:
                    for patch_key in file_cluster_level:
                        patch_idx = int(patch_key.replace('patch.', ''))
                        global_patch_idx = patch_level_start_idx + patch_idx
                        file_cluster_patch = file_cluster_level[patch_key]

                        lo_patch = patch_extents[global_patch_idx][0] - lo_level
                        hi_patch = patch_extents[global_patch_idx][1] - lo_level
                        patch_shape = hi_patch[0] - lo_patch[0] + 1

                        x_start_idx = lo_patch[0]
                        x_end_idx = hi_patch[0] + 1

                        for component_idx in range(0, var_num_components[var_name]):
                            if self._data_order == 'C':
                                self._data[var_name][component_idx, x_start_idx:x_end_idx] = \
                                    file_cluster_patch[var_component_names[var_name][component_idx]][()]
                            else:
                                self._data[var_name][x_start_idx:x_end_idx, component_idx] = \
                                    file_cluster_patch[var_component_names[var_name][component_idx]][()]

                f_input.close()

        elif dim == 2:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(self._processor_zero_padding_length) + '.samrai'
                full_path = self._full_viz_folder_paths[self._step] + '/' + file_name
                f_input = h5py.File(full_path, 'r')

                file_cluster = f_input['processor.' + str(process_idx).zfill(self._processor_zero_padding_length)]
                file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]

                for var_name in var_names:
                    for patch_key in file_cluster_level:
                        patch_idx = int(patch_key.replace('patch.', ''))
                        global_patch_idx = patch_level_start_idx + patch_idx
                        file_cluster_patch = file_cluster_level[patch_key]

                        lo_patch = self._patch_extents[global_patch_idx][0] - lo_level
                        hi_patch = self._patch_extents[global_patch_idx][1] - lo_level
                        patch_shape = hi_patch[0:2] - lo_patch[0:2] + numpy.ones(2, dtype = numpy.int)

                        x_start_idx = lo_patch[0]
                        x_end_idx = hi_patch[0] + 1

                        y_start_idx = lo_patch[1]
                        y_end_idx = hi_patch[1] + 1

                        for component_idx in range(0, var_num_components[var_name]):
                            if self._data_order == 'C':
                                self._data[var_name][component_idx, x_start_idx:x_end_idx, y_start_idx:y_end_idx] = \
                                    file_cluster_patch[var_component_names[var_name][component_idx]][()].reshape( \
                                        patch_shape, order = 'F')
                            else:
                                self._data[var_name][x_start_idx:x_end_idx, y_start_idx:y_end_idx, component_idx] = \
                                    file_cluster_patch[var_component_names[var_name][component_idx]][()].reshape( \
                                        patch_shape, order = 'F')

                f_input.close()

        elif dim == 3:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(self._processor_zero_padding_length) + '.samrai'
                full_path = self._full_viz_folder_paths[self._step] + '/' + file_name
                f_input = h5py.File(full_path, 'r')

                file_cluster = f_input['processor.' + str(process_idx).zfill(self._processor_zero_padding_length)]
                file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]

                for var_name in var_names:
                    for patch_key in file_cluster_level:
                        patch_idx = int(patch_key.replace('patch.', ''))
                        global_patch_idx = patch_level_start_idx + patch_idx
                        file_cluster_patch = file_cluster_level[patch_key]

                        lo_patch = self._patch_extents[global_patch_idx][0] - lo_level
                        hi_patch = self._patch_extents[global_patch_idx][1] - lo_level
                        patch_shape = hi_patch[0:3] - lo_patch[0:3] + numpy.ones(3, dtype = numpy.int)

                        x_start_idx = lo_patch[0]
                        x_end_idx = hi_patch[0] + 1

                        y_start_idx = lo_patch[1]
                        y_end_idx = hi_patch[1] + 1

                        z_start_idx = lo_patch[2]
                        z_end_idx = hi_patch[2] + 1

                        for component_idx in range(0, var_num_components[var_name]):
                            if self._data_order == 'C':
                                self._data[var_name] \
                                    [component_idx, x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx] = \
                                        file_cluster_patch[var_component_names[var_name][component_idx]][()].reshape( \
                                            patch_shape, order = 'F')
                            else:
                                self._data[var_name] \
                                    [x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx, component_idx] = \
                                        file_cluster_patch[var_component_names[var_name][component_idx]][()].reshape( \
                                            patch_shape, order = 'F')

                f_input.close()

        else:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')

        self._data_loaded = True


    def getCombinedCoordinatesInSubdomainFromAllLevels(self, num_ghosts = None):
        """
        Get coordinates in a sub-domain from all levels, refine the coordinates to the finest level
        and combine the coordinates from different levels.
        """

        lo_subdomain = numpy.asarray(self._lo_subdomain)
        hi_subdomain = numpy.asarray(self._hi_subdomain)

        # Get the dimension of the problem, number of levels and number of patches.

        dim = self._basic_info['dim']
        num_levels = self._basic_info['num_levels']
        num_patches = self._basic_info['num_patches']
        num_patches_root_level = num_patches[0]

        if num_ghosts is None:
            if dim == 1:
                num_ghosts = (0,)

            elif dim == 2:
                num_ghosts = (0, 0)

            elif dim == 3:
                num_ghosts = (0, 0, 0)

        else:
            # Check the length of num_ghosts is correct.
            if dim == 1:
                if len(num_ghosts) < 1:
                    raise RuntimeError('Dimension of num_ghosts is not correct!')

            elif dim == 2:
                if len(num_ghosts) < 2:
                    raise RuntimeError('Dimension of num_ghosts is not correct!')

            elif dim == 3:
                if len(num_ghosts) < 3:
                    raise RuntimeError('Dimension of num_ghosts is not correct!')

        num_ghosts = numpy.asarray(num_ghosts[0:dim])

        # Get the refinement ratios from different levels to finest level.

        ratios_to_coarser_levels = self._basic_info['ratios_to_coarser_levels']
        ratios_to_finest_level = numpy.empty(ratios_to_coarser_levels.shape, dtype = ratios_to_coarser_levels.dtype)
        ratios_to_finest_level[num_levels - 1] = -ratios_to_coarser_levels[0]
        for level_idx in range(num_levels - 2, -1, -1):
            ratios_to_finest_level[level_idx] = numpy.multiply(ratios_to_coarser_levels[level_idx + 1], \
                                                ratios_to_finest_level[level_idx + 1])

        # Get the lower and upper indices of the domain.

        lo_root_level = self._patch_extents[0][0]
        hi_root_level = self._patch_extents[0][1]

        for patch_idx in range(1, num_patches_root_level):
            lo_root_level = numpy.minimum(lo_root_level, self._patch_extents[patch_idx][0])
            hi_root_level = numpy.maximum(hi_root_level, self._patch_extents[patch_idx][1])

        # Refine the the lower and upper indices of the domain to the highest level.

        lo_root_level_refined = numpy.multiply(lo_root_level[0:dim], ratios_to_finest_level[0][0:dim])
        hi_root_level_refined = numpy.multiply(hi_root_level[0:dim] + numpy.ones(dim, dtype = numpy.int), \
            ratios_to_finest_level[0][0:dim]) \
            - numpy.ones(dim, dtype = numpy.int)

        # Compute the shape of the domain refined to the highest level.

        domain_shape = hi_root_level_refined[0:dim] - lo_root_level_refined[0:dim] \
            + numpy.ones(dim, dtype = lo_root_level_refined.dtype)

        # Compute the shape of the domain refined to the highest level with ghost cells.

        domain_shape_ghosts = domain_shape + 2*num_ghosts

        # Get the grid spacings at different levels.

        dx = self._basic_info['dx']

        # Get the coordinates of the corners at this level.

        x_lo_root_level = self._patch_extents[0][2]
        x_hi_root_level = self._patch_extents[0][3]

        for patch_idx in range(1, num_patches_root_level):
            x_lo_root_level = numpy.minimum(x_lo_root_level, self._patch_extents[patch_idx][2])
            x_hi_root_level = numpy.maximum(x_hi_root_level, self._patch_extents[patch_idx][3])

        # Compute the coordinates of the full domain refined to the finest level.

        x_coords = []
        y_coords = []
        z_coords = []

        # Include the ghost cells in the sub-domain.

        lo_subdomain = lo_subdomain[0:dim] - num_ghosts[0:dim]
        hi_subdomain = hi_subdomain[0:dim] + num_ghosts[0:dim]

        if dim == 1:
            x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_hi_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = domain_shape_ghosts[0])

            x_coords = x_coords[lo_subdomain[0] + num_ghosts[0]:hi_subdomain[0] + 1 + num_ghosts[0]]

            return x_coords

        elif dim == 2:
            x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_hi_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = domain_shape_ghosts[0])
            y_coords = numpy.linspace(x_lo_root_level[1] + (0.5 - num_ghosts[1])*dx[-1][1],
                x_hi_root_level[1] + (num_ghosts[1] - 0.5)*dx[-1][1], \
                num = domain_shape_ghosts[1])

            x_coords = x_coords[lo_subdomain[0] + num_ghosts[0]:hi_subdomain[0] + 1 + num_ghosts[0]]
            y_coords = y_coords[lo_subdomain[1] + num_ghosts[1]:hi_subdomain[1] + 1 + num_ghosts[1]]

            return x_coords, y_coords

        elif dim == 3:
            x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_hi_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = domain_shape_ghosts[0])
            y_coords = numpy.linspace(x_lo_root_level[1] + (0.5 - num_ghosts[1])*dx[-1][1], \
                x_hi_root_level[1] + (num_ghosts[1] - 0.5)*dx[-1][1], \
                num = domain_shape_ghosts[1])
            z_coords = numpy.linspace(x_lo_root_level[2] + (0.5 - num_ghosts[2])*dx[-1][2], \
                x_hi_root_level[2] + (num_ghosts[2] - 0.5)*dx[-1][2], \
                num = domain_shape_ghosts[2])

            x_coords = x_coords[lo_subdomain[0] + num_ghosts[0]:hi_subdomain[0] + 1 + num_ghosts[0]]
            y_coords = y_coords[lo_subdomain[1] + num_ghosts[1]:hi_subdomain[1] + 1 + num_ghosts[1]]
            z_coords = z_coords[lo_subdomain[2] + num_ghosts[2]:hi_subdomain[2] + 1 + num_ghosts[2]]

            return x_coords, y_coords, z_coords


    def readCombinedDataInSubdomainFromAllLevels(self, \
            var_names, \
            num_ghosts = None):
        """
        Read data in a sub-domain from all levels, refine the data to the finest level
        and combine the data from different levels.
        """

        lo_subdomain = numpy.asarray(self._lo_subdomain)
        hi_subdomain = numpy.asarray(self._hi_subdomain)

        periodic_dimensions = self._periodic_dimensions

        # Get the number of file clusters.

        num_file_clusters = self._basic_info['num_file_clusters']

        # Get the dimension of the problem, number of levels and number of patches.

        dim = self._basic_info['dim']
        num_levels = self._basic_info['num_levels']
        num_patches = self._basic_info['num_patches']
        num_patches_root_level = num_patches[0]

        # Get the number of ghost cells.

        if num_ghosts is None:
            if dim == 1:
                num_ghosts = (0,)

            elif dim == 2:
                num_ghosts = (0, 0)

            elif dim == 3:
                num_ghosts = (0, 0, 0)

        else:
            # Check the length of num_ghosts is correct.

            if dim == 1:
                if len(num_ghosts) < 1:
                    raise RuntimeError('Dimension of num_ghosts is not correct!')

            elif dim == 2:
                if len(num_ghosts) < 2:
                    raise RuntimeError('Dimension of num_ghosts is not correct!')

            elif dim == 3:
                if len(num_ghosts) < 3:
                    raise RuntimeError('Dimension of num_ghosts is not correct!')

        num_ghosts = numpy.asarray(num_ghosts[0:dim])

        # Get the variable names.

        var_num_components = {}
        var_component_names = {}

        for var_name in var_names:
            var_idx = numpy.where(self._basic_info['var_names'] == var_name)[0][0]
            var_num_components[var_name] = self._basic_info['num_var_components'][var_idx]

            var_component_names[var_name] = [None]*var_num_components[var_name]
            if var_num_components[var_name] == 1:
                var_component_names[var_name] = [var_name]
            else:
                for component_idx in range(0, var_num_components[var_name]):
                    var_component_names[var_name][component_idx] = var_name + '.' + str(component_idx).zfill(2)

        # Get the refinement ratios from different levels to finest level.

        ratios_to_coarser_levels = self._basic_info['ratios_to_coarser_levels']
        ratios_to_finest_level = numpy.empty(ratios_to_coarser_levels.shape, dtype = ratios_to_coarser_levels.dtype)
        ratios_to_finest_level[num_levels - 1] = -ratios_to_coarser_levels[0]
        for level_idx in range(num_levels - 2, -1, -1):
            ratios_to_finest_level[level_idx] = numpy.multiply(ratios_to_coarser_levels[level_idx + 1], \
                                                ratios_to_finest_level[level_idx + 1])

        # Get the lower and upper indices of the domain.

        lo_root_level = self._patch_extents[0][0]
        hi_root_level = self._patch_extents[0][1]

        for patch_idx in range(1, num_patches_root_level):
            lo_root_level = numpy.minimum(lo_root_level, self._patch_extents[patch_idx][0])
            hi_root_level = numpy.maximum(hi_root_level, self._patch_extents[patch_idx][1])

        # Refine the the lower and upper indices of the domain to the highest level.

        lo_root_level_refined = []
        hi_root_level_refined = []

        if num_levels == 1:
            lo_root_level_refined = lo_root_level[0:dim]
            hi_root_level_refined = hi_root_level[0:dim]

        else:
            lo_root_level_refined = numpy.multiply(lo_root_level[0:dim], ratios_to_finest_level[0][0:dim])
            hi_root_level_refined = numpy.multiply(hi_root_level[0:dim] + numpy.ones(dim, dtype = numpy.int), \
                ratios_to_finest_level[0][0:dim]) \
                - numpy.ones(dim, dtype = numpy.int)

        # Compute the shape of the domain refined to the highest level.

        domain_shape = hi_root_level_refined[0:dim] - lo_root_level_refined[0:dim] + numpy.ones(dim, dtype = numpy.int)

        # Check whether the requested sub-domain is inside the computational domain.

        if numpy.all(numpy.greater_equal(lo_subdomain[0:dim], lo_root_level_refined)) == False:
            raise RuntimeError('Input sub-domain not inside the computational domain!')
        if numpy.all(numpy.less_equal(hi_subdomain[0:dim], hi_root_level_refined)) == False:
            raise RuntimeError('Input sub-domain not inside the computational domain!')

        # Compute the shape of the domain refined to the highest level.

        domain_shape = hi_root_level_refined[0:dim] - lo_root_level_refined[0:dim] \
            + numpy.ones(dim, dtype = lo_root_level_refined.dtype)

        # Compute the domain shape at each level.

        domain_shape_level = numpy.empty((num_levels, dim), dtype = domain_shape.dtype)

        for level_num in range(num_levels - 1):
            domain_shape_level[level_num] = numpy.divide(domain_shape, ratios_to_finest_level[level_num][0:dim])

        domain_shape_level[-1] = domain_shape

        # Get the number of ghost cells required for upsampling.

        num_ghosts_upsampling = self._upsampler.getNumberOfGhostCells()*numpy.ones(dim, dtype = num_ghosts.dtype)

        # Compute the lower and upper indices of the sub-domain coarsen to any level.
        # (including ghost cells requested by user and those for upsampling)

        lo_subdomain_level = numpy.empty((num_levels, dim), dtype = lo_subdomain.dtype)
        hi_subdomain_level = numpy.empty((num_levels, dim), dtype = hi_subdomain.dtype)

        for level_num in range(num_levels - 1):
            num_ghosts_level = (num_ghosts + (ratios_to_finest_level[level_num][0:dim] \
                - numpy.ones(dim, dtype = num_ghosts.dtype))) / ratios_to_finest_level[level_num][0:dim] \
                + num_ghosts_upsampling

            lo_subdomain_level[level_num] = lo_subdomain / ratios_to_finest_level[level_num][0:dim] \
                - num_ghosts_level
            hi_subdomain_level[level_num] = hi_subdomain / ratios_to_finest_level[level_num][0:dim] \
                + num_ghosts_level

        lo_subdomain_level[-1] = lo_subdomain - num_ghosts[0:dim]
        hi_subdomain_level[-1] = hi_subdomain + num_ghosts[0:dim]

        # Include the ghost cells in the sub-domain.

        lo_subdomain = lo_subdomain[0:dim] - num_ghosts[0:dim]
        hi_subdomain = hi_subdomain[0:dim] + num_ghosts[0:dim]

        # Determine which file clusters to load.

        file_clusters_to_load = []

        for level_num in range(num_levels):
            patch_level_start_idx = 0
            for level_idx in range(0, level_num):
                patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]

            for patch_idx in range(num_patches[level_num]):
                global_patch_idx = patch_level_start_idx + patch_idx
                lo_patch = self._patch_extents[global_patch_idx][0]
                hi_patch = self._patch_extents[global_patch_idx][1]

                file_cluster_num = self._patch_map[global_patch_idx][1]

                # Determine whether the sub-domain overlaps with the patches.

                if dim == 1:
                    load_file_cluster = False

                    if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                       (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                        load_file_cluster = True

                    # Check whether the patches overlap with the ghost cell regions if the domain is periodic.

                    if periodic_dimensions[0] == True:
                        # Check the left boundary.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]):
                            load_file_cluster = True

                        # Check the right boundary.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]):
                            load_file_cluster = True

                    if load_file_cluster and (file_cluster_num not in file_clusters_to_load):
                        file_clusters_to_load.append(file_cluster_num)

                elif dim == 2:
                    load_file_cluster = False

                    if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                       (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                       (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                       (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                        load_file_cluster = True

                    # Check whether the patches overlap with the ghost cell regions if the domain is periodic.

                    if periodic_dimensions[0] == True:
                        # Check the left edge.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]):
                            if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                load_file_cluster = True

                        # Check the right edge.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]):
                            if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                load_file_cluster = True

                    if periodic_dimensions[1] == True:
                        # Check the bottom edge.

                        if (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                load_file_cluster = True

                        # Check the top edge.

                        if (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                load_file_cluster = True

                    if (periodic_dimensions[0] == True) and (periodic_dimensions[1] == True):
                        # Check the left-bottom corner.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                            load_file_cluster = True

                        # Check the left-top corner.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                            load_file_cluster = True

                        # Check the right-bottom corner.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                            load_file_cluster = True

                        # Check the right-top corner.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                            load_file_cluster = True

                    if load_file_cluster and (file_cluster_num not in file_clusters_to_load):
                        file_clusters_to_load.append(file_cluster_num)

                elif dim == 3:
                    load_file_cluster = False

                    if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                       (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                       (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                       (hi_patch[1] >= lo_subdomain_level[level_num][1]) and \
                       (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                       (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                        load_file_cluster = True

                    # Check whether the patches overlap with the ghost cell regions if the domain is periodic.

                    if periodic_dimensions[0] == True:
                        # Check the left face.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]):
                            if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]) and \
                               (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                               (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                load_file_cluster = True

                        # Check the right face.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]):
                            if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]) and \
                               (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                               (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                load_file_cluster = True

                    if periodic_dimensions[1] == True:
                        # Check the bottom face.

                        if (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                               (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                               (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                load_file_cluster = True

                        # Check the top face.

                        if (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                               (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                               (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                load_file_cluster = True

                    if periodic_dimensions[2] == True:
                        # Check the back face.

                        if (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                               (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                load_file_cluster = True

                        # Check the front face.

                        if (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                               (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                load_file_cluster = True

                    if (periodic_dimensions[0] == True) and (periodic_dimensions[1] == True):
                        # Check the left-bottom edge.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                            if (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                               (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                load_file_cluster = True

                        # Check the left-top edge.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                            if (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                               (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                load_file_cluster = True

                        # Check the right-bottom edge.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                            if (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                               (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                load_file_cluster = True

                        # Check the right-top edge.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                            if (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                               (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                load_file_cluster = True

                    if (periodic_dimensions[0] == True) and (periodic_dimensions[2] == True):
                        # Check the left-back edge.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                load_file_cluster = True

                        # Check the left-front edge.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                load_file_cluster = True

                        # Check the right-back edge.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                load_file_cluster = True

                        # Check the right-front edge.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                               (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                load_file_cluster = True

                    if (periodic_dimensions[1] == True) and (periodic_dimensions[2] == True):
                        # Check the bottom-back edge.

                        if (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) and \
                           (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                load_file_cluster = True

                        # Check the bottom-front edge.

                        if (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) and \
                           (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                load_file_cluster = True

                        # Check the top-back edge.

                        if (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) and \
                           (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                load_file_cluster = True

                        # Check the top-front edge.

                        if (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) and \
                           (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                               (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                load_file_cluster = True

                    if (periodic_dimensions[0] == True) and \
                       (periodic_dimensions[1] == True) and \
                       (periodic_dimensions[2] == True):
                        # Check the left-bottom-back corner.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) and \
                           (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            load_file_cluster = True

                        # Check the left-top-back corner.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) and \
                           (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            load_file_cluster = True

                        # Check the right-bottom-back corner.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) and \
                           (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            load_file_cluster = True

                        # Check the right-top-back corner.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) and \
                           (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                            load_file_cluster = True

                        # Check the left-bottom-front corner.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) and \
                           (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            load_file_cluster = True

                        # Check the left-top-front corner.

                        if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) and \
                           (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            load_file_cluster = True

                        # Check the right-bottom-front corner.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) and \
                           (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            load_file_cluster = True

                        # Check the right-top-front corner.

                        if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) and \
                           (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) and \
                           (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                            load_file_cluster = True

                    if load_file_cluster and (file_cluster_num not in file_clusters_to_load):
                        file_clusters_to_load.append(file_cluster_num)

                else:
                    raise RuntimeError('Problem dimension < 1 or > 3 not supported!')

        # Initialize containers to store the data at different levels. The elements in the containers
        # are initialized as NAN values.

        level_data = {}

        for var_name in var_names:
            level_data[var_name] = []

            for level_num in range(num_levels):
                data_shape = hi_subdomain_level[level_num][0:dim] - lo_subdomain_level[level_num][0:dim] \
                    + numpy.ones(dim, dtype = lo_subdomain_level.dtype)

                if self._data_order == 'C':
                    data_shape = numpy.insert(data_shape, 0, var_num_components[var_name])
                else:
                    data_shape = numpy.append(data_shape, var_num_components[var_name])

                data = numpy.empty(data_shape, dtype = numpy.float64, order = self._data_order)
                data[:] = numpy.NAN

                level_data[var_name].append(data)

        if dim == 1:
            for process_idx in file_clusters_to_load:
                file_name = 'processor_cluster.' + str(process_idx).zfill(self._processor_zero_padding_length) + '.samrai'
                full_path = self._full_viz_folder_paths[self._step] + '/' + file_name
                f_input = h5py.File(full_path, 'r')

                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]

                    file_cluster = f_input['processor.' + str(process_idx).zfill(self._processor_zero_padding_length)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]

                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]

                            # Get the lower and upper indices of the current patch.

                            lo_patch = self._patch_extents[global_patch_idx][0]
                            hi_patch = self._patch_extents[global_patch_idx][1]

                            for component_idx in range(0, var_num_components[var_name]):
                                # Get the patch data.
                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]][()]

                                if self._data_order == 'C':
                                    self._loadDataFromPatchToSubdomain( \
                                        lo_subdomain_level[level_num], hi_subdomain_level[level_num], \
                                        lo_patch, hi_patch, \
                                        level_data[var_name][level_num][component_idx, :], patch_data)
                                else:
                                    self._loadDataFromPatchToSubdomain( \
                                        lo_subdomain_level[level_num], hi_subdomain_level[level_num], \
                                        lo_patch, hi_patch, \
                                        level_data[var_name][level_num][:, component_idx], patch_data)

                                # Check whether the patches overlap with the ghost cell regions if the domain is periodic.

                                lo_subdomain_shifted = numpy.empty(1, dtype = lo_subdomain_level.dtype)
                                hi_subdomain_shifted = numpy.empty(1, dtype = hi_subdomain_level.dtype)

                                if periodic_dimensions[0] == True:
                                    # Check the left boundary.
                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, component_idx], patch_data)

                                    # Check the right boundary.
                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, component_idx], patch_data)

                f_input.close()

        elif dim == 2:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(self._processor_zero_padding_length) + '.samrai'
                full_path = self._full_viz_folder_paths[self._step] + '/' + file_name
                f_input = h5py.File(full_path, 'r')

                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]

                    file_cluster = f_input['processor.' + str(process_idx).zfill(self._processor_zero_padding_length)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]

                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]

                            # Get the lower and upper indices of the current patch.

                            lo_patch = self._patch_extents[global_patch_idx][0]
                            hi_patch = self._patch_extents[global_patch_idx][1]

                            # Get the shape of the patch.

                            patch_shape = hi_patch[0:2] - lo_patch[0:2] + numpy.ones(2, dtype = lo_patch.dtype)

                            for component_idx in range(0, var_num_components[var_name]):
                                # Get the patch data.

                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]][()].reshape( \
                                    patch_shape, order = 'F')

                                if self._data_order == 'C':
                                    patch_data  = patch_data.copy(order = 'C')

                                    self._loadDataFromPatchToSubdomain( \
                                        lo_subdomain_level[level_num], hi_subdomain_level[level_num], \
                                        lo_patch, hi_patch, \
                                        level_data[var_name][level_num][component_idx, :, :], patch_data)
                                else:
                                    self._loadDataFromPatchToSubdomain( \
                                        lo_subdomain_level[level_num], hi_subdomain_level[level_num], \
                                        lo_patch, hi_patch, \
                                        level_data[var_name][level_num][:, :, component_idx], patch_data)

                                # Check whether the patches overlap with the ghost cell regions if the domain is periodic.

                                lo_subdomain_shifted = numpy.empty(2, dtype = lo_subdomain_level.dtype)
                                hi_subdomain_shifted = numpy.empty(2, dtype = hi_subdomain_level.dtype)

                                if periodic_dimensions[0] == True:
                                    # Check the left edge.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]):
                                        if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, component_idx], patch_data)

                                    # Check the right edge.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]):
                                        if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, component_idx], patch_data)

                                if periodic_dimensions[1] == True:
                                    # Check the bottom edge.

                                    if (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, component_idx], patch_data)

                                    # Check the top edge.

                                    if (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, component_idx], patch_data)

                                if (periodic_dimensions[0] == True) and (periodic_dimensions[1] == True):
                                    # Check the left-bottom corner.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, component_idx], patch_data)

                                    # Check the left-top corner.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, component_idx], patch_data)

                                    # Check the right-bottom corner.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, component_idx], patch_data)

                                    # Check the right-top corner.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                                        load_file_cluster = True
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, component_idx], patch_data)

                f_input.close()

        elif dim == 3:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(self._processor_zero_padding_length) + '.samrai'
                full_path = self._full_viz_folder_paths[self._step] + '/' + file_name
                f_input = h5py.File(full_path, 'r')

                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]

                    file_cluster = f_input['processor.' + str(process_idx).zfill(self._processor_zero_padding_length)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]

                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]

                            # Get the lower and upper indices of the current patch.

                            lo_patch = self._patch_extents[global_patch_idx][0]
                            hi_patch = self._patch_extents[global_patch_idx][1]

                            # Get the shape of the patch.

                            patch_shape = hi_patch - lo_patch + numpy.ones(3, dtype = lo_patch.dtype)

                            for component_idx in range(0, var_num_components[var_name]):
                                # Get the patch data and upsample the data to the finest resolution.

                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]][()].reshape( \
                                    patch_shape, order = 'F')

                                if self._data_order == 'C':
                                    patch_data  = patch_data.copy(order = 'C')

                                    self._loadDataFromPatchToSubdomain( \
                                        lo_subdomain_level[level_num], hi_subdomain_level[level_num], \
                                        lo_patch, hi_patch, \
                                        level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                else:
                                    self._loadDataFromPatchToSubdomain( \
                                        lo_subdomain_level[level_num], hi_subdomain_level[level_num], \
                                        lo_patch, hi_patch, \
                                        level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                # Check whether the patches overlap with the ghost cell regions if the domain is periodic.

                                lo_subdomain_shifted = numpy.empty(3, dtype = lo_subdomain_level.dtype)
                                hi_subdomain_shifted = numpy.empty(3, dtype = hi_subdomain_level.dtype)

                                if periodic_dimensions[0] == True:
                                    # Check the left face.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]):
                                        if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]) and \
                                           (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                                           (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right face.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]):
                                        if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]) and \
                                           (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                                           (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                if periodic_dimensions[1] == True:
                                    # Check the bottom face.

                                    if (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                                           (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                                           (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the top face.

                                    if (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                                           (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                                           (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                if periodic_dimensions[2] == True:
                                    # Check the back face.

                                    if (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                                           (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the front face.

                                    if (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]) and \
                                           (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                if (periodic_dimensions[0] == True) and (periodic_dimensions[1] == True):
                                    # Check the left-bottom edge.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                                        if (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                                           (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the left-top edge.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                                        if (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                                           (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right-bottom edge.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]):
                                        if (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                                           (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right-top edge.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]):
                                        if (lo_patch[2] <= hi_subdomain_level[level_num][2]) and \
                                           (hi_patch[2] >= lo_subdomain_level[level_num][2]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                if (periodic_dimensions[0] == True) and (periodic_dimensions[2] == True):
                                    # Check the left-back edge.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the left-front edge.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                + domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right-back edge.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right-front edge.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        if (lo_patch[1] <= hi_subdomain_level[level_num][1]) and \
                                           (hi_patch[1] >= lo_subdomain_level[level_num][1]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                                - domain_shape_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                if (periodic_dimensions[1] == True) and (periodic_dimensions[2] == True):
                                    # Check the bottom-back edge.

                                    if (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) \
                                       and \
                                       (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the bottom-front edge.

                                    if (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) \
                                       and \
                                       (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                + domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the top-back edge.

                                    if (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) \
                                       and \
                                       (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                + domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the top-front edge.

                                    if (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) \
                                       and \
                                       (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        if (lo_patch[0] <= hi_subdomain_level[level_num][0]) and \
                                           (hi_patch[0] >= lo_subdomain_level[level_num][0]):
                                            lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0]
                                            hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0]

                                            lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]
                                            hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                                - domain_shape_level[level_num][1]

                                            lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]
                                            hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                                - domain_shape_level[level_num][2]

                                            if self._data_order == 'C':
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                            else:
                                                self._loadDataFromPatchToSubdomain( \
                                                    lo_subdomain_shifted, hi_subdomain_shifted, \
                                                    lo_patch, hi_patch, \
                                                    level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                if (periodic_dimensions[0] == True) and \
                                   (periodic_dimensions[1] == True) and \
                                   (periodic_dimensions[2] == True):
                                    # Check the left-bottom-back corner.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) \
                                       and \
                                       (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]

                                        lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                            + domain_shape_level[level_num][2]
                                        hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                            + domain_shape_level[level_num][2]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the left-top-back corner.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) \
                                       and \
                                       (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]

                                        lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                            + domain_shape_level[level_num][2]
                                        hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                            + domain_shape_level[level_num][2]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right-bottom-back corner.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) \
                                       and \
                                       (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]

                                        lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                            + domain_shape_level[level_num][2]
                                        hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                            + domain_shape_level[level_num][2]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right-top-back corner.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) \
                                       and \
                                       (lo_subdomain_level[level_num][2] <= hi_patch[2] - domain_shape_level[level_num][2]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]

                                        lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                            + domain_shape_level[level_num][2]
                                        hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                            + domain_shape_level[level_num][2]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the left-bottom-front corner.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) \
                                       and \
                                       (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]

                                        lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                            - domain_shape_level[level_num][2]
                                        hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                            - domain_shape_level[level_num][2]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the left-top-front corner.

                                    if (lo_subdomain_level[level_num][0] <= hi_patch[0] - domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) \
                                       and \
                                       (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            + domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]

                                        lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                            - domain_shape_level[level_num][2]
                                        hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                            - domain_shape_level[level_num][2]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right-bottom-front corner.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (lo_subdomain_level[level_num][1] <= hi_patch[1] - domain_shape_level[level_num][1]) \
                                       and \
                                       (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            + domain_shape_level[level_num][1]

                                        lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                            - domain_shape_level[level_num][2]
                                        hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                            - domain_shape_level[level_num][2]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                                    # Check the right-top-front corner.

                                    if (hi_subdomain_level[level_num][0] >= lo_patch[0] + domain_shape_level[level_num][0]) \
                                       and \
                                       (hi_subdomain_level[level_num][1] >= lo_patch[1] + domain_shape_level[level_num][1]) \
                                       and \
                                       (hi_subdomain_level[level_num][2] >= lo_patch[2] + domain_shape_level[level_num][2]):
                                        lo_subdomain_shifted[0] = lo_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]
                                        hi_subdomain_shifted[0] = hi_subdomain_level[level_num][0] \
                                            - domain_shape_level[level_num][0]

                                        lo_subdomain_shifted[1] = lo_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]
                                        hi_subdomain_shifted[1] = hi_subdomain_level[level_num][1] \
                                            - domain_shape_level[level_num][1]

                                        lo_subdomain_shifted[2] = lo_subdomain_level[level_num][2] \
                                            - domain_shape_level[level_num][2]
                                        hi_subdomain_shifted[2] = hi_subdomain_level[level_num][2] \
                                            - domain_shape_level[level_num][2]

                                        if self._data_order == 'C':
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                        else:
                                            self._loadDataFromPatchToSubdomain(lo_subdomain_shifted, hi_subdomain_shifted, \
                                                lo_patch, hi_patch, \
                                                level_data[var_name][level_num][:, :, :, component_idx], patch_data)

                f_input.close()

        else:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')

        # Combine data at all levels.

        for var_name in var_names:

            data_shape = hi_subdomain_level[-1][0:dim] - lo_subdomain_level[-1][0:dim] \
                + numpy.ones(dim, dtype = lo_subdomain_level.dtype)

            if self._data_order == 'C':
                data_shape = numpy.insert(data_shape, 0, var_num_components[var_name])
            else:
                data_shape = numpy.append(data_shape, var_num_components[var_name])

            self._data[var_name] = numpy.empty(data_shape, dtype = numpy.float64, order = self._data_order)
            self._data[var_name][:] = numpy.NAN

            if self._data_order == 'C':
                data_shape = numpy.array(data_shape[1:])
            else:
                data_shape = numpy.array(data_shape[:-1])

            if dim == 1:
                lo_root_refined = numpy.empty(dim, dtype = lo_subdomain_level.dtype)

                lo_root_refined[0] = lo_subdomain_level[0][0]*ratios_to_finest_level[0][0]

                x_start_idx = lo_subdomain[0] - lo_root_refined[0]
                x_end_idx = data_shape[0] + x_start_idx

                for component_idx in range(0, var_num_components[var_name]):
                    root_data_component = self._upsampler_constant.upsample(level_data[var_name][0], \
                        ratios_to_finest_level[0], \
                        component_idx)

                    if self._data_order == 'C':
                        root_data_component = root_data_component[x_start_idx:x_end_idx]
                        is_finite_idx = numpy.isfinite(root_data_component)
                        self._data[var_name][component_idx, :][is_finite_idx] = root_data_component[is_finite_idx]

                    else:
                        root_data_component = root_data_component[x_start_idx:x_end_idx]
                        is_finite_idx = numpy.isfinite(root_data_component)
                        self._data[var_name][:, component_idx][is_finite_idx] = root_data_component[is_finite_idx]

                for level_num in range(0, num_levels):
                    lo_level_refined = numpy.empty(dim, dtype = lo_subdomain_level.dtype)

                    lo_level_refined[0] = lo_subdomain_level[level_num][0]*ratios_to_finest_level[level_num][0]

                    x_start_idx = lo_subdomain[0] - lo_level_refined[0]
                    x_end_idx = data_shape[0] + x_start_idx

                    for component_idx in range(0, var_num_components[var_name]):
                        if level_num != num_levels - 1:
                            level_data_component = self._upsampler.upsample(level_data[var_name][level_num], \
                                ratios_to_finest_level[level_num], \
                                component_idx)
                        else:
                            if self._data_order == 'C':
                                level_data_component = level_data[var_name][level_num][component_idx, :]
                            else:
                                level_data_component = level_data[var_name][level_num][:, component_idx]

                        if self._data_order == 'C':
                            level_data_component = level_data_component[x_start_idx:x_end_idx]
                            is_finite_idx = numpy.isfinite(level_data_component)
                            self._data[var_name][component_idx, :][is_finite_idx] = level_data_component[is_finite_idx]

                        else:
                            level_data_component = level_data_component[x_start_idx:x_end_idx]
                            is_finite_idx = numpy.isfinite(level_data_component)
                            self._data[var_name][:, component_idx][is_finite_idx] = level_data_component[is_finite_idx]

            elif dim == 2:
                lo_root_refined = numpy.empty(dim, dtype = lo_subdomain_level.dtype)

                lo_root_refined[0] = lo_subdomain_level[0][0]*ratios_to_finest_level[0][0]
                lo_root_refined[1] = lo_subdomain_level[0][1]*ratios_to_finest_level[0][1]

                x_start_idx = lo_subdomain[0] - lo_root_refined[0]
                x_end_idx = data_shape[0] + x_start_idx

                y_start_idx = lo_subdomain[1] - lo_root_refined[1]
                y_end_idx = data_shape[1] + y_start_idx

                for component_idx in range(0, var_num_components[var_name]):
                    root_data_component = self._upsampler_constant.upsample(level_data[var_name][0], \
                        ratios_to_finest_level[0], \
                        component_idx)

                    if self._data_order == 'C':
                        root_data_component = root_data_component[x_start_idx:x_end_idx, y_start_idx:y_end_idx]
                        is_finite_idx = numpy.isfinite(root_data_component)
                        self._data[var_name][component_idx, :, :][is_finite_idx] = root_data_component[is_finite_idx]

                    else:
                        root_data_component = root_data_component[x_start_idx:x_end_idx, y_start_idx:y_end_idx]
                        is_finite_idx = numpy.isfinite(root_data_component)
                        self._data[var_name][:, :, component_idx][is_finite_idx] = root_data_component[is_finite_idx]

                for level_num in range(0, num_levels):
                    lo_level_refined = numpy.empty(dim, dtype = lo_subdomain_level.dtype)

                    lo_level_refined[0] = lo_subdomain_level[level_num][0]*ratios_to_finest_level[level_num][0]
                    lo_level_refined[1] = lo_subdomain_level[level_num][1]*ratios_to_finest_level[level_num][1]

                    x_start_idx = lo_subdomain[0] - lo_level_refined[0]
                    x_end_idx = data_shape[0] + x_start_idx

                    y_start_idx = lo_subdomain[1] - lo_level_refined[1]
                    y_end_idx = data_shape[1] + y_start_idx

                    for component_idx in range(0, var_num_components[var_name]):
                        if level_num != num_levels - 1:
                            level_data_component = self._upsampler.upsample(level_data[var_name][level_num], \
                                ratios_to_finest_level[level_num], \
                                component_idx)
                        else:
                            if self._data_order == 'C':
                                level_data_component = level_data[var_name][level_num][component_idx, :, :]
                            else:
                                level_data_component = level_data[var_name][level_num][:, :, component_idx]

                        if self._data_order == 'C':
                            level_data_component = level_data_component[x_start_idx:x_end_idx, y_start_idx:y_end_idx]
                            is_finite_idx = numpy.isfinite(level_data_component)
                            self._data[var_name][component_idx, :, :][is_finite_idx] = level_data_component[is_finite_idx]

                        else:
                            level_data_component = level_data_component[x_start_idx:x_end_idx, y_start_idx:y_end_idx]
                            is_finite_idx = numpy.isfinite(level_data_component)
                            self._data[var_name][:, :, component_idx][is_finite_idx] = level_data_component[is_finite_idx]

            elif dim == 3:
                lo_root_refined = numpy.empty(dim, dtype = lo_subdomain_level.dtype)

                lo_root_refined[0] = lo_subdomain_level[0][0]*ratios_to_finest_level[0][0]
                lo_root_refined[1] = lo_subdomain_level[0][1]*ratios_to_finest_level[0][1]
                lo_root_refined[2] = lo_subdomain_level[0][2]*ratios_to_finest_level[0][2]

                x_start_idx = lo_subdomain[0] - lo_root_refined[0]
                x_end_idx = data_shape[0] + x_start_idx

                y_start_idx = lo_subdomain[1] - lo_root_refined[1]
                y_end_idx = data_shape[1] + y_start_idx

                z_start_idx = lo_subdomain[2] - lo_root_refined[2]
                z_end_idx = data_shape[2] + z_start_idx

                for component_idx in range(0, var_num_components[var_name]):
                    root_data_component = self._upsampler_constant.upsample(level_data[var_name][0], \
                        ratios_to_finest_level[0], \
                        component_idx)

                    if self._data_order == 'C':
                        root_data_component = \
                            root_data_component[x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx]
                        is_finite_idx = numpy.isfinite(root_data_component)
                        self._data[var_name][component_idx, :, :, :][is_finite_idx] = root_data_component[is_finite_idx]

                    else:
                        root_data_component = \
                            root_data_component[x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx]
                        is_finite_idx = numpy.isfinite(root_data_component)
                        self._data[var_name][:, :, :, component_idx][is_finite_idx] = root_data_component[is_finite_idx]

                for level_num in range(0, num_levels):
                    lo_level_refined = numpy.empty(dim, dtype = lo_subdomain_level.dtype)

                    lo_level_refined[0] = lo_subdomain_level[level_num][0]*ratios_to_finest_level[level_num][0]
                    lo_level_refined[1] = lo_subdomain_level[level_num][1]*ratios_to_finest_level[level_num][1]
                    lo_level_refined[2] = lo_subdomain_level[level_num][2]*ratios_to_finest_level[level_num][2]

                    x_start_idx = lo_subdomain[0] - lo_level_refined[0]
                    x_end_idx = data_shape[0] + x_start_idx

                    y_start_idx = lo_subdomain[1] - lo_level_refined[1]
                    y_end_idx = data_shape[1] + y_start_idx

                    z_start_idx = lo_subdomain[2] - lo_level_refined[2]
                    z_end_idx = data_shape[2] + z_start_idx

                    for component_idx in range(0, var_num_components[var_name]):
                        if level_num != num_levels - 1:
                            level_data_component = self._upsampler.upsample(level_data[var_name][level_num], \
                                ratios_to_finest_level[level_num], \
                                component_idx)
                        else:
                            if self._data_order == 'C':
                                level_data_component = level_data[var_name][level_num][component_idx, :, :, :]
                            else:
                                level_data_component = level_data[var_name][level_num][:, :, :, component_idx]

                        if self._data_order == 'C':
                            level_data_component = \
                                level_data_component[x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx]
                            is_finite_idx = numpy.isfinite(level_data_component)
                            self._data[var_name][component_idx, :, :, :][is_finite_idx] = level_data_component[is_finite_idx]

                        else:
                            level_data_component = \
                                level_data_component[x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx]
                            is_finite_idx = numpy.isfinite(level_data_component)
                            self._data[var_name][:, :, :, component_idx][is_finite_idx] = level_data_component[is_finite_idx]

        self._data_loaded = True


    def _loadDataFromPatchToSubdomain(self, \
            lo_subdomain, \
            hi_subdomain, \
            lo_patch, \
            hi_patch, \
            subdomain_data, \
            patch_data):
        """
        Private method to load data from patch to sub-domain.
        """

        dim = lo_subdomain.shape[0]

        if dim == 1:
            # Get the global start and end indices.

            x_global_start_idx = max(lo_patch[0], lo_subdomain[0])
            x_global_end_idx = min(hi_patch[0] + 1, hi_subdomain[0] + 1)

            # Get the local start and end indices in the sub-domain.

            x_local_start_idx = x_global_start_idx - lo_subdomain[0]
            x_local_start_idx = max(0, x_local_start_idx)

            x_local_end_idx = x_global_end_idx - lo_subdomain[0]
            x_local_end_idx = min(hi_subdomain[0] - lo_subdomain[0] + 1, x_local_end_idx)

            # Get the local start and end indices in the current patch.

            x_patch_start_idx = x_global_start_idx - lo_patch[0]
            x_patch_start_idx = max(0, x_patch_start_idx)

            x_patch_end_idx = patch_data.shape[0] - ((hi_patch[0] + 1) - x_global_end_idx)
            x_patch_end_idx = min(patch_data.shape[0], x_patch_end_idx)

            if x_local_start_idx < x_local_end_idx:
                subdomain_data[x_local_start_idx:x_local_end_idx] = \
                    patch_data[x_patch_start_idx:x_patch_end_idx]

        elif dim == 2:
            # Get the global start and end indices.

            x_global_start_idx = max(lo_patch[0], lo_subdomain[0])
            x_global_end_idx = min(hi_patch[0] + 1, hi_subdomain[0] + 1)

            y_global_start_idx = max(lo_patch[1], lo_subdomain[1])
            y_global_end_idx = min(hi_patch[1] + 1, hi_subdomain[1] + 1)

            # Get the local start and end indices in the sub-domain.

            x_local_start_idx = x_global_start_idx - lo_subdomain[0]
            x_local_start_idx = max(0, x_local_start_idx)

            x_local_end_idx = x_global_end_idx - lo_subdomain[0]
            x_local_end_idx = min(hi_subdomain[0] - lo_subdomain[0] + 1, x_local_end_idx)

            y_local_start_idx = y_global_start_idx - lo_subdomain[1]
            y_local_start_idx = max(0, y_local_start_idx)

            y_local_end_idx = y_global_end_idx - lo_subdomain[1]
            y_local_end_idx = min(hi_subdomain[1] - lo_subdomain[1] + 1, y_local_end_idx)

            # Get the local start and end indices in the current patch.

            x_patch_start_idx = x_global_start_idx - lo_patch[0]
            x_patch_start_idx = max(0, x_patch_start_idx)

            x_patch_end_idx = patch_data.shape[0] - ((hi_patch[0] + 1) - x_global_end_idx)
            x_patch_end_idx = min(patch_data.shape[0], x_patch_end_idx)

            y_patch_start_idx = y_global_start_idx - lo_patch[1]
            y_patch_start_idx = max(0, y_patch_start_idx)

            y_patch_end_idx = patch_data.shape[1] - ((hi_patch[1] + 1) - y_global_end_idx)
            y_patch_end_idx = min(patch_data.shape[1], y_patch_end_idx)

            if (x_local_start_idx < x_local_end_idx) and (y_local_start_idx < y_local_end_idx):
                subdomain_data[x_local_start_idx:x_local_end_idx, y_local_start_idx:y_local_end_idx] = \
                    patch_data[x_patch_start_idx:x_patch_end_idx, y_patch_start_idx:y_patch_end_idx]

        elif dim == 3:
            # Get the global start and end indices.

            x_global_start_idx = max(lo_patch[0], lo_subdomain[0])
            x_global_end_idx = min(hi_patch[0] + 1, hi_subdomain[0] + 1)

            y_global_start_idx = max(lo_patch[1], lo_subdomain[1])
            y_global_end_idx = min(hi_patch[1] + 1, hi_subdomain[1] + 1)

            z_global_start_idx = max(lo_patch[2], lo_subdomain[2])
            z_global_end_idx = min(hi_patch[2] + 1, hi_subdomain[2] + 1)

            # Get the local start and end indices in the sub-domain.

            x_local_start_idx = x_global_start_idx - lo_subdomain[0]
            x_local_start_idx = max(0, x_local_start_idx)

            x_local_end_idx = x_global_end_idx - lo_subdomain[0]
            x_local_end_idx = min(hi_subdomain[0] - lo_subdomain[0] + 1, x_local_end_idx)

            y_local_start_idx = y_global_start_idx - lo_subdomain[1]
            y_local_start_idx = max(0, y_local_start_idx)

            y_local_end_idx = y_global_end_idx - lo_subdomain[1]
            y_local_end_idx = min(hi_subdomain[1] - lo_subdomain[1] + 1, y_local_end_idx)

            z_local_start_idx = z_global_start_idx - lo_subdomain[2]
            z_local_start_idx = max(0, z_local_start_idx)

            z_local_end_idx = z_global_end_idx - lo_subdomain[2]
            z_local_end_idx = min(hi_subdomain[2] - lo_subdomain[2] + 1, z_local_end_idx)

            # Get the local start and end indices in the current patch.

            x_patch_start_idx = x_global_start_idx - lo_patch[0]
            x_patch_start_idx = max(0, x_patch_start_idx)

            x_patch_end_idx = patch_data.shape[0] - ((hi_patch[0] + 1) - x_global_end_idx)
            x_patch_end_idx = min(patch_data.shape[0], x_patch_end_idx)

            y_patch_start_idx = y_global_start_idx - lo_patch[1]
            y_patch_start_idx = max(0, y_patch_start_idx)

            y_patch_end_idx = patch_data.shape[1] - ((hi_patch[1] + 1) - y_global_end_idx)
            y_patch_end_idx = min(patch_data.shape[1], y_patch_end_idx)

            z_patch_start_idx = z_global_start_idx - lo_patch[2]
            z_patch_start_idx = max(0, z_patch_start_idx)

            z_patch_end_idx = patch_data.shape[2] - ((hi_patch[2] + 1) - z_global_end_idx)
            z_patch_end_idx = min(patch_data.shape[2], z_patch_end_idx)

            if (x_local_start_idx < x_local_end_idx) and \
               (y_local_start_idx < y_local_end_idx) and \
               (z_local_start_idx < z_local_end_idx):
                subdomain_data[x_local_start_idx:x_local_end_idx, \
                                      y_local_start_idx:y_local_end_idx, \
                                      z_local_start_idx:z_local_end_idx] = \
                    patch_data[x_patch_start_idx:x_patch_end_idx, \
                               y_patch_start_idx:y_patch_end_idx, \
                               z_patch_start_idx:z_patch_end_idx]

        else:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')


    def readCoordinates(self):
        """
        Get the coordinates of the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        """

        # Get the dimension of the problem.

        dim = self._basic_info['dim']

        x_c = []
        y_c = []
        z_c = []

        if dim == 1:
            return self.getCombinedCoordinatesInSubdomainFromAllLevels()

        elif dim == 2:
            x_coords, y_coords = self.getCombinedCoordinatesInSubdomainFromAllLevels()

            if self._data_order == 'C':
                x_c, y_c = numpy.meshgrid(x_coords, y_coords, indexing = 'ij', sparse = False, copy=True)
            else:
                x_c, y_c = numpy.meshgrid(x_coords, y_coords, indexing = 'ij', sparse = False, copy=False)
                x_c = numpy.asfortranarray(x_c)
                y_c = numpy.asfortranarray(y_c)

            return x_c, y_c

        elif dim == 3:
            x_coords, y_coords, z_coords = self.getCombinedCoordinatesInSubdomainFromAllLevels()

            if self._data_order == 'C':
                x_c, y_c, z_c = numpy.meshgrid(x_coords, y_coords, z_coords, indexing = 'ij', sparse = False, copy=True)
            else:
                x_c, y_c, z_c = numpy.meshgrid(x_coords, y_coords, z_coords, indexing = 'ij', sparse = False, copy=False)
                x_c = numpy.asfortranarray(x_c)
                y_c = numpy.asfortranarray(y_c)
                z_c = numpy.asfortranarray(z_c)

            return x_c, y_c, z_c


    def readData(self, var_names, data=None):
        """
        Read the data of several variables in the stored sub-domain.
        Default to the full domain when the sub-domain is not set.
        """
        if self._data_loaded == True:
            self.clearData()

        # If a simple string is passed in, convert to a tuple.
        if isinstance(var_names, str):
            var_names = (var_names,)

        self.readCombinedDataInSubdomainFromAllLevels(var_names)

        dim = self._basic_info['dim']

        if data == None:
            _data = []
            for i in range(len(var_names)):
                if self._data_order == 'C':
                    if self._data[var_names[i]].shape[0] == 1:
                        _data.append(numpy.squeeze(self._data[var_names[i]], 0))
                    else:
                        _data.append(self._data[var_names[i]])
                else:
                    if self._data[var_names[i]].shape[-1] == 1:
                        _data.append(numpy.squeeze(self._data[var_names[i]], dim))
                    else:
                        _data.append(self._data[var_names[i]])
            return tuple(_data)
        else:
            for i in range(len(var_names)):
                if self._data_order == 'C':
                    if self._data[var_names[i]].shape[0] == 1:
                        data[i] = numpy.squeeze(self._data[var_names[i]], 0)
                    else:
                        data[i] = self._data[var_names[i]]
                else:
                    if self._data[var_names[i]].shape[-1] == 1:
                        data[i] = numpy.squeeze(self._data[var_names[i]], dim)
                    else:
                        data[i] = self._data[var_names[i]]

BaseReader.register(SamraiDataReader)

if __name__ == '__main__':
    print('Subclass: %s' % issubclass(SamraiDataReader, BaseReader))
    print('Instance: %s' % isinstance(SamraiDataReader("../tests/test_data_samrai_AMR"), BaseReader))
