import numpy

import explicit.first as first_der
import explicit.second as second_der

class Differentiator(object):
    """
    Class to perform derivatives with explicit finite difference schemes.
    """

    def __init__(self, grid_spacing, order, dimension=3, data_order='F'):
        """
        Constructor of the class.

        grid_spacing : iterable of floats for the grid spacing in each direction
        order : integer iterable with each value in {2, 4, 6} representing the order of accuracy of the derivatives
                in each direction
        dimension : dimension of problem
        data_order : a string {'F', 'C'} describing whether the data is in Fortran or C contiguous layout
        """

        if dimension < 1 or dimension > 3:
            raise RuntimeError("Class only works with data with number of dimensions between 1 and 3!")

        self._dim = dimension

        if len(grid_spacing) < self._dim:
            raise RuntimeError("Size of 'grid_spacing' is smaller than problem dimension!")

        if len(order) < self._dim:
            raise RuntimeError("Size of 'order' is smaller than problem dimension!")

        for i in range(self._dim):
            if order[i] not in [2, 4, 6]:
                raise RuntimeError("order[%d] has to be one of {2, 4, 6}" %i)

        self._order = tuple(order)

        self._dx = grid_spacing[0]
        if self._dim > 1:
            self._dy = grid_spacing[1]
        if self._dim > 2:
            self._dz = grid_spacing[2]

        if data_order != 'C' and data_order != 'F':
            raise RuntimeError("Invalid data order! Data order can only be 'C' or 'F'.")

        self._data_order = data_order


    def getNumberOfGhostCells(self):
        """
        Determine the number of ghost cells needed for the chosen explicit finite difference scheme.
        """

        n_ghosts = numpy.empty(self._dim, dtype=numpy.int32)

        if self._order[0] == 2:
            n_ghosts[0] = 1
        elif self._order[0] == 4:
            n_ghosts[0] = 2
        elif self._order[0] == 6:
            n_ghosts[0] = 3

        if self._dim > 1:
            if self._order[1] == 2:
                n_ghosts[1] = 1
            elif self._order[1] == 4:
                n_ghosts[1] = 2
            elif self._order[1] == 6:
                n_ghosts[1] = 3

        if self._dim > 2:
            if self._order[2] == 2:
                n_ghosts[2] = 1
            elif self._order[2] == 4:
                n_ghosts[2] = 2
            elif self._order[2] == 6:
                n_ghosts[2] = 3

        return n_ghosts


    @property
    def num_ghosts(self):
        """
        Return the number of ghost cells needed for the chosen explicit finite difference scheme.
        """

        return self.getNumberOfGhostCells()


    def ddx(self, data, der=None, component_idx=None, use_one_sided=False):
        """
        Method to compute the first order derivative of data in first direction.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem dimension
        der : optional output numpy array in the chosen contiguous layout. This array must be consistent with the problem
              dimension. This method will return der if der is None
        component_idx : index of component in data for taking derivative. None if there is only one component in the
                        data
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            if len(data_shape) != 1:
                raise RuntimeError("Make sure data is 1D!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")

        return_der = True
        if der is None:
            der = numpy.empty(data_shape, dtype=numpy.float64, order='F')
        else:
            if der.shape != data_shape:
                raise RuntimeError("Make sure shape of der is consistent with that of data!")
            return_der = False

        if self._order[0] == 2:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateSecondOrderFiniteDifference(\
                data, self._dx, 0, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[0] == 4:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateFourthOrderFiniteDifference(\
                data, self._dx, 0, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[0] == 6:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateSixthOrderFiniteDifference(\
                data, self._dx, 0, component_idx, use_one_sided, self._dim, self._data_order)

        if return_der:
            return der


    def ddy(self, data, der=None, component_idx=None, use_one_sided=False):
        """
        Method to compute the first order derivative of data in second direction.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem dimension
        der : optional output numpy array in the chosen contiguous layout. This array must be consistent with the problem
              dimension. This method will return der if der is None
        component_idx : index of component in data for taking derivative. None if there is only one component in the
                        data
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            if len(data_shape) != 1:
                raise RuntimeError("Make sure data is 1D!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")

        return_der = True
        if der is None:
            der = numpy.empty(data_shape, dtype=numpy.float64, order='F')
        else:
            if der.shape != data_shape:
                raise RuntimeError("Make sure shape of der is consistent with that of data!")
            return_der = False

        if self._order[1] == 2:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateSecondOrderFiniteDifference(\
                data, self._dy, 1, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[1] == 4:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateFourthOrderFiniteDifference(\
                data, self._dy, 1, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[1] == 6:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateSixthOrderFiniteDifference(\
                data, self._dy, 1, component_idx, use_one_sided, self._dim, self._data_order)

        if return_der:
            return der


    def ddz(self, data, der=None, component_idx=None, use_one_sided=False):
        """
        Method to compute the first order derivative of data in third direction.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem dimension
        der : optional output numpy array in the chosen contiguous layout. This array must be consistent with the problem
              dimension. This method will return der if der is None
        component_idx : index of component in data for taking derivative. None if there is only one component in the
                        data
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            if len(data_shape) != 1:
                raise RuntimeError("Make sure data is 1D!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")

        return_der = True
        if der is None:
            der = numpy.empty(data_shape, dtype=numpy.float64, order='F')
        else:
            if der.shape != data_shape:
                raise RuntimeError("Make sure shape of der is consistent with that of data!")
            return_der = False

        if self._order[2] == 2:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateSecondOrderFiniteDifference(\
                data, self._dz, 2, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[2] == 4:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateFourthOrderFiniteDifference(\
                data, self._dz, 2, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[2] == 6:
            der[tuple([slice(None)]*self._dim)] = first_der.differentiateSixthOrderFiniteDifference(\
                data, self._dz, 2, component_idx, use_one_sided, self._dim, self._data_order)

        if return_der:
            return der


    def d2dx2(self, data, der=None, component_idx=None, use_one_sided=False):
        """
        Method to compute the second order derivative of data in first direction.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem dimension
        der : optional output numpy array in the chosen contiguous layout. This array must be consistent with the problem
              dimension. This method will return der if der is None
        component_idx : index of component in data for taking derivative. None if there is only one component in the
                        data
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            if len(data_shape) != 1:
                raise RuntimeError("Make sure data is 1D!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")

        return_der = True
        if der is None:
            der = numpy.empty(data_shape, dtype=numpy.float64, order='F')
        else:
            if der.shape != data_shape:
                raise RuntimeError("Make sure shape of der is consistent with that of data!")
            return_der = False

        if self._order[0] == 2:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateSecondOrderFiniteDifference(\
                data, self._dx, 0, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[0] == 4:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateFourthOrderFiniteDifference(\
                data, self._dx, 0, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[0] == 6:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateSixthOrderFiniteDifference(\
                data, self._dx, 0, component_idx, use_one_sided, self._dim, self._data_order)

        if return_der:
            return der


    def d2dy2(self, data, der=None, component_idx=None, use_one_sided=False):
        """
        Method to compute the second order derivative of data in second direction.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem dimension
        der : optional output numpy array in the chosen contiguous layout. This array must be consistent with the problem
              dimension. This method will return der if der is None
        component_idx : index of component in data for taking derivative. None if there is only one component in the
                        data
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            if len(data_shape) != 1:
                raise RuntimeError("Make sure data is 1D!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")

        return_der = True
        if der is None:
            der = numpy.empty(data_shape, dtype=numpy.float64, order='F')
        else:
            if der.shape != data_shape:
                raise RuntimeError("Make sure shape of der is consistent with that of data!")
            return_der = False

        if self._order[1] == 2:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateSecondOrderFiniteDifference(\
                data, self._dy, 1, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[1] == 4:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateFourthOrderFiniteDifference(\
                data, self._dy, 1, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[1] == 6:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateSixthOrderFiniteDifference(\
                data, self._dy, 1, component_idx, use_one_sided, self._dim, self._data_order)

        if return_der:
            return der


    def d2dz2(self, data, der=None, component_idx=None, use_one_sided=False):
        """
        Method to compute the second order derivative of data in third direction.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem dimension
        der : optional output numpy array in the chosen contiguous layout. This array must be consistent with the problem
              dimension. This method will return der if der is None
        component_idx : index of component in data for taking derivative. None if there is only one component in the
                        data
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        """

        data_shape = data.shape

        if component_idx is not None:
            data_shape = data_shape[0:-1]

        if self._dim == 1:
            if len(data_shape) != 1:
                raise RuntimeError("Make sure data is 1D!")

        elif self._dim == 2:
            if len(data_shape) != 2:
                raise RuntimeError("Make sure data is 2D!")

        else:
            if len(data_shape) != 3:
                raise RuntimeError("Make sure data is 3D!")

        return_der = True
        if der is None:
            der = numpy.empty(data_shape, dtype=numpy.float64, order='F')
        else:
            if der.shape != data_shape:
                raise RuntimeError("Make sure shape of der is consistent with that of data!")
            return_der = False

        if self._order[2] == 2:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateSecondOrderFiniteDifference(\
                data, self._dz, 2, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[2] == 4:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateFourthOrderFiniteDifference(\
                data, self._dz, 2, component_idx, use_one_sided, self._dim, self._data_order)
        elif self._order[2] == 6:
            der[tuple([slice(None)]*self._dim)] = second_der.differentiateSixthOrderFiniteDifference(\
                data, self._dz, 2, component_idx, use_one_sided, self._dim, self._data_order)

        if return_der:
            return der


    def gradient(self, data, component_idx=None, use_one_sided=False):
        """
        Method to compute the gradient of data.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem
               dimension
        component_idx : index of component in data for taking derivative. None if there is only one component in the
                        data
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        gradient_* : returned output numpy array in the chosen contiguous layout. This array must be consistent with
                     the problem dimension
        """

        if self._dim == 1:
            return self.ddx(data, component_idx=component_idx, use_one_sided=use_one_sided)

        elif self._dim == 2:
            return self.ddx(data, component_idx=component_idx,use_one_sided=use_one_sided), \
                   self.ddy(data, component_idx=component_idx, use_one_sided=use_one_sided)

        else:
            return self.ddx(data, component_idx=component_idx, use_one_sided=use_one_sided), \
                   self.ddy(data, component_idx=component_idx, use_one_sided=use_one_sided), \
                   self.ddz(data, component_idx=component_idx, use_one_sided=use_one_sided)


    def divergence(self, data, use_one_sided=False):
        """
        Method to compute the gradient of a vector.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem
               dimension. The number of components should be as same as the number of dimensions
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        divergence : output numpy array in the chosen contiguous layout. This array is consistent with the problem
                     dimension
        """

        data_shape = data.shape

        if self._dim == 1 and len(data_shape) != 2:
            raise RuntimeError("Make sure data is 1D and has enough number of components!")

        if self._dim == 2 and len(data_shape) != 3:
            raise RuntimeError("Make sure data is 2D and has enough number of components!")

        if self._dim == 3 and len(data_shape) != 4:
            raise RuntimeError("Make sure data is 3D and has enough number of components!")

        divergence = self.ddx(data, component_idx=0, use_one_sided=use_one_sided)
        if self._dim >= 2:
            divergence = divergence + self.ddy(data, component_idx=1, use_one_sided=use_one_sided)
        if self._dim == 3:
            divergence = divergence + self.ddz(data, component_idx=2, use_one_sided=use_one_sided)

        return divergence


    def curl(self, data, use_one_sided=False):
        """
        Method to compute the curl of a vector.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem
               dimension. The number of components should be as same as the number of dimensions
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        curl : output numpy array in the chosen contiguous layout. This array is consistent with the problem
               dimension
        """

        data_shape = data.shape

        if self._dim == 1:
            raise RuntimeError("There is no curl for 1D problem!")

        if self._dim == 2 and len(data_shape) != 3:
            raise RuntimeError("Make sure data is 2D and has enough number of components!")

        if self._dim == 3 and len(data_shape) != 4:
            raise RuntimeError("Make sure data is 3D and has enough number of components!")

        curl = None
        if self._dim == 2:
            dvdx = self.ddx(data, component_idx=1, use_one_sided=use_one_sided)
            dudy = self.ddy(data, component_idx=0, use_one_sided=use_one_sided)

            curl = dvdx - dudy

        if self._dim == 3:
            dvdx = self.ddx(data, component_idx=1, use_one_sided=use_one_sided)
            dwdx = self.ddx(data, component_idx=2, use_one_sided=use_one_sided)

            dudy = self.ddy(data, component_idx=0, use_one_sided=use_one_sided)
            dwdy = self.ddy(data, component_idx=2, use_one_sided=use_one_sided)

            dudz = self.ddz(data, component_idx=0, use_one_sided=use_one_sided)
            dvdz = self.ddz(data, component_idx=1, use_one_sided=use_one_sided)

            if self._data_order == 'C':
                curl = numpy.empty( (3, data_shape[1], data_shape[2], data_shape[3]), dtype=numpy.float64, order='C' )

                curl[0, :, :, :] = dwdy - dvdz
                curl[1, :, :, :] = dudz - dwdx
                curl[2, :, :, :] = dvdx - dudy

            else:
                curl = numpy.empty( (data_shape[0], data_shape[1], data_shape[2], 3), dtype=numpy.float64, order='F' )

                curl[:, :, :, 0] = dwdy - dvdz
                curl[:, :, :, 1] = dudz - dwdx
                curl[:, :, :, 2] = dvdx - dudy

        return curl


    def laplacian(self, data, component_idx=None, use_one_sided=False):
        """
        Method to compute the laplacian of data.

        data : input numpy array in the chosen contiguous layout. This array must be consistent with the problem
               dimension
        component_idx : index of component in data for taking derivative. None if there is only one component in the
                        data
        use_one_sided : boolean to decide whether to use one-sided scheme at the boundaries
        laplacian : returned output numpy array in the chosen contiguous layout. This array must be consistent with
                      the problem dimension
        """

        laplacian = self.d2dx2(data, component_idx=component_idx, use_one_sided=use_one_sided)
        if self._dim >= 2:
            laplacian = laplacian + self.d2dy2(data, component_idx=component_idx, use_one_sided=use_one_sided)
        if self._dim == 3:
            laplacian = laplacian + self.d2dz2(data, component_idx=component_idx, use_one_sided=use_one_sided)

        return laplacian