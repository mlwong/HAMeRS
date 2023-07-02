"""
Functions for computing second order deriatives with explicit finite differencing.
"""

import numpy

def differentiateSecondOrderFiniteDifference(data, dx, direction, component_idx, use_one_sided, dimension, data_order):
    """
    Compute second order derivative using explicit second order finite differencing.
    """

    # Get the shape of data.

    data_shape = numpy.array(data.shape)

    # Check whether the direction is valid.

    if direction < 0 or direction > 2:
        raise RuntimeError('Direction < 0 or > 2 is invalid!')

    # Check whether the shape of data is valid.

    if component_idx is None:
        if data.ndim < 1 or data.ndim > 3:
            raise RuntimeError('Shape of data is invalid!')

    else:
        if data.ndim < 2 or data.ndim > 4:
            raise RuntimeError('Shape of data is invalid!')

        # Check whether the component_idx is valid and get the shape of the component's data.

        if data_order == 'C':
            if component_idx >= data.shape[0] or component_idx < 0:
                raise RuntimeError('Component index is invalid!')

            data_shape = numpy.array(data_shape[1:])

        else:
            if component_idx >= data.shape[-1] or component_idx < 0:
                raise RuntimeError('Component index is invalid!')

            data_shape = numpy.array(data_shape[:-1])

    # Get the dimension of data.

    dim = dimension

    # Check whether data size is large enough for second order second derivative.

    if use_one_sided == True:
        if direction == 0:
            if data_shape[0] < 4:
                raise RuntimeError('First dimension of data is not large enough!')

        elif direction == 1:
            if data_shape[1] < 4:
                raise RuntimeError('Second dimension of data is not large enough!')

        elif direction == 2:
            if data_shape[2] < 4:
                raise RuntimeError('Third dimension of data is not large enough!')

    else:
        if direction == 0:
            if data_shape[0] < 3:
                raise RuntimeError('First dimension of data is not large enough!')

        elif direction == 1:
            if data_shape[1] < 3:
                raise RuntimeError('Second dimension of data is not large enough!')

        elif direction == 2:
            if data_shape[2] < 3:
                raise RuntimeError('Third dimension of data is not large enough!')

    # Initialize container to store the derivatives. The elements in the container
    # are initialized as NAN values.

    diff_data = numpy.empty(data_shape, dtype=data.dtype, order=data_order)
    diff_data[:] = numpy.NAN

    # Get the component's data.

    data_component = None

    if component_idx is None:
        data_component = data
    else:
        if data_order == 'C':
            if dim == 1:
                data_component = data[component_idx, :]
            elif dim == 2:
                data_component = data[component_idx, :, :]
            elif dim == 3:
                data_component = data[component_idx, :, :, :]
        else:
            if dim == 1:
                data_component = data[:, component_idx]
            elif dim == 2:
                data_component = data[:, :, component_idx]
            elif dim == 3:
                data_component = data[:, :, :, component_idx]

    # Compute the derivatives in the interior of the domain.

    if direction == 0:
        if dim == 1:
            diff_data[1:-1] = (data_component[0:-2] - 2.0*data_component[1:-1] \
                               + data_component[2:])/(dx*dx)

        elif dim == 2:
            diff_data[1:-1, :] = (data_component[0:-2, :] - 2.0*data_component[1:-1, :] \
                                  + data_component[2:, :])/(dx*dx)

        elif dim == 3:
            diff_data[1:-1, :, :] = (data_component[0:-2, :, :] - 2.0*data_component[1:-1, :, :] \
                                  + data_component[2:, :, :])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    elif direction == 1:
        if dim < 2:
            raise IOError('There is no second direction in data with less than two dimensions!')

        elif dim == 2:
            diff_data[:, 1:-1] = (data_component[:, 0:-2] - 2.0*data_component[:, 1:-1] \
                                  + data_component[:, 2:])/(dx*dx)

        elif dim == 3:
            diff_data[:, 1:-1, :] = (data_component[:, 0:-2, :] - 2.0*data_component[:, 1:-1, :] \
                                     + data_component[:, 2:, :])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    elif direction == 2:
        if dim < 3:
            raise IOError('There is no third direction in data with less than three dimensions!')

        elif dim == 3:
            diff_data[:, :, 1:-1] = (data_component[:, :, 0:-2] - 2.0*data_component[:, :, 1:-1] \
                                     + data_component[:, :, 2:])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    # Compute the derivatives at the boundaries.

    if use_one_sided == True:
        if direction == 0:
            if dim == 1:
                diff_data[0] = (2.0*data_component[0] - 5.0*data_component[1] \
                                + 4.0*data_component[2] - data_component[3])/(dx*dx)

                diff_data[-1] = (-data_component[-4] + 4.0*data_component[-3] \
                                 - 5.0*data_component[-2] + 2.0*data_component[-1])/(dx*dx)

            elif dim == 2:
                diff_data[0, :] = (2.0*data_component[0, :] - 5.0*data_component[1, :] \
                                   + 4.0*data_component[2, :] - data_component[3, :])/(dx*dx)

                diff_data[-1, :] = (-data_component[-4, :] + 4.0*data_component[-3, :] \
                                    - 5.0*data_component[-2, :] + 2.0*data_component[-1, :])/(dx*dx)


            elif dim == 3:
                diff_data[0, :, :] = (2.0*data_component[0, :, :] - 5.0*data_component[1, :, :] \
                                      + 4.0*data_component[2, :, :] - data_component[3, :, :])/(dx*dx)

                diff_data[-1, :, :] = (-data_component[-4, :, :] + 4.0*data_component[-3, :, :] \
                                       - 5.0*data_component[-2, :, :] + 2.0*data_component[-1, :, :])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

        elif direction == 1:
            if dim < 2:
                raise RuntimeError('There is no second direction in data with less than two dimensions!')

            elif dim == 2:
                diff_data[:, 0] = (2.0*data_component[:, 0] - 5.0*data_component[:, 1] \
                                   + 4.0*data_component[:, 2] - data_component[:, 3])/(dx*dx)

                diff_data[:, -1] = (-data_component[:, -4] + 4.0*data_component[:, -3] \
                                    - 5.0*data_component[:, -2] + 2.0*data_component[:, -1])/(dx*dx)

            elif dim == 3:
                diff_data[:, 0, :] = (2.0*data_component[:, 0, :] - 5.0*data_component[:, 1, :] \
                                   + 4.0*data_component[:, 2, :] - data_component[:, 3, :])/(dx*dx)

                diff_data[:, -1, :] = (-data_component[:, -4, :] + 4.0*data_component[:, -3, :] \
                                    - 5.0*data_component[:, -2, :] + 2.0*data_component[:, -1, :])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

        elif direction == 2:
            if dim < 3:
                raise IOError('There is no third direction in data with less than three dimensions!')

            elif dim == 3:
                diff_data[:, :, 0] = (2.0*data_component[:, :, 0] - 5.0*data_component[:, :, 1] \
                                      + 4.0*data_component[:, :, 2] - data_component[:, :, 3])/(dx*dx)

                diff_data[:, :, -1] = (-data_component[:, :, -4] + 4.0*data_component[:, :, -3] \
                                       - 5.0*data_component[:, :, -2] + 2.0*data_component[:, :, -1])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

    return diff_data


def differentiateFourthOrderFiniteDifference(data, dx, direction, component_idx, use_one_sided, dimension, data_order):
    """
    Compute second order derivative using explicit fourth order finite differencing.
    """

    # Get the shape of data.

    data_shape = numpy.array(data.shape)

    # Check whether the direction is valid.

    if direction < 0 or direction > 2:
        raise RuntimeError('Direction < 0 or > 2 is invalid!')

    # Check whether the shape of data is valid.

    if component_idx is None:
        if data.ndim < 1 or data.ndim > 3:
            raise RuntimeError('Shape of data is invalid!')

    else:
        if data.ndim < 2 or data.ndim > 4:
            raise RuntimeError('Shape of data is invalid!')

        # Check whether the component_idx is valid and get the shape of the component's data.

        if data_order == 'C':
            if component_idx >= data.shape[0] or component_idx < 0:
                raise RuntimeError('Component index is invalid!')

            data_shape = numpy.array(data_shape[1:])

        else:
            if component_idx >= data.shape[-1] or component_idx < 0:
                raise RuntimeError('Component index is invalid!')

            data_shape = numpy.array(data_shape[:-1])

    # Get the dimension of data.

    dim = dimension

    # Check whether data size is large enough for fourth order second derivative.

    if use_one_sided == True:
        if direction == 0:
            if data_shape[0] < 6:
                raise RuntimeError('First dimension of data is not large enough!')

        elif direction == 1:
            if data_shape[1] < 6:
                raise RuntimeError('Second dimension of data is not large enough!')

        elif direction == 2:
            if data_shape[2] < 6:
                raise RuntimeError('Third dimension of data is not large enough!')

    else:
        if direction == 0:
            if data_shape[0] < 5:
                raise RuntimeError('First dimension of data is not large enough!')

        elif direction == 1:
            if data_shape[1] < 5:
                raise RuntimeError('Second dimension of data is not large enough!')

        elif direction == 2:
            if data_shape[2] < 5:
                raise RuntimeError('Third dimension of data is not large enough!')

    # Initialize container to store the derivatives. The elements in the container
    # are initialized as NAN values.

    diff_data = numpy.empty(data_shape, dtype=data.dtype, order=data_order)
    diff_data[:] = numpy.NAN

    # Get the component's data.

    data_component = None

    if component_idx is None:
        data_component = data
    else:
        if data_order == 'C':
            if dim == 1:
                data_component = data[component_idx, :]
            elif dim == 2:
                data_component = data[component_idx, :, :]
            elif dim == 3:
                data_component = data[component_idx, :, :, :]
        else:
            if dim == 1:
                data_component = data[:, component_idx]
            elif dim == 2:
                data_component = data[:, :, component_idx]
            elif dim == 3:
                data_component = data[:, :, :, component_idx]

    # Compute the derivatives in the interior of the domain.

    if direction == 0:
        if dim == 1:
            diff_data[2:-2] = (-1.0/12.0*data_component[0:-4] + 4.0/3.0*data_component[1:-3] \
                               - 5.0/2.0*data_component[2:-2] + 4.0/3.0*data_component[3:-1] \
                               - 1.0/12.0*data_component[4:])/(dx*dx)

        elif dim == 2:
            diff_data[2:-2, :] = (-1.0/12.0*data_component[0:-4, :] + 4.0/3.0*data_component[1:-3, :] \
                                  - 5.0/2.0*data_component[2:-2, :] + 4.0/3.0*data_component[3:-1, :] \
                                  - 1.0/12.0*data_component[4:, :])/(dx*dx)

        elif dim == 3:
            diff_data[2:-2, :, :] = (-1.0/12.0*data_component[0:-4, :, :] + 4.0/3.0*data_component[1:-3, :, :] \
                                     - 5.0/2.0*data_component[2:-2, :, :] + 4.0/3.0*data_component[3:-1, :, :] \
                                     - 1.0/12.0*data_component[4:, :, :])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    elif direction == 1:
        if dim < 2:
            raise IOError('There is no second direction in data with less than two dimensions!')

        elif dim == 2:
            diff_data[:, 2:-2] = (-1.0/12.0*data_component[:, 0:-4] + 4.0/3.0*data_component[:, 1:-3] \
                                  - 5.0/2.0*data_component[:, 2:-2] + 4.0/3.0*data_component[:, 3:-1] \
                                  - 1.0/12.0*data_component[:, 4:])/(dx*dx)

        elif dim == 3:
            diff_data[:, 2:-2, :] = (-1.0/12.0*data_component[:, 0:-4, :] + 4.0/3.0*data_component[:, 1:-3, :] \
                                     - 5.0/2.0*data_component[:, 2:-2, :] + 4.0/3.0*data_component[:, 3:-1, :] \
                                     - 1.0/12.0*data_component[:, 4:, :])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    elif direction == 2:
        if dim < 3:
            raise IOError('There is no third direction in data with less than three dimensions!')

        elif dim == 3:
            diff_data[:, :, 2:-2] = (-1.0/12.0*data_component[:, :, 0:-4] + 4.0/3.0*data_component[:, :, 1:-3] \
                                     - 5.0/2.0*data_component[:, :, 2:-2] + 4.0/3.0*data_component[:, :, 3:-1] \
                                     - 1.0/12.0*data_component[:, :, 4:])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    # Compute the derivatives at the boundaries.

    if use_one_sided == True:
        if direction == 0:
            if dim == 1:
                diff_data[0] = (15.0/4.0*data_component[0] - 77.0/6.0*data_component[1] \
                                + 107.0/6.0*data_component[2] - 13.0*data_component[3] \
                                + 61.0/12.0*data_component[4] - 5.0/6.0*data_component[5])/(dx*dx)

                diff_data[1] = (5.0/6.0*data_component[0] - 5.0/4.0*data_component[1] \
                                - 1.0/3.0*data_component[2] + 7.0/6.0*data_component[3] \
                                - 1.0/2.0*data_component[4] + 1.0/12.0*data_component[5])/(dx*dx)

                diff_data[-2] = (1.0/12.0*data_component[-6] - 1.0/2.0*data_component[-5] \
                                 + 7.0/6.0*data_component[-4] - 1.0/3.0*data_component[-3] \
                                 - 5.0/4.0*data_component[-2] + 5.0/6.0*data_component[-1])/(dx*dx)

                diff_data[-1] = (-5.0/6.0*data_component[-6] + 61.0/12.0*data_component[-5] \
                                 - 13.0*data_component[-4] + 107.0/6.0*data_component[-3] \
                                 - 77.0/6.0*data_component[-2] + 15.0/4.0*data_component[-1])/(dx*dx)

            elif dim == 2:
                diff_data[0, :] = (15.0/4.0*data_component[0, :] - 77.0/6.0*data_component[1, :] \
                                   + 107.0/6.0*data_component[2, :] - 13.0*data_component[3, :] \
                                   + 61.0/12.0*data_component[4, :] - 5.0/6.0*data_component[5, :])/(dx*dx)

                diff_data[1, :] = (5.0/6.0*data_component[0, :] - 5.0/4.0*data_component[1, :] \
                                   - 1.0/3.0*data_component[2, :] + 7.0/6.0*data_component[3, :] \
                                   - 1.0/2.0*data_component[4, :] + 1.0/12.0*data_component[5, :])/(dx*dx)

                diff_data[-2, :] = (1.0/12.0*data_component[-6, :] - 1.0/2.0*data_component[-5, :] \
                                    + 7.0/6.0*data_component[-4, :] - 1.0/3.0*data_component[-3, :] \
                                    - 5.0/4.0*data_component[-2, :] + 5.0/6.0*data_component[-1, :])/(dx*dx)

                diff_data[-1, :] = (-5.0/6.0*data_component[-6, :] + 61.0/12.0*data_component[-5, :] \
                                    - 13.0*data_component[-4, :] + 107.0/6.0*data_component[-3, :] \
                                    - 77.0/6.0*data_component[-2, :] + 15.0/4.0*data_component[-1, :])/(dx*dx)

            elif dim == 3:
                diff_data[0, :, :] = (15.0/4.0*data_component[0, :, :] - 77.0/6.0*data_component[1, :, :] \
                                      + 107.0/6.0*data_component[2, :, :] - 13.0*data_component[3, :, :] \
                                      + 61.0/12.0*data_component[4, :, :] - 5.0/6.0*data_component[5, :, :])/(dx*dx)

                diff_data[1, :, :] = (5.0/6.0*data_component[0, :, :] - 5.0/4.0*data_component[1, :, :] \
                                      - 1.0/3.0*data_component[2, :, :] + 7.0/6.0*data_component[3, :, :] \
                                      - 1.0/2.0*data_component[4, :, :] + 1.0/12.0*data_component[5, :, :])/(dx*dx)

                diff_data[-2, :, :] = (1.0/12.0*data_component[-6, :, :] - 1.0/2.0*data_component[-5, :, :] \
                                       + 7.0/6.0*data_component[-4, :, :] - 1.0/3.0*data_component[-3, :, :] \
                                       - 5.0/4.0*data_component[-2, :, :] + 5.0/6.0*data_component[-1, :, :])/(dx*dx)

                diff_data[-1, :, :] = (-5.0/6.0*data_component[-6, :, :] + 61.0/12.0*data_component[-5, :, :] \
                                       - 13.0*data_component[-4, :, :] + 107.0/6.0*data_component[-3, :, :] \
                                       - 77.0/6.0*data_component[-2, :, :] + 15.0/4.0*data_component[-1, :, :])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

        elif direction == 1:
            if dim < 2:
                raise RuntimeError('There is no second direction in data with less than two dimensions!')

            elif dim == 2:
                diff_data[:, 0] = (15.0/4.0*data_component[:, 0] - 77.0/6.0*data_component[:, 1] \
                                   + 107.0/6.0*data_component[:, 2] - 13.0*data_component[:, 3] \
                                   + 61.0/12.0*data_component[:, 4] - 5.0/6.0*data_component[:, 5])/(dx*dx)

                diff_data[:, 1] = (5.0/6.0*data_component[:, 0] - 5.0/4.0*data_component[:, 1] \
                                   - 1.0/3.0*data_component[:, 2] + 7.0/6.0*data_component[:, 3] \
                                   - 1.0/2.0*data_component[:, 4] + 1.0/12.0*data_component[:, 5])/(dx*dx)

                diff_data[:, -2] = (1.0/12.0*data_component[:, -6] - 1.0/2.0*data_component[:, -5] \
                                    + 7.0/6.0*data_component[:, -4] - 1.0/3.0*data_component[:, -3] \
                                    - 5.0/4.0*data_component[:, -2] + 5.0/6.0*data_component[:, -1])/(dx*dx)

                diff_data[:, -1] = (-5.0/6.0*data_component[:, -6] + 61.0/12.0*data_component[:, -5] \
                                    - 13.0*data_component[:, -4] + 107.0/6.0*data_component[:, -3] \
                                    - 77.0/6.0*data_component[:, -2] + 15.0/4.0*data_component[:, -1])/(dx*dx)

            elif dim == 3:
                diff_data[:, 0, :] = (15.0/4.0*data_component[:, 0, :] - 77.0/6.0*data_component[:, 1, :] \
                                      + 107.0/6.0*data_component[:, 2, :] - 13.0*data_component[:, 3, :] \
                                      + 61.0/12.0*data_component[:, 4, :] - 5.0/6.0*data_component[:, 5, :])/(dx*dx)

                diff_data[:, 1, :] = (5.0/6.0*data_component[:, 0, :] - 5.0/4.0*data_component[:, 1, :] \
                                      - 1.0/3.0*data_component[:, 2, :] + 7.0/6.0*data_component[:, 3, :] \
                                      - 1.0/2.0*data_component[:, 4, :] + 1.0/12.0*data_component[:, 5, :])/(dx*dx)

                diff_data[:, -2, :] = (1.0/12.0*data_component[:, -6, :] - 1.0/2.0*data_component[:, -5, :] \
                                       + 7.0/6.0*data_component[:, -4, :] - 1.0/3.0*data_component[:, -3, :] \
                                       - 5.0/4.0*data_component[:, -2, :] + 5.0/6.0*data_component[:, -1, :])/(dx*dx)

                diff_data[:, -1, :] = (-5.0/6.0*data_component[:, -6, :] + 61.0/12.0*data_component[:, -5, :] \
                                       - 13.0*data_component[:, -4, :] + 107.0/6.0*data_component[:, -3, :] \
                                       - 77.0/6.0*data_component[:, -2, :] + 15.0/4.0*data_component[:, -1, :])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

        elif direction == 2:
            if dim < 3:
                raise IOError('There is no third direction in data with less than three dimensions!')

            elif dim == 3:
                diff_data[:, :, 0] = (15.0/4.0*data_component[:, :, 0] - 77.0/6.0*data_component[:, :, 1] \
                                      + 107.0/6.0*data_component[:, :, 2] - 13.0*data_component[:, :, 3] \
                                      + 61.0/12.0*data_component[:, :, 4] - 5.0/6.0*data_component[:, :, 5])/(dx*dx)

                diff_data[:, :, 1] = (5.0/6.0*data_component[:, :, 0] - 5.0/4.0*data_component[:, :, 1] \
                                      - 1.0/3.0*data_component[:, :, 2] + 7.0/6.0*data_component[:, :, 3] \
                                      - 1.0/2.0*data_component[:, :, 4] + 1.0/12.0*data_component[:, :, 5])/(dx*dx)

                diff_data[:, :, -2] = (1.0/12.0*data_component[:, :, -6] - 1.0/2.0*data_component[:, :, -5] \
                                       + 7.0/6.0*data_component[:, :, -4] - 1.0/3.0*data_component[:, :, -3] \
                                       - 5.0/4.0*data_component[:, :, -2] + 5.0/6.0*data_component[:, :, -1])/(dx*dx)

                diff_data[:, :, -1] = (-5.0/6.0*data_component[:, :, -6] + 61.0/12.0*data_component[:, :, -5] \
                                       - 13.0*data_component[:, :, -4] + 107.0/6.0*data_component[:, :, -3] \
                                       - 77.0/6.0*data_component[:, :, -2] + 15.0/4.0*data_component[:, :, -1])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

    return diff_data


def differentiateSixthOrderFiniteDifference(data, dx, direction, component_idx, use_one_sided, dimension, data_order):
    """
    Compute second order derivative using explicit sixth order finite differencing.
    """

    # Get the shape of data.

    data_shape = numpy.array(data.shape)

    # Check whether the direction is valid.

    if direction < 0 or direction > 2:
        raise RuntimeError('Direction < 0 or > 2 is invalid!')

    # Check whether the shape of data is valid.

    if component_idx is None:
        if data.ndim < 1 or data.ndim > 3:
            raise RuntimeError('Shape of data is invalid!')

    else:
        if data.ndim < 2 or data.ndim > 4:
            raise RuntimeError('Shape of data is invalid!')

        # Check whether the component_idx is valid and get the shape of the component's data.

        if data_order == 'C':
            if component_idx >= data.shape[0] or component_idx < 0:
                raise RuntimeError('Component index is invalid!')

            data_shape = numpy.array(data_shape[1:])

        else:
            if component_idx >= data.shape[-1] or component_idx < 0:
                raise RuntimeError('Component index is invalid!')

            data_shape = numpy.array(data_shape[:-1])

    # Get the dimension of data.

    dim = dimension

    # Check whether data size is large enough for sixth order second derivative.

    if use_one_sided == True:
        if direction == 0:
            if data_shape[0] < 8:
                raise RuntimeError('First dimension of data is not large enough!')

        elif direction == 1:
            if data_shape[1] < 8:
                raise RuntimeError('Second dimension of data is not large enough!')

        elif direction == 2:
            if data_shape[2] < 8:
                raise RuntimeError('Third dimension of data is not large enough!')

    else:
        if direction == 0:
            if data_shape[0] < 7:
                raise RuntimeError('First dimension of data is not large enough!')

        elif direction == 1:
            if data_shape[1] < 7:
                raise RuntimeError('Second dimension of data is not large enough!')

        elif direction == 2:
            if data_shape[2] < 7:
                raise RuntimeError('Third dimension of data is not large enough!')

    # Initialize container to store the derivatives. The elements in the container
    # are initialized as NAN values.

    diff_data = numpy.empty(data_shape, dtype=data.dtype, order=data_order)
    diff_data[:] = numpy.NAN

    # Get the component's data.

    data_component = None

    if component_idx is None:
        data_component = data
    else:
        if data_order == 'C':
            if dim == 1:
                data_component = data[component_idx, :]
            elif dim == 2:
                data_component = data[component_idx, :, :]
            elif dim == 3:
                data_component = data[component_idx, :, :, :]
        else:
            if dim == 1:
                data_component = data[:, component_idx]
            elif dim == 2:
                data_component = data[:, :, component_idx]
            elif dim == 3:
                data_component = data[:, :, :, component_idx]

    # Compute the derivatives in the interior of the domain.

    if direction == 0:
        if dim == 1:
            diff_data[3:-3] = (1.0/90.0*data_component[0:-6] - 3.0/20.0*data_component[1:-5] \
                               + 3.0/2.0*data_component[2:-4] - 49.0/18.0*data_component[3:-3] \
                               + 3.0/2.0*data_component[4:-2] - 3.0/20.0*data_component[5:-1] \
                               + 1.0/90.0*data_component[6:])/(dx*dx)

        elif dim == 2:
            diff_data[3:-3, :] = (1.0/90.0*data_component[0:-6, :] - 3.0/20.0*data_component[1:-5, :] \
                                  + 3.0/2.0*data_component[2:-4, :] - 49.0/18.0*data_component[3:-3, :] \
                                  + 3.0/2.0*data_component[4:-2, :] - 3.0/20.0*data_component[5:-1, :] \
                                  + 1.0/90.0*data_component[6:, :])/(dx*dx)

        elif dim == 3:
            diff_data[3:-3, :, :] = (1.0/90.0*data_component[0:-6, :, :] - 3.0/20.0*data_component[1:-5, :, :] \
                                     + 3.0/2.0*data_component[2:-4, :, :] - 49.0/18.0*data_component[3:-3, :, :] \
                                     + 3.0/2.0*data_component[4:-2, :, :] - 3.0/20.0*data_component[5:-1, :, :] \
                                     + 1.0/90.0*data_component[6:, :, :])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    elif direction == 1:
        if dim < 2:
            raise IOError('There is no second direction in data with less than two dimensions!')

        elif dim == 2:
            diff_data[:, 3:-3] = (1.0/90.0*data_component[:, 0:-6] - 3.0/20.0*data_component[:, 1:-5] \
                                  + 3.0/2.0*data_component[:, 2:-4] - 49.0/18.0*data_component[:, 3:-3] \
                                  + 3.0/2.0*data_component[:, 4:-2] - 3.0/20.0*data_component[:, 5:-1] \
                                  + 1.0/90.0*data_component[:, 6:])/(dx*dx)

        elif dim == 3:
            diff_data[:, 3:-3, :] = (1.0/90.0*data_component[:, 0:-6, :] - 3.0/20.0*data_component[:, 1:-5, :] \
                                     + 3.0/2.0*data_component[:, 2:-4, :] - 49.0/18.0*data_component[:, 3:-3, :] \
                                     + 3.0/2.0*data_component[:, 4:-2, :] - 3.0/20.0*data_component[:, 5:-1, :] \
                                     + 1.0/90.0*data_component[:, 6:, :])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    elif direction == 2:
        if dim < 3:
            raise IOError('There is no third direction in data with less than three dimensions!')

        elif dim == 3:
            diff_data[:, :, 3:-3] = (1.0/90.0*data_component[:, :, 0:-6] - 3.0/20.0*data_component[:, :, 1:-5] \
                                  + 3.0/2.0*data_component[:, :, 2:-4] - 49.0/18.0*data_component[:, :, 3:-3] \
                                  + 3.0/2.0*data_component[:, :, 4:-2] - 3.0/20.0*data_component[:, :, 5:-1] \
                                  + 1.0/90.0*data_component[:, :, 6:])/(dx*dx)

        else:
            raise RuntimeError('Data dimension > 3 not supported!')

    # Compute the derivatives at the boundaries.

    if use_one_sided == True:
        if direction == 0:
            if dim == 1:
                diff_data[0] = (469.0/90.0*data_component[0] - 223.0/10.0*data_component[1] \
                                + 879.0/20.0*data_component[2] - 949.0/18.0*data_component[3] \
                                + 41.0*data_component[4] - 201.0/10.0*data_component[5] \
                                + 1019.0/180.0*data_component[6] - 7.0/10.0*data_component[7])/(dx*dx)

                diff_data[1] = (7.0/10.0*data_component[0] - 7.0/18.0*data_component[1] \
                                - 27.0/10.0*data_component[2] + 19.0/4.0*data_component[3] \
                                - 67.0/18.0*data_component[4] + 9.0/5.0*data_component[5] \
                                - 1.0/2.0*data_component[6] + 11.0/180.0*data_component[7])/(dx*dx)

                diff_data[2] = (-11.0/180.0*data_component[0] + 107.0/90.0*data_component[1] \
                                - 21.0/10.0*data_component[2] + 13.0/18.0*data_component[3] \
                                + 17.0/36.0*data_component[4] - 3.0/10.0*data_component[5] \
                                + 4.0/45.0*data_component[6] - 1.0/90.0*data_component[7])/(dx*dx)

                diff_data[-3] = (-1.0/90.0*data_component[-8] + 4.0/45.0*data_component[-7] \
                                 - 3.0/10.0*data_component[-6] + 17.0/36.0*data_component[-5] \
                                 + 13.0/18.0*data_component[-4] - 21.0/10.0*data_component[-3] \
                                 + 107.0/90.0*data_component[-2] - 11.0/180.0*data_component[-1])/(dx*dx)

                diff_data[-2] = (11.0/180.0*data_component[-8] - 1.0/2.0*data_component[-7] \
                                 + 9.0/5.0*data_component[-6] - 67.0/18.0*data_component[-5] \
                                 + 19.0/4.0*data_component[-4] - 27.0/10.0*data_component[-3] \
                                 - 7.0/18.0*data_component[-2] + 7.0/10.0*data_component[-1])/(dx*dx)

                diff_data[-1] = (-7.0/10.0*data_component[-8] + 1019.0/180.0*data_component[-7] \
                                 - 201.0/10.0*data_component[-6] + 41.0*data_component[-5] \
                                 - 949.0/18.0*data_component[-4] + 879.0/20.0*data_component[-3] \
                                 - 223.0/10.0*data_component[-2] + 469.0/90.0*data_component[-1])/(dx*dx)

            elif dim == 2:
                diff_data[0, :] = (469.0/90.0*data_component[0, :] - 223.0/10.0*data_component[1, :] \
                                   + 879.0/20.0*data_component[2, :] - 949.0/18.0*data_component[3, :] \
                                   + 41.0*data_component[4, :] - 201.0/10.0*data_component[5, :] \
                                   + 1019.0/180.0*data_component[6, :] - 7.0/10.0*data_component[7, :])/(dx*dx)

                diff_data[1, :] = (7.0/10.0*data_component[0, :] - 7.0/18.0*data_component[1, :] \
                                   - 27.0/10.0*data_component[2, :] + 19.0/4.0*data_component[3, :] \
                                   - 67.0/18.0*data_component[4, :] + 9.0/5.0*data_component[5, :] \
                                   - 1.0/2.0*data_component[6, :] + 11.0/180.0*data_component[7, :])/(dx*dx)

                diff_data[2, :] = (-11.0/180.0*data_component[0, :] + 107.0/90.0*data_component[1, :] \
                                   - 21.0/10.0*data_component[2, :] + 13.0/18.0*data_component[3, :] \
                                   + 17.0/36.0*data_component[4, :] - 3.0/10.0*data_component[5, :] \
                                   + 4.0/45.0*data_component[6, :] - 1.0/90.0*data_component[7, :])/(dx*dx)

                diff_data[-3, :] = (-1.0/90.0*data_component[-8, :] + 4.0/45.0*data_component[-7, :] \
                                    - 3.0/10.0*data_component[-6, :] + 17.0/36.0*data_component[-5, :] \
                                    + 13.0/18.0*data_component[-4, :] - 21.0/10.0*data_component[-3, :] \
                                    + 107.0/90.0*data_component[-2, :] - 11.0/180.0*data_component[-1, :])/(dx*dx)

                diff_data[-2, :] = (11.0/180.0*data_component[-8, :] - 1.0/2.0*data_component[-7, :] \
                                    + 9.0/5.0*data_component[-6, :] - 67.0/18.0*data_component[-5, :] \
                                    + 19.0/4.0*data_component[-4, :] - 27.0/10.0*data_component[-3, :] \
                                    - 7.0/18.0*data_component[-2, :] + 7.0/10.0*data_component[-1, :])/(dx*dx)

                diff_data[-1, :] = (-7.0/10.0*data_component[-8, :] + 1019.0/180.0*data_component[-7, :] \
                                    - 201.0/10.0*data_component[-6, :] + 41.0*data_component[-5, :] \
                                    - 949.0/18.0*data_component[-4, :] + 879.0/20.0*data_component[-3, :] \
                                    - 223.0/10.0*data_component[-2, :] + 469.0/90.0*data_component[-1, :])/(dx*dx)

            elif dim == 3:
                diff_data[0, :, :] = (469.0/90.0*data_component[0, :, :] - 223.0/10.0*data_component[1, :, :] \
                                      + 879.0/20.0*data_component[2, :, :] - 949.0/18.0*data_component[3, :, :] \
                                      + 41.0*data_component[4, :, :] - 201.0/10.0*data_component[5, :, :] \
                                      + 1019.0/180.0*data_component[6, :, :] - 7.0/10.0*data_component[7, :, :])/(dx*dx)

                diff_data[1, :, :] = (7.0/10.0*data_component[0, :, :] - 7.0/18.0*data_component[1, :, :] \
                                      - 27.0/10.0*data_component[2, :, :] + 19.0/4.0*data_component[3, :, :] \
                                      - 67.0/18.0*data_component[4, :, :] + 9.0/5.0*data_component[5, :, :] \
                                      - 1.0/2.0*data_component[6, :, :] + 11.0/180.0*data_component[7, :, :])/(dx*dx)

                diff_data[2, :, :] = (-11.0/180.0*data_component[0, :, :] + 107.0/90.0*data_component[1, :, :] \
                                      - 21.0/10.0*data_component[2, :, :] + 13.0/18.0*data_component[3, :, :] \
                                      + 17.0/36.0*data_component[4, :, :] - 3.0/10.0*data_component[5, :, :] \
                                      + 4.0/45.0*data_component[6, :, :] - 1.0/90.0*data_component[7, :, :])/(dx*dx)

                diff_data[-3, :, :] = (-1.0/90.0*data_component[-8, :, :] + 4.0/45.0*data_component[-7, :, :] \
                                       - 3.0/10.0*data_component[-6, :, :] + 17.0/36.0*data_component[-5, :, :] \
                                       + 13.0/18.0*data_component[-4, :, :] - 21.0/10.0*data_component[-3, :, :] \
                                       + 107.0/90.0*data_component[-2, :, :] - 11.0/180.0*data_component[-1, :, :])/(dx*dx)

                diff_data[-2, :, :] = (11.0/180.0*data_component[-8, :, :] - 1.0/2.0*data_component[-7, :, :] \
                                       + 9.0/5.0*data_component[-6, :, :] - 67.0/18.0*data_component[-5, :, :] \
                                       + 19.0/4.0*data_component[-4, :, :] - 27.0/10.0*data_component[-3, :, :] \
                                       - 7.0/18.0*data_component[-2, :, :] + 7.0/10.0*data_component[-1, :, :])/(dx*dx)

                diff_data[-1, :, :] = (-7.0/10.0*data_component[-8, :, :] + 1019.0/180.0*data_component[-7, :, :] \
                                       - 201.0/10.0*data_component[-6, :, :] + 41.0*data_component[-5, :, :] \
                                       - 949.0/18.0*data_component[-4, :, :] + 879.0/20.0*data_component[-3, :, :] \
                                       - 223.0/10.0*data_component[-2, :, :] + 469.0/90.0*data_component[-1, :, :])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

        elif direction == 1:
            if dim < 2:
                raise RuntimeError('There is no second direction in data with less than two dimensions!')

            elif dim == 2:
                diff_data[:, 0] = (469.0/90.0*data_component[:, 0] - 223.0/10.0*data_component[:, 1] \
                                   + 879.0/20.0*data_component[:, 2] - 949.0/18.0*data_component[:, 3] \
                                   + 41.0*data_component[:, 4] - 201.0/10.0*data_component[:, 5] \
                                   + 1019.0/180.0*data_component[:, 6] - 7.0/10.0*data_component[:, 7])/(dx*dx)

                diff_data[:, 1] = (7.0/10.0*data_component[:, 0] - 7.0/18.0*data_component[:, 1] \
                                   - 27.0/10.0*data_component[:, 2] + 19.0/4.0*data_component[:, 3] \
                                   - 67.0/18.0*data_component[:, 4] + 9.0/5.0*data_component[:, 5] \
                                   - 1.0/2.0*data_component[:, 6] + 11.0/180.0*data_component[:, 7])/(dx*dx)

                diff_data[:, 2] = (-11.0/180.0*data_component[:, 0] + 107.0/90.0*data_component[:, 1] \
                                   - 21.0/10.0*data_component[:, 2] + 13.0/18.0*data_component[:, 3] \
                                   + 17.0/36.0*data_component[:, 4] - 3.0/10.0*data_component[:, 5] \
                                   + 4.0/45.0*data_component[:, 6] - 1.0/90.0*data_component[:, 7])/(dx*dx)

                diff_data[:, -3] = (-1.0/90.0*data_component[:, -8] + 4.0/45.0*data_component[:, -7] \
                                    - 3.0/10.0*data_component[:, -6] + 17.0/36.0*data_component[:, -5] \
                                    + 13.0/18.0*data_component[:, -4] - 21.0/10.0*data_component[:, -3] \
                                    + 107.0/90.0*data_component[:, -2] - 11.0/180.0*data_component[:, -1])/(dx*dx)

                diff_data[:, -2] = (11.0/180.0*data_component[:, -8] - 1.0/2.0*data_component[:, -7] \
                                    + 9.0/5.0*data_component[:, -6] - 67.0/18.0*data_component[:, -5] \
                                    + 19.0/4.0*data_component[:, -4] - 27.0/10.0*data_component[:, -3] \
                                    - 7.0/18.0*data_component[:, -2] + 7.0/10.0*data_component[:, -1])/(dx*dx)

                diff_data[:, -1] = (-7.0/10.0*data_component[:, -8] + 1019.0/180.0*data_component[:, -7] \
                                    - 201.0/10.0*data_component[:, -6] + 41.0*data_component[:, -5] \
                                    - 949.0/18.0*data_component[:, -4] + 879.0/20.0*data_component[:, -3] \
                                    - 223.0/10.0*data_component[:, -2] + 469.0/90.0*data_component[:, -1])/(dx*dx)

            elif dim == 3:
                diff_data[:, 0, :] = (469.0/90.0*data_component[:, 0, :] - 223.0/10.0*data_component[:, 1, :] \
                                      + 879.0/20.0*data_component[:, 2, :] - 949.0/18.0*data_component[:, 3, :] \
                                      + 41.0*data_component[:, 4, :] - 201.0/10.0*data_component[:, 5, :] \
                                      + 1019.0/180.0*data_component[:, 6, :] - 7.0/10.0*data_component[:, 7, :])/(dx*dx)

                diff_data[:, 1, :] = (7.0/10.0*data_component[:, 0, :] - 7.0/18.0*data_component[:, 1, :] \
                                      - 27.0/10.0*data_component[:, 2, :] + 19.0/4.0*data_component[:, 3, :] \
                                      - 67.0/18.0*data_component[:, 4, :] + 9.0/5.0*data_component[:, 5, :] \
                                      - 1.0/2.0*data_component[:, 6, :] + 11.0/180.0*data_component[:, 7, :])/(dx*dx)

                diff_data[:, 2, :] = (-11.0/180.0*data_component[:, 0, :] + 107.0/90.0*data_component[:, 1, :] \
                                      - 21.0/10.0*data_component[:, 2, :] + 13.0/18.0*data_component[:, 3, :] \
                                      + 17.0/36.0*data_component[:, 4, :] - 3.0/10.0*data_component[:, 5, :] \
                                      + 4.0/45.0*data_component[:, 6, :] - 1.0/90.0*data_component[:, 7, :])/(dx*dx)

                diff_data[:, -3, :] = (-1.0/90.0*data_component[:, -8, :] + 4.0/45.0*data_component[:, -7, :] \
                                       - 3.0/10.0*data_component[:, -6, :] + 17.0/36.0*data_component[:, -5, :] \
                                       + 13.0/18.0*data_component[:, -4, :] - 21.0/10.0*data_component[:, -3, :] \
                                       + 107.0/90.0*data_component[:, -2, :] - 11.0/180.0*data_component[:, -1, :])/(dx*dx)

                diff_data[:, -2, :] = (11.0/180.0*data_component[:, -8, :] - 1.0/2.0*data_component[:, -7, :] \
                                       + 9.0/5.0*data_component[:, -6, :] - 67.0/18.0*data_component[:, -5, :] \
                                       + 19.0/4.0*data_component[:, -4, :] - 27.0/10.0*data_component[:, -3, :] \
                                       - 7.0/18.0*data_component[:, -2, :] + 7.0/10.0*data_component[:, -1, :])/(dx*dx)

                diff_data[:, -1, :] = (-7.0/10.0*data_component[:, -8, :] + 1019.0/180.0*data_component[:, -7, :] \
                                       - 201.0/10.0*data_component[:, -6, :] + 41.0*data_component[:, -5, :] \
                                       - 949.0/18.0*data_component[:, -4, :] + 879.0/20.0*data_component[:, -3, :] \
                                       - 223.0/10.0*data_component[:, -2, :] + 469.0/90.0*data_component[:, -1, :])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

        elif direction == 2:
            if dim < 3:
                raise IOError('There is no third direction in data with less than three dimensions!')

            elif dim == 3:
                diff_data[:, :, 0] = (469.0/90.0*data_component[:, :, 0] - 223.0/10.0*data_component[:, :, 1] \
                                      + 879.0/20.0*data_component[:, :, 2] - 949.0/18.0*data_component[:, :, 3] \
                                      + 41.0*data_component[:, :, 4] - 201.0/10.0*data_component[:, :, 5] \
                                      + 1019.0/180.0*data_component[:, :, 6] - 7.0/10.0*data_component[:, :, 7])/(dx*dx)

                diff_data[:, :, 1] = (7.0/10.0*data_component[:, :, 0] - 7.0/18.0*data_component[:, :, 1] \
                                      - 27.0/10.0*data_component[:, :, 2] + 19.0/4.0*data_component[:, :, 3] \
                                      - 67.0/18.0*data_component[:, :, 4] + 9.0/5.0*data_component[:, :, 5] \
                                      - 1.0/2.0*data_component[:, :, 6] + 11.0/180.0*data_component[:, :, 7])/(dx*dx)

                diff_data[:, :, 2] = (-11.0/180.0*data_component[:, :, 0] + 107.0/90.0*data_component[:, :, 1] \
                                      - 21.0/10.0*data_component[:, :, 2] + 13.0/18.0*data_component[:, :, 3] \
                                      + 17.0/36.0*data_component[:, :, 4] - 3.0/10.0*data_component[:, :, 5] \
                                      + 4.0/45.0*data_component[:, :, 6] - 1.0/90.0*data_component[:, :, 7])/(dx*dx)

                diff_data[:, :, -3] = (-1.0/90.0*data_component[:, :, -8] + 4.0/45.0*data_component[:, :, -7] \
                                       - 3.0/10.0*data_component[:, :, -6] + 17.0/36.0*data_component[:, :, -5] \
                                       + 13.0/18.0*data_component[:, :, -4] - 21.0/10.0*data_component[:, :, -3] \
                                       + 107.0/90.0*data_component[:, :, -2] - 11.0/180.0*data_component[:, :, -1])/(dx*dx)

                diff_data[:, :, -2] = (11.0/180.0*data_component[:, :, -8] - 1.0/2.0*data_component[:, :, -7] \
                                       + 9.0/5.0*data_component[:, :, -6] - 67.0/18.0*data_component[:, :, -5] \
                                       + 19.0/4.0*data_component[:, :, -4] - 27.0/10.0*data_component[:, :, -3] \
                                       - 7.0/18.0*data_component[:, :, -2] + 7.0/10.0*data_component[:, :, -1])/(dx*dx)

                diff_data[:, :, -1] = (-7.0/10.0*data_component[:, :, -8] + 1019.0/180.0*data_component[:, :, -7] \
                                       - 201.0/10.0*data_component[:, :, -6] + 41.0*data_component[:, :, -5] \
                                       - 949.0/18.0*data_component[:, :, -4] + 879.0/20.0*data_component[:, :, -3] \
                                       - 223.0/10.0*data_component[:, :, -2] + 469.0/90.0*data_component[:, :, -1])/(dx*dx)

            else:
                raise RuntimeError('Data dimension > 3 not supported!')

    return diff_data
