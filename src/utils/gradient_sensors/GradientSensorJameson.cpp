#include "utils/gradient_sensors/GradientSensorJameson.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#define EPSILON 1e-40

GradientSensorJameson::GradientSensorJameson(
    const std::string& object_name,
    const tbox::Dimension& dim):
        GradientSensor(
            object_name,
            dim)
{
    d_num_gradient_ghosts = hier::IntVector::getOne(d_dim);
}


/*
 * Compute the gradient with the given cell data.
 */
boost::shared_ptr<pdat::CellData<double> >
GradientSensorJameson::ComputeGradient(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<double> > cell_data)
{
    // Get the depth of the given cell data.
    int data_depth = cell_data->getDepth();
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box dummy_box = patch.getBox();
    const hier::Box interior_box = dummy_box;
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of Patch plus
    // ghost cells.
    dummy_box.grow(d_num_ghosts);
    const hier::Box ghost_box = dummy_box;
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    // Allocate sensor values.
    boost::shared_ptr<pdat::CellData<double> > sensor_value(
        new pdat::CellData<double>(interior_box, data_depth, d_num_ghosts));
    
    for (int di = 0; di < data_depth; di++)
    {
        // Get the pointer to the current depth component of the given cell data.
        double* f = cell_data->getPointer(di);
        
        // Get the pointer to the current depth component of the sensor cell data.
        double* psi = sensor_value->getPointer(di);
        
        if (d_dim == tbox::Dimension(1))
        {
            // NOT YET IMPLEMENTED
        }
        else if (d_dim == tbox::Dimension(2))
        {
            // Allocate sensor values in different dimensions.
            boost::shared_ptr<pdat::CellData<double> > sensor_value_x(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
            boost::shared_ptr<pdat::CellData<double> > sensor_value_y(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
            
            double* psi_x = sensor_value_x->getPointer(0);
            double* psi_y = sensor_value_y->getPointer(0);
            
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*ghostcell_dims[0];
                    
                    const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*ghostcell_dims[0];
                    
                    const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*ghostcell_dims[0];
                    
                    psi_x[idx] = fabs(f[idx_x_R] - 2*f[idx] + f[idx_x_L])/
                        (f[idx_x_R] + 2*f[idx] + f[idx_x_L] + EPSILON);
                }
            }
            
            for (int i = 0; i < interior_dims[0]; i++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    // Compute indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*ghostcell_dims[0];
                    
                    const int idx_y_B = (i + d_num_ghosts[0]) +
                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                    
                    const int idx_y_T = (i + d_num_ghosts[0]) +
                        (j + 1 + d_num_ghosts[1])*ghostcell_dims[0];
                    
                    psi_y[idx] = fabs(f[idx_y_T] - 2*f[idx] + f[idx_y_B])/
                        (f[idx_y_T] + 2*f[idx] + f[idx_y_B] + EPSILON);
                }
            }
            
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute indices.
                    const int idx = (i + d_num_ghosts[0]) +
                        (j + d_num_ghosts[1])*ghostcell_dims[0];
                    
                    psi[idx] = fmax(psi_x[idx], psi_y[idx]);
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            // Allocate sensor values in different dimensions.
            boost::shared_ptr<pdat::CellData<double> > sensor_value_x(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
            boost::shared_ptr<pdat::CellData<double> > sensor_value_y(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
            boost::shared_ptr<pdat::CellData<double> > sensor_value_z(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
            
            double* psi_x = sensor_value_x->getPointer(0);
            double* psi_y = sensor_value_y->getPointer(0);
            double* psi_z = sensor_value_z->getPointer(0);
            
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        const int idx_x_L = (i - 1 + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        const int idx_x_R = (i + 1 + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        psi_x[idx] = fabs(f[idx_x_R] - 2*f[idx] + f[idx_x_L])/
                            (f[idx_x_R] + 2*f[idx] + f[idx_x_L] + EPSILON);
                    }
                }
            }
            
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int i = 0; i < interior_dims[1]; i++)
                {
                    for (int j = 0; j < interior_dims[0]; j++)
                    {
                        // Compute indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        const int idx_y_B = (i + d_num_ghosts[0]) +
                            (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        const int idx_y_T = (i + d_num_ghosts[0]) +
                            (j + 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        psi_y[idx] = fabs(f[idx_y_T] - 2*f[idx] + f[idx_y_B])/
                            (f[idx_y_T] + 2*f[idx] + f[idx_y_B] + EPSILON);
                    }
                }
            }
            
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        // Compute indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        const int idx_z_B = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        const int idx_z_F = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        psi_z[idx] = fabs(f[idx_z_F] - 2*f[idx] + f[idx_z_B])/
                            (f[idx_z_F] + 2*f[idx] + f[idx_z_B] + EPSILON);
                    }
                }
            }
            
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute indices.
                        const int idx = (i + d_num_ghosts[0]) +
                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                        
                        psi[idx] = fmax(fmax(psi_x[idx], psi_y[idx]), psi_z[idx]);
                    }
                }
            }
        }
    }
    
    return sensor_value;
}
