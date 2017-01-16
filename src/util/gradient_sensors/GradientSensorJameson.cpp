#include "util/gradient_sensors/GradientSensorJameson.hpp"

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
void
GradientSensorJameson::computeGradient(
    hier::Patch& patch,
    boost::shared_ptr<pdat::CellData<double> > cell_data,
    boost::shared_ptr<pdat::CellData<double> > gradient,
    int depth)
{
    if (cell_data->getGhostCellWidth() < d_num_gradient_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The ghost cell width is smaller than required."
            << std::endl);
    }
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    // Get the number of ghost cells of the cell data and gradient data.
    const hier::IntVector num_ghosts_cell_data = cell_data->getGhostCellWidth();
    const hier::IntVector num_ghosts_gradient = gradient->getGhostCellWidth();
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus ghost cells.
    const hier::Box ghost_box_cell_data = cell_data->getGhostBox();
    const hier::IntVector ghostcell_dims_cell_data = ghost_box_cell_data.numberCells();
    
    const hier::Box ghost_box_gradient = gradient->getGhostBox();
    const hier::IntVector ghostcell_dims_gradient = ghost_box_gradient.numberCells();
    
    // Get the pointer to the current depth component of the given cell data.
    double* f = cell_data->getPointer(depth);
    
    // Get the pointer to the gradient.
    double* psi = gradient->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        // Allocate memory.
        boost::shared_ptr<pdat::CellData<double> > gradient_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        
        double* psi_x = gradient_x->getPointer(0);
        double* mean_x = local_mean_value_x->getPointer(0);
        
        for (int i = 0; i < interior_dims[0]; i++)
        {
            // Compute the linear indices.
            const int idx = i + num_ghosts_gradient[0];
            const int idx_x_L = i - 1 + num_ghosts_cell_data[0];
            const int idx_x   = i + num_ghosts_cell_data[0];
            const int idx_x_R = i + 1 + num_ghosts_cell_data[0];
            
            psi_x[idx] = f[idx_x_R] - 2*f[idx_x] + f[idx_x_L];
            mean_x[idx] = f[idx_x_R] + 2*f[idx_x] + f[idx_x_L];
        }
        
        for (int i = 0; i < interior_dims[0]; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_gradient[0];
            
            psi[idx] = fabs(psi_x[idx])/(mean_x[idx] + EPSILON);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        // Allocate memory in different dimensions.
        boost::shared_ptr<pdat::CellData<double> > gradient_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > gradient_y(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_y(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        
        double* psi_x = gradient_x->getPointer(0);
        double* psi_y = gradient_y->getPointer(0);
        double* mean_x = local_mean_value_x->getPointer(0);
        double* mean_y = local_mean_value_y->getPointer(0);
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_gradient[0]) +
                    (j + num_ghosts_gradient[1])*ghostcell_dims_gradient[0];
                
                const int idx_x_L = (i - 1 + num_ghosts_cell_data[0]) +
                    (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                
                const int idx_x   = (i + num_ghosts_cell_data[0]) +
                    (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                
                const int idx_x_R = (i + 1 + num_ghosts_cell_data[0]) +
                    (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                
                psi_x[idx] = f[idx_x_R] - 2*f[idx_x] + f[idx_x_L];
                mean_x[idx] = f[idx_x_R] + 2*f[idx_x] + f[idx_x_L];
            }
        }
        
        for (int i = 0; i < interior_dims[0]; i++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                // Compute the linear indices.
                const int idx = (i + num_ghosts_gradient[0]) +
                    (j + num_ghosts_gradient[1])*ghostcell_dims_gradient[0];
                
                const int idx_y_B = (i + num_ghosts_cell_data[0]) +
                    (j - 1 + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                
                const int idx_y   = (i + num_ghosts_cell_data[0]) +
                    (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                
                const int idx_y_T = (i + num_ghosts_cell_data[0]) +
                    (j + 1 + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0];
                
                psi_y[idx] = f[idx_y_T] - 2*f[idx_y] + f[idx_y_B];
                mean_y[idx] = f[idx_y_T] + 2*f[idx_y] + f[idx_y_B];
            }
        }
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                // Compute the linear index.
                const int idx = (i + num_ghosts_gradient[0]) +
                    (j + num_ghosts_gradient[1])*ghostcell_dims_gradient[0];
                
                psi[idx] = sqrt(psi_x[idx]*psi_x[idx] + psi_y[idx]*psi_y[idx])/
                    (sqrt(mean_x[idx]*mean_x[idx] + mean_y[idx]*mean_y[idx]) + EPSILON);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        // Allocate memory in different dimensions.
        boost::shared_ptr<pdat::CellData<double> > gradient_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > gradient_y(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > gradient_z(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_x(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_y(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        boost::shared_ptr<pdat::CellData<double> > local_mean_value_z(
            new pdat::CellData<double>(interior_box, 1, num_ghosts_gradient));
        
        double* psi_x = gradient_x->getPointer(0);
        double* psi_y = gradient_y->getPointer(0);
        double* psi_z = gradient_z->getPointer(0);
        double* mean_x = local_mean_value_x->getPointer(0);
        double* mean_y = local_mean_value_y->getPointer(0);
        double* mean_z = local_mean_value_z->getPointer(0);
        
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_gradient[0]) +
                        (j + num_ghosts_gradient[1])*ghostcell_dims_gradient[0] +
                        (k + num_ghosts_gradient[2])*ghostcell_dims_gradient[0]*
                            ghostcell_dims_gradient[1];
                    
                    const int idx_x_L = (i - 1 + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    const int idx_x = (i + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    const int idx_x_R = (i + 1 + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    psi_x[idx] = f[idx_x_R] - 2*f[idx_x] + f[idx_x_L];
                    mean_x[idx] = f[idx_x_R] + 2*f[idx_x] + f[idx_x_L];
                }
            }
        }
        
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_gradient[0]) +
                        (j + num_ghosts_gradient[1])*ghostcell_dims_gradient[0] +
                        (k + num_ghosts_gradient[2])*ghostcell_dims_gradient[0]*
                            ghostcell_dims_gradient[1];
                    
                    const int idx_y_B = (i + num_ghosts_cell_data[0]) +
                        (j - 1 + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    const int idx_y = (i + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    const int idx_y_T = (i + num_ghosts_cell_data[0]) +
                        (j + 1 + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    psi_y[idx] = f[idx_y_T] - 2*f[idx_y] + f[idx_y_B];
                    mean_y[idx] = f[idx_y_T] + 2*f[idx_y] + f[idx_y_B];
                }
            }
        }
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                for (int k = 0; k < interior_dims[2]; k++)
                {
                    // Compute the linear indices.
                    const int idx = (i + num_ghosts_gradient[0]) +
                        (j + num_ghosts_gradient[1])*ghostcell_dims_gradient[0] +
                        (k + num_ghosts_gradient[2])*ghostcell_dims_gradient[0]*
                            ghostcell_dims_gradient[1];
                    
                    const int idx_z_B = (i + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k - 1 + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    const int idx_z = (i + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    const int idx_z_F = (i + num_ghosts_cell_data[0]) +
                        (j + num_ghosts_cell_data[1])*ghostcell_dims_cell_data[0] +
                        (k + 1 + num_ghosts_cell_data[2])*ghostcell_dims_cell_data[0]*
                            ghostcell_dims_cell_data[1];
                    
                    psi_z[idx] = f[idx_z_F] - 2*f[idx_z] + f[idx_z_B];
                    mean_z[idx] = f[idx_z_F] + 2*f[idx_z] + f[idx_z_B];
                }
            }
        }
        
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the index.
                    const int idx = (i + num_ghosts_gradient[0]) +
                        (j + num_ghosts_gradient[1])*ghostcell_dims_gradient[0] +
                        (k + num_ghosts_gradient[2])*ghostcell_dims_gradient[0]*
                            ghostcell_dims_gradient[1];
                    
                    psi[idx] = sqrt(psi_x[idx]*psi_x[idx] + psi_y[idx]*psi_y[idx] + psi_z[idx]*psi_z[idx])/
                        (sqrt(mean_x[idx]*mean_x[idx] + mean_y[idx]*mean_y[idx] + mean_z[idx]*mean_z[idx]) +
                         EPSILON);
                }
            }
        }
    }
}
