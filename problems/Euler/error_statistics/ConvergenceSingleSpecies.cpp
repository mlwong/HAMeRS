#include "apps/Euler/EulerErrorStatistics.hpp"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

/*
 * Set the data on the patch interior to some initial values.
 */
void
EulerErrorStatistics::printErrorStatistics(
    std::ostream& os,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& variable_context) const
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, 0, 0);
    
    std::vector<std::string> variable_names = d_flow_model->getNamesOfConservativeVariables();
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<double> > > variables =
        d_flow_model->getConservativeVariables();
    
   if (d_project_name == "2D convergence test single-species")
   {
        for (int vi = 0; vi < static_cast<int>(variables.size()); vi++)
        {
            if (variable_names[vi] == "density")
            {
                HAMERS_SHARED_PTR<hier::PatchLevel> level_root(
                    patch_hierarchy->getPatchLevel(0));
                
                double error_sum_local = 0.0;
                double error_squared_sum_local = 0.0;
                double error_max_local = 0.0;
                
                for (hier::PatchLevel::iterator ip(level_root->begin());
                     ip != level_root->end();
                     ip++)
                {
                    const HAMERS_SHARED_PTR<hier::Patch>& patch = *ip;
                    
                    // Get the dimensions of box that covers the interior of Patch.
                    hier::Box patch_box = patch->getBox();
                    const hier::IntVector patch_dims = patch_box.numberCells();
                    
                    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch->getPatchGeometry()));
                    
                    const double* const dx = patch_geom->getDx();
                    const double* const patch_xlo = patch_geom->getXLower();
                    
                    HAMERS_SHARED_PTR<pdat::CellData<double> > density(
                        HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch->getPatchData(variables[vi], variable_context)));
                    
                    double* rho = density->getPointer(0);
                    
                    for (int j = 0; j < patch_dims[1]; j++)
                    {
                        for (int i = 0; i < patch_dims[0]; i++)
                        {
                            // Compute index into linear data array.
                            int idx_cell = i + j*patch_dims[0];
                            
                            // Compute the coordinates.
                            double x[2];
                            x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                            x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                            
                            const double rho_exact   = 1.0 + 0.5*sin(M_PI*(x[0] + x[1]));
                            const double error = fabs(rho_exact - rho[idx_cell]);
                            
                            error_sum_local += dx[0]*dx[0]*error;
                            error_squared_sum_local += dx[0]*dx[0]*error*error;
                            error_max_local = fmax(error, error_max_local);
                        }
                    }
                }
                
                double error_sum_global = 0.0;
                double error_squared_sum_global = 0.0;
                double error_max_global = 0.0;
                
                mpi.Allreduce(
                    &error_sum_local,
                    &error_sum_global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM);
                
                mpi.Allreduce(
                    &error_squared_sum_local,
                    &error_squared_sum_global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM);
                
                mpi.Allreduce(
                    &error_max_local,
                    &error_max_global,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX);
                
                
                os.precision(17);
                os << "L1_error: " << std::scientific << error_sum_global << std::endl;
                os << "L2_error: " << std::scientific << sqrt(error_squared_sum_global) << std::endl;
                os << "Linf_error: " << std::scientific << error_max_global << std::endl;
            }
        }
    }
}
