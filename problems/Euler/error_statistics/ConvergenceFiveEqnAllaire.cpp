#include "apps/Euler/EulerErrorStatistics.hpp"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

void
EulerErrorStatistics::printErrorStatistics(
    std::ostream& os,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& variable_context,
    const double time) const
{
    if ((d_project_name != "2D convergence test five-eqn by Allaire") &&
        (d_project_name != "3D convergence test five-eqn by Allaire"))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D convergence test five-eqn by Allaire' or "
            << "'3D convergence test five-eqn by Allaire'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    math::HierarchyCellDataOpsReal<double> cell_double_operator(patch_hierarchy, 0, 0);
    
    std::vector<std::string> variable_names = d_flow_model->getNamesOfConservativeVariables();
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<double> > > variables =
        d_flow_model->getConservativeVariables();
    
    if (d_project_name == "2D convergence test five-eqn by Allaire")
    {
        for (int vi = 0; vi < static_cast<int>(variables.size()); vi++)
        {
            if (variable_names[vi] == "volume fractions")
            {
                HAMERS_SHARED_PTR<hier::PatchLevel> level_root(
                patch_hierarchy->getPatchLevel(0));
                
                double dvol_sum_local = double(0);
                
                double error_sum_local         = double(0);
                double error_squared_sum_local = double(0);
                double error_max_local         = double(0);
                
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
                    
                    HAMERS_SHARED_PTR<pdat::CellData<double> > volume_fraction(
                        HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch->getPatchData(variables[vi], variable_context)));
                    
                    double* Z_1 = volume_fraction->getPointer(0);
                    
                    for (int j = 0; j < patch_dims[1]; j++)
                    {
                        for (int i = 0; i < patch_dims[0]; i++)
                        {
                            // Compute index into linear data array.
                            int idx_cell = i + j*patch_dims[0];
                            
                            // Compute the coordinates.
                            double x[2];
                            x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                            x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                            
                            const double Z_1_exact =
                                double(1)/double(2) + double(1)/double(4)*sin(M_PI*(x[0] + x[1] - double(2)*time));
                            
                            const double error = fabs(Z_1_exact - Z_1[idx_cell]);
                            
                            dvol_sum_local += dx[0]*dx[1];
                            
                            error_sum_local         += dx[0]*dx[0]*error;
                            error_squared_sum_local += dx[0]*dx[0]*error*error;
                            error_max_local         = fmax(error, error_max_local);
                        }
                    }
                }
                
                double dvol_sum_global = double(0);
                
                double error_sum_global         = double(0);
                double error_squared_sum_global = double(0);
                double error_max_global         = double(0);
                
                mpi.Allreduce(
                    &dvol_sum_local,
                    &dvol_sum_global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM);
                
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
                
                const double L1_error   = error_sum_global/dvol_sum_global;
                const double L2_error   = sqrt(error_squared_sum_global/dvol_sum_global);
                const double Linf_error = error_max_global;
                
                os.precision(17);
                os << "L1_error: " << std::scientific << L1_error << std::endl;
                os << "L2_error: " << std::scientific << L2_error << std::endl;
                os << "Linf_error: " << std::scientific << Linf_error << std::endl;
            }
        }
    }
    else if (d_project_name == "3D convergence test five-eqn by Allaire")
    {
        for (int vi = 0; vi < static_cast<int>(variables.size()); vi++)
        {
            if (variable_names[vi] == "volume fractions")
            {
                HAMERS_SHARED_PTR<hier::PatchLevel> level_root(
                patch_hierarchy->getPatchLevel(0));
                
                double dvol_sum_local = double(0);
                
                double error_sum_local         = double(0);
                double error_squared_sum_local = double(0);
                double error_max_local         = double(0);
                
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
                    
                    HAMERS_SHARED_PTR<pdat::CellData<double> > volume_fraction(
                        HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch->getPatchData(variables[vi], variable_context)));
                    
                    double* Z_1 = volume_fraction->getPointer(0);
                    
                    for (int k = 0; k < patch_dims[2]; k++)
                    {
                        for (int j = 0; j < patch_dims[1]; j++)
                        {
                            for (int i = 0; i < patch_dims[0]; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i +
                                    j*patch_dims[0] +
                                    k*patch_dims[0]*patch_dims[1];
                                
                                // Compute the coordinates.
                                double x[3];
                                x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                                x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                                x[2] = patch_xlo[2] + (double(k) + double(1)/double(2))*dx[2];
                                
                                const double Z_1_exact =
                                    double(1)/double(2) + double(1)/double(4)*sin(M_PI*(x[0] + x[1] + x[2] - double(3)*time));
                                
                                const double error = fabs(Z_1_exact - Z_1[idx_cell]);
                                
                                dvol_sum_local += dx[0]*dx[1]*dx[2];
                                
                                error_sum_local         += dx[0]*dx[1]*dx[2]*error;
                                error_squared_sum_local += dx[0]*dx[1]*dx[2]*error*error;
                                error_max_local         = fmax(error, error_max_local);
                            }
                        }
                    }
                }
                
                double dvol_sum_global = double(0);
                
                double error_sum_global         = double(0);
                double error_squared_sum_global = double(0);
                double error_max_global         = double(0);
                
                mpi.Allreduce(
                    &dvol_sum_local,
                    &dvol_sum_global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM);
                
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
                
                const double L1_error   = error_sum_global/dvol_sum_global;
                const double L2_error   = sqrt(error_squared_sum_global/dvol_sum_global);
                const double Linf_error = error_max_global;
                
                os.precision(17);
                os << "L1_error: " << std::scientific << L1_error << std::endl;
                os << "L2_error: " << std::scientific << L2_error << std::endl;
                os << "Linf_error: " << std::scientific << Linf_error << std::endl;
            }
        }
    }
}
