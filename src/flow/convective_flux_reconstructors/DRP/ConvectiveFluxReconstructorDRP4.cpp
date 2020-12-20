#include "flow/convective_flux_reconstructors/DRP/ConvectiveFluxReconstructorDRP4.hpp"

/*
 * Timers interspersed throughout the class.
 */

HAMERS_SHARED_PTR<tbox::Timer> ConvectiveFluxReconstructorDRP4::t_reconstruct_flux;
HAMERS_SHARED_PTR<tbox::Timer> ConvectiveFluxReconstructorDRP4::t_compute_source;


ConvectiveFluxReconstructorDRP4::ConvectiveFluxReconstructorDRP4(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructor(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            convective_flux_reconstructor_db)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
    
    d_stencil_width = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("stencil_width", 9);
    d_stencil_width = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("d_stencil_width", d_stencil_width);
    
    d_eqn_form = d_flow_model->getEquationsForm();
    d_has_advective_eqn_form = false;
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
        {
            d_has_advective_eqn_form = true;
        }
    }
    
    t_reconstruct_flux = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorDRP4::t_reconstruct_flux");
    
    t_compute_source = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorDRP4::t_compute_source");
}


ConvectiveFluxReconstructorDRP4::~ConvectiveFluxReconstructorDRP4()
{
    t_reconstruct_flux.reset();
    t_compute_source.reset();
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorDRP4::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorDRP4 object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorDRP4: this = "
       << (ConvectiveFluxReconstructorDRP4 *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_stencil_width = "
       << d_stencil_width
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorDRP4::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putDouble("d_stencil_width", d_stencil_width);
}


/*
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorDRP4::computeConvectiveFluxAndSourceOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::SideVariable<double> >& variable_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
    if (d_dim == tbox::Dimension(1))
    {
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
    } // if (d_dim == tbox::Dimension(3))
}
