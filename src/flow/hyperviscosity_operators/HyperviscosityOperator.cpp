#include "flow/hyperviscosity_operators/HyperviscosityOperator.hpp"


HyperviscosityOperator::HyperviscosityOperator(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& hyperviscosity_operator_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_hyperviscosity_op_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_eqn(num_eqn),
        d_flow_model(flow_model),
        d_hyperviscosity_operator_db(
            hyperviscosity_operator_db)
{
    d_lap_order = d_hyperviscosity_operator_db->
        getIntegerWithDefault("lap_order", 6);
    d_lap_order = d_hyperviscosity_operator_db->
        getIntegerWithDefault("d_lap_order", d_lap_order);
    
    d_accuracy_order = d_hyperviscosity_operator_db->
        getIntegerWithDefault("accuracy_order", 6);
    d_accuracy_order = d_hyperviscosity_operator_db->
        getIntegerWithDefault("d_accuracy_order", d_accuracy_order);
    
    d_use_conservative_form = d_hyperviscosity_operator_db->
        getBoolWithDefault("use_conservative_form", false);
    d_use_conservative_form = d_hyperviscosity_operator_db->
        getBoolWithDefault("d_use_conservative_form", d_use_conservative_form);
    
    if (d_accuracy_order != 2 && d_accuracy_order != 4  && d_accuracy_order != 6)
    {
        TBOX_ERROR("HyperviscosityOperator::HyperviscosityOperator:"
            " Only 2nd, 4th, 6th order accurate schemes are implemented!");
    }
    
    if (d_lap_order == 2)
    {
        if (d_accuracy_order == 2)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim);
        }
        else if (d_accuracy_order == 4)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*2;
        }
        else if (d_accuracy_order == 6)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*3;
        }
    }
    else if (d_lap_order == 4)
    {
        if (d_accuracy_order == 2)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*2;
        }
        else if (d_accuracy_order == 4)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*3;
        }
        else if (d_accuracy_order == 6)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*4;
        }
    }
    else if (d_lap_order == 6)
    {
        if (d_accuracy_order == 2)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*3;
        }
        else if (d_accuracy_order == 4)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*4;
        }
        else if (d_accuracy_order == 6)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*5;
        }
    }
    else
    {
        TBOX_ERROR("HyperviscosityOperator::HyperviscosityOperator:"
            " Only 2nd, 4th, 6th order Laplacian are implemented!");
    }
}


/*
 * Print all characteristics of the non-conservative diffusive flux divergence operator class.
 */
void
HyperviscosityOperator::printClassData(
    std::ostream& os) const
{
    os << "\nPrint HyperviscosityOperator object..."
       << std::endl;
    
    os << std::endl;
    
    os << "HyperviscosityOperator: this = "
       << (HyperviscosityOperator *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_lap_order = "
       << d_lap_order
       << std::endl;
    os << "d_accuracy_order = "
       << d_accuracy_order
       << std::endl;
    os << "d_use_conservative_form = "
       << d_use_conservative_form
       << std::endl;
}


/*
 * Put the characteristics of the hyperviscosity operator into the restart database.
 */
void
HyperviscosityOperator::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_lap_order", d_lap_order);
    restart_db->putInteger("d_accuracy_order", d_accuracy_order);
}


/*
 * Perform the hyperviscosity operator on a patch.
 */
void
HyperviscosityOperator::performHyperviscosityOperatorOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::CoarseFineBoundary> coarse_fine_bdry,
    const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(patch);
    NULL_USE(coarse_fine_bdry);
    NULL_USE(variable_convective_flux);
    NULL_USE(variable_source);
    NULL_USE(data_context);
    NULL_USE(time);
    NULL_USE(dt);
    NULL_USE(RK_step_number);
    
    Real a_n = Real(0);
    Real b_n = Real(0);
    Real c_n = Real(0);
    Real d_n = Real(0);
    Real e_n = Real(0);
    Real f_n = Real(0);
    
    if (d_lap_order == 2)
    {
        if (d_accuracy_order == 2)
        {
            a_n = -Real(2);
            b_n =  Real(1);
        }
        else if (d_accuracy_order == 4)
        {
            a_n = -Real(5)/Real(2);
            b_n =  Real(4)/Real(3);
            c_n = -Real(1)/Real(12);
        }
        else if (d_accuracy_order == 6)
        {
            a_n = -Real(49)/Real(18);
            b_n =  Real(3)/Real(2);
            c_n = -Real(3)/Real(20);
            d_n =  Real(1)/Real(90);
        }
    }
    else if (d_lap_order == 4)
    {
        if (d_accuracy_order == 2)
        {
            a_n =  Real(6);
            b_n = -Real(4);
            c_n =  Real(1);
        }
        else if (d_accuracy_order == 4)
        {
            a_n =  Real(28)/Real(3);
            b_n = -Real(13)/Real(2);
            c_n =  Real(2);
            d_n = -Real(1)/Real(6);
        }
        else if (d_accuracy_order == 6)
        {
            a_n =  Real(91)/Real(8);
            b_n = -Real(122)/Real(15);
            c_n =  Real(169)/Real(60);
            d_n = -Real(2)/Real(5);
            e_n =  Real(7)/Real(240);
        }
    }
    else if (d_lap_order == 6)
    {
        if (d_accuracy_order == 2)
        {
            a_n = -Real(20);
            b_n =  Real(15);
            c_n = -Real(6);
            d_n =  Real(1);
        }
        else if (d_accuracy_order == 4)
        {
            a_n = -Real(75)/Real(2);
            b_n =  Real(29);
            c_n = -Real(13);
            d_n =  Real(3);
            e_n = -Real(1)/Real(4);
        }
        else if (d_accuracy_order == 6)
        {
            a_n = -Real(1023)/Real(20);
            b_n =  Real(323)/Real(8);
            c_n = -Real(39)/Real(2);
            d_n =  Real(87)/Real(16);
            e_n = -Real(19)/Real(24);
            f_n =  Real(13)/Real(240);
        }
    }
    
    const Real a_m = b_n + c_n + d_n + e_n + f_n;
    const Real b_m = c_n + d_n + e_n + f_n;
    const Real c_m = d_n + e_n + f_n;
    const Real d_m = e_n + f_n;
    const Real e_m = f_n;
}
