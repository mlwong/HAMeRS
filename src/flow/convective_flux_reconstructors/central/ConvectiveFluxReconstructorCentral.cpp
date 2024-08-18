#include "flow/convective_flux_reconstructors/central/ConvectiveFluxReconstructorCentral.hpp"

#define EPSILON HAMERS_EPSILON


/*
 * Timers interspersed throughout the class.
 */

HAMERS_SHARED_PTR<tbox::Timer> ConvectiveFluxReconstructorCentral::t_reconstruct_flux;
HAMERS_SHARED_PTR<tbox::Timer> ConvectiveFluxReconstructorCentral::t_compute_source;


/*
 * Interger based power function.
 */
static inline __attribute__((always_inline)) Real ipow(Real base, const int& exp)
{
    Real result = base;
    for (int i = 1; i < exp; i++)
    {
        result *= base;
    }

    return result;
}


/*
 * Compute local sigma.
 */
static inline __attribute__((always_inline)) void computeLocalSigma(
    Real& sigma,
    Real** U_array,
    const int& idx_side)
{
    /*
     * Compute the sigma.
     */
    
    const Real alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
    const Real alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
    const Real alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
    
    const Real theta_1 = std::abs(alpha_1 - alpha_2)/(std::abs(alpha_1) + std::abs(alpha_2) + EPSILON);
    const Real theta_2 = std::abs(alpha_2 - alpha_3)/(std::abs(alpha_2) + std::abs(alpha_3) + EPSILON);
    
    sigma = std::max(theta_1, theta_2);
}


/*
 * Compute local beta's.
 */
static inline __attribute__((always_inline)) void computeLocalBeta(
    Real& beta_0,
    Real& beta_1,
    Real& beta_2,
    Real& beta_3,
    Real** U_array,
    const int& idx_side)
{
    beta_0 = Real(1)/Real(3)*(U_array[0][idx_side]*(Real(4)*U_array[0][idx_side] -
         Real(19)*U_array[1][idx_side] + Real(11)*U_array[2][idx_side]) +
         U_array[1][idx_side]*(Real(25)*U_array[1][idx_side] - Real(31)*U_array[2][idx_side]) +
         Real(10)*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_1 = Real(1)/Real(3)*(U_array[1][idx_side]*(Real(4)*U_array[1][idx_side] -
         Real(13)*U_array[2][idx_side] + Real(5)*U_array[3][idx_side]) +
         Real(13)*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
         Real(4)*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_2 = Real(1)/Real(3)*(U_array[2][idx_side]*(Real(10)*U_array[2][idx_side] -
         Real(31)*U_array[3][idx_side] + Real(11)*U_array[4][idx_side]) +
         U_array[3][idx_side]*(Real(25)*U_array[3][idx_side] - Real(19)*U_array[4][idx_side]) +
         Real(4)*U_array[4][idx_side]*U_array[4][idx_side]);
    
    beta_3 = Real(1)/Real(232243200)*(U_array[0][idx_side]*(Real(525910327)*U_array[0][idx_side] -
         Real(4562164630)*U_array[1][idx_side] + Real(7799501420)*U_array[2][idx_side] -
         Real(6610694540)*U_array[3][idx_side] + Real(2794296070)*U_array[4][idx_side] -
         Real(472758974)*U_array[5][idx_side]) + Real(5)*U_array[1][idx_side]*
        (Real(2146987907)*U_array[1][idx_side] - Real(7722406988)*U_array[2][idx_side] +
         Real(6763559276)*U_array[3][idx_side] - Real(2926461814)*U_array[4][idx_side] +
         Real(503766638)*U_array[5][idx_side]) + Real(20)*U_array[2][idx_side]*
        (Real(1833221603)*U_array[2][idx_side] - Real(3358664662)*U_array[3][idx_side] +
         Real(1495974539)*U_array[4][idx_side] - Real(263126407)*U_array[5][idx_side]) +
        Real(20)*U_array[3][idx_side]*(Real(1607794163)*U_array[3][idx_side] -
         Real(1486026707)*U_array[4][idx_side] + Real(268747951)*U_array[5][idx_side]) +
        Real(5)*U_array[4][idx_side]*(Real(1432381427)*U_array[4][idx_side] -
         Real(536951582)*U_array[5][idx_side]) +
        Real(263126407)*U_array[5][idx_side]*U_array[5][idx_side]);
}


/*
 * Compute local beta_tilde's.
 */
static inline __attribute__((always_inline)) void computeLocalBetaTilde(
    Real& beta_tilde_0,
    Real& beta_tilde_1,
    Real& beta_tilde_2,
    Real& beta_tilde_3,
    Real** U_array,
    const int& idx_side)
{
    beta_tilde_0 = Real(1)/Real(3)*(U_array[5][idx_side]*(Real(4)*U_array[5][idx_side] -
         Real(19)*U_array[4][idx_side] + Real(11)*U_array[3][idx_side]) +
         U_array[4][idx_side]*(Real(25)*U_array[4][idx_side] - Real(31)*U_array[3][idx_side]) +
         Real(10)*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_tilde_1 = Real(1)/Real(3)*(U_array[4][idx_side]*(Real(4)*U_array[4][idx_side] -
         Real(13)*U_array[3][idx_side] + Real(5)*U_array[2][idx_side]) +
         Real(13)*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
         Real(4)*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_tilde_2 = Real(1)/Real(3)*(U_array[3][idx_side]*(Real(10)*U_array[3][idx_side] -
         Real(31)*U_array[2][idx_side] + Real(11)*U_array[1][idx_side]) +
         U_array[2][idx_side]*(Real(25)*U_array[2][idx_side] - Real(19)*U_array[1][idx_side]) +
         Real(4)*U_array[1][idx_side]*U_array[1][idx_side]);
    
    beta_tilde_3 = Real(1)/Real(232243200)*(U_array[5][idx_side]*(Real(525910327)*U_array[5][idx_side] -
         Real(4562164630)*U_array[4][idx_side] + Real(7799501420)*U_array[3][idx_side] -
         Real(6610694540)*U_array[2][idx_side] + Real(2794296070)*U_array[1][idx_side] -
         Real(472758974)*U_array[0][idx_side]) + Real(5)*U_array[4][idx_side]*
        (Real(2146987907)*U_array[4][idx_side] - Real(7722406988)*U_array[3][idx_side] +
         Real(6763559276)*U_array[2][idx_side] - Real(2926461814)*U_array[1][idx_side] +
         Real(503766638)*U_array[0][idx_side]) + Real(20)*U_array[3][idx_side]*
        (Real(1833221603)*U_array[3][idx_side] - Real(3358664662)*U_array[2][idx_side] +
         Real(1495974539)*U_array[1][idx_side] - Real(263126407)*U_array[0][idx_side]) +
        Real(20)*U_array[2][idx_side]*(Real(1607794163)*U_array[2][idx_side] -
         Real(1486026707)*U_array[1][idx_side] + Real(268747951)*U_array[0][idx_side]) +
        Real(5)*U_array[1][idx_side]*(Real(1432381427)*U_array[1][idx_side] -
         Real(536951582)*U_array[0][idx_side]) +
        Real(263126407)*U_array[0][idx_side]*U_array[0][idx_side]);
}


/*
 * Perform local WENO interpolation of U_minus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationMinus(
    Real* U_minus,
    Real** U_array,
    const int& idx_side,
    const int& p,
    const int& q,
    const Real& C,
    const Real& alpha_tau)
{
    /*
     * Compute sigma.
     */
    
    Real sigma;
    
    computeLocalSigma(sigma, U_array, idx_side);
    
    /*
     * Compute beta's.
     */
    
    Real beta_0, beta_1, beta_2, beta_3;
    
    computeLocalBeta(beta_0, beta_1, beta_2, beta_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind.
     */
    
    Real omega_upwind_0, omega_upwind_1, omega_upwind_2;
    
    Real tau_5 = std::abs(beta_0 - beta_2);
    
    omega_upwind_0 = Real(1)/Real(16)*(Real(1) + ipow(tau_5/(beta_0 + EPSILON), p));
    omega_upwind_1 = Real(5)/Real(8)*(Real(1) + ipow(tau_5/(beta_1 + EPSILON), p));
    omega_upwind_2 = Real(5)/Real(16)*(Real(1) + ipow(tau_5/(beta_2 + EPSILON), p));
    
    Real omega_upwind_sum = omega_upwind_0 + omega_upwind_1 + omega_upwind_2;
    
    omega_upwind_0 = omega_upwind_0/omega_upwind_sum;
    omega_upwind_1 = omega_upwind_1/omega_upwind_sum;
    omega_upwind_2 = omega_upwind_2/omega_upwind_sum;
    
    /*
     * Compute the weights omega_central (store in omega first).
     */
    
    Real omega_0, omega_1, omega_2, omega_3;
    
    Real beta_avg = Real(1)/Real(8)*(beta_0 + beta_2 + Real(6)*beta_1);
    Real tau_6 = std::abs(beta_3 - beta_avg);
    
    omega_0 = Real(1)/Real(32)*(C + ipow(tau_6/(beta_0 + EPSILON), q));
    omega_1 = Real(15)/Real(32)*(C + ipow(tau_6/(beta_1 + EPSILON), q));
    omega_2 = Real(15)/Real(32)*(C + ipow(tau_6/(beta_2 + EPSILON), q));
    omega_3 = Real(1)/Real(32)*(C + ipow(tau_6/(beta_3 + EPSILON), q));
    
    Real omega_sum = omega_0 + omega_1 + omega_2 + omega_3;
    
    omega_0 = omega_0/omega_sum;
    omega_1 = omega_1/omega_sum;
    omega_2 = omega_2/omega_sum;
    omega_3 = omega_3/omega_sum;
    
    /*
     * Compute the weights omega.
     */
    
    Real R_tau = tau_6/(beta_avg + EPSILON);
    
    if (R_tau > alpha_tau)
    {
        omega_0 = sigma*omega_upwind_0 + (Real(1) - sigma)*omega_0;
        omega_1 = sigma*omega_upwind_1 + (Real(1) - sigma)*omega_1;
        omega_2 = sigma*omega_upwind_2 + (Real(1) - sigma)*omega_2;
        omega_3 = (Real(1) - sigma)*omega_3;
    }
    
    /*
     * Compute U_minus.
     */
    
    U_minus[idx_side] = Real(3)/Real(8)*omega_0*U_array[0][idx_side] +
        (-Real(10)/Real(8)*omega_0 - Real(1)/Real(8)*omega_1)*U_array[1][idx_side] +
        (Real(15)/Real(8)*omega_0 + Real(6)/Real(8)*omega_1 + Real(3)/Real(8)*omega_2)*
            U_array[2][idx_side] +
        (Real(3)/Real(8)*omega_1 + Real(6)/Real(8)*omega_2 + Real(15)/Real(8)*omega_3)*
            U_array[3][idx_side] +
        (-Real(1)/Real(8)*omega_2 - Real(10)/Real(8)*omega_3)*U_array[4][idx_side] +
        Real(3)/Real(8)*omega_3*U_array[5][idx_side];
}


/*
 * Perform local WENO interpolation of U_plus.
 */
static inline __attribute__((always_inline)) void performLocalWENOInterpolationPlus(
    Real* U_plus,
    Real** U_array,
    const int& idx_side,
    const int& p,
    const int& q,
    const Real& C,
    const Real& alpha_tau)
{
    /*
     * Compute sigma.
     */
    
    Real sigma;
    
    computeLocalSigma(sigma, U_array, idx_side);
    
    /*
     * Compute beta_tilde's.
     */
    
    Real beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3;
    
    computeLocalBetaTilde(beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind_tilde.
     */
    
    Real omega_upwind_tilde_0, omega_upwind_tilde_1, omega_upwind_tilde_2;
    
    Real tau_5_tilde = std::abs(beta_tilde_0 - beta_tilde_2);
    
    omega_upwind_tilde_0 = Real(1)/Real(16)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_0 + EPSILON), p));
    omega_upwind_tilde_1 = Real(5)/Real(8)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_1 + EPSILON), p));
    omega_upwind_tilde_2 = Real(5)/Real(16)*(Real(1) + ipow(tau_5_tilde/(beta_tilde_2 + EPSILON), p));
    
    Real omega_upwind_tilde_sum = omega_upwind_tilde_0 + omega_upwind_tilde_1 + omega_upwind_tilde_2;
    
    omega_upwind_tilde_0 = omega_upwind_tilde_0/omega_upwind_tilde_sum;
    omega_upwind_tilde_1 = omega_upwind_tilde_1/omega_upwind_tilde_sum;
    omega_upwind_tilde_2 = omega_upwind_tilde_2/omega_upwind_tilde_sum;
    
    /*
     * Compute the weights omega_central_tilde (store in omega_tilde first).
     */
    
    Real omega_tilde_0, omega_tilde_1, omega_tilde_2, omega_tilde_3;
    
    Real beta_avg_tilde = Real(1)/Real(8)*(beta_tilde_0 + beta_tilde_2 + Real(6)*beta_tilde_1);
    Real tau_6_tilde = std::abs(beta_tilde_3 - beta_avg_tilde);
    
    omega_tilde_0 = Real(1)/Real(32)*(C + ipow(tau_6_tilde/(beta_tilde_0 + EPSILON), q));
    omega_tilde_1 = Real(15)/Real(32)*(C + ipow(tau_6_tilde/(beta_tilde_1 + EPSILON), q));
    omega_tilde_2 = Real(15)/Real(32)*(C + ipow(tau_6_tilde/(beta_tilde_2 + EPSILON), q));
    omega_tilde_3 = Real(1)/Real(32)*(C + ipow(tau_6_tilde/(beta_tilde_3 + EPSILON), q));
    
    Real omega_tilde_sum = omega_tilde_0 + omega_tilde_1 + omega_tilde_2 + omega_tilde_3;
    
    omega_tilde_0 = omega_tilde_0/omega_tilde_sum;
    omega_tilde_1 = omega_tilde_1/omega_tilde_sum;
    omega_tilde_2 = omega_tilde_2/omega_tilde_sum;
    omega_tilde_3 = omega_tilde_3/omega_tilde_sum;
    
    /*
     * Compute the weights omega_tilde.
     */
    
    Real R_tau_tilde = tau_6_tilde/(beta_avg_tilde + EPSILON);
    
    if (R_tau_tilde > alpha_tau)
    {
        omega_tilde_0 = sigma*omega_upwind_tilde_0 + (Real(1) - sigma)*omega_tilde_0;
        omega_tilde_1 = sigma*omega_upwind_tilde_1 + (Real(1) - sigma)*omega_tilde_1;
        omega_tilde_2 = sigma*omega_upwind_tilde_2 + (Real(1) - sigma)*omega_tilde_2;
        omega_tilde_3 = (Real(1) - sigma)*omega_tilde_3;
    }
    
    /*
     * Compute U_plus.
     */
    
    U_plus[idx_side] = Real(3)/Real(8)*omega_tilde_0*U_array[5][idx_side] +
        (-Real(10)/Real(8)*omega_tilde_0 - Real(1)/Real(8)*omega_tilde_1)*U_array[4][idx_side] +
        (Real(15)/Real(8)*omega_tilde_0 + Real(6)/Real(8)*omega_tilde_1 + Real(3)/Real(8)*omega_tilde_2)*
            U_array[3][idx_side] +
        (Real(3)/Real(8)*omega_tilde_1 + Real(6)/Real(8)*omega_tilde_2 + Real(15)/Real(8)*omega_tilde_3)*
            U_array[2][idx_side] +
        (-Real(1)/Real(8)*omega_tilde_2 - Real(10)/Real(8)*omega_tilde_3)*U_array[1][idx_side] +
        Real(3)/Real(8)*omega_tilde_3*U_array[0][idx_side];
}


ConvectiveFluxReconstructorCentral::ConvectiveFluxReconstructorCentral(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const FLOW_MODEL::TYPE& flow_model_type,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructor(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model_type,
            flow_model,
            convective_flux_reconstructor_db)
{
    d_order = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("order", 12);
    d_order = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("d_order", d_order);
    
    
    if (d_order == 2)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim);
    }
    else if (d_order == 4)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*2;
    }
    else if (d_order == 6)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*3;
    }
    else if (d_order == 8)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
    }
    else if (d_order == 10)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*5;
    }
    else if (d_order == 12)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*6;
    }
    else
    {
        TBOX_ERROR("ConvectiveFluxReconstructorCentral::ConvectiveFluxReconstructorCentral:"
            " Only 2nd, 4th, 6th, 8th, 10th, 12th order central schemes are implemented!");
    }
    
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
        getTimer("ConvectiveFluxReconstructorCentral::t_reconstruct_flux");
    
    t_compute_source = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorCentral::t_compute_source");
}


ConvectiveFluxReconstructorCentral::~ConvectiveFluxReconstructorCentral()
{
    t_reconstruct_flux.reset();
    t_compute_source.reset();
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorCentral::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorCentral object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorCentral: this = "
       << (ConvectiveFluxReconstructorCentral *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_order = "
       << d_order
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorCentral::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_order", d_order);
}


/*
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorCentral::computeConvectiveFluxAndSourceOnPatch(
    hier::Patch& patch,
    const int level_number,
    const HAMERS_SHARED_PTR<hier::CoarseFineBoundary> coarse_fine_bdry,
    const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(level_number);
    NULL_USE(coarse_fine_bdry);
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
    Real a_n = Real(0);
    Real b_n = Real(0);
    Real c_n = Real(0);
    Real d_n = Real(0);
    Real e_n = Real(0);
    Real f_n = Real(0);
    
    Real a_m = Real(0);
    Real b_m = Real(0);
    Real c_m = Real(0);
    Real d_m = Real(0);
    Real e_m = Real(0);
    Real f_m = Real(0);
    
    if (d_order == 2)
    {
        a_n =  Real(1)/Real(2);
        
        a_m = a_n;
    }
    else if (d_order == 4)
    {
        a_n =  Real(2)/Real(3);
        b_n = -Real(1)/Real(12);
        
        a_m = a_n + b_n;
        b_m = b_n;
    }
    else if (d_order == 6)
    {
        a_n =  Real(3)/Real(4);
        b_n = -Real(3)/Real(20);
        c_n =  Real(1)/Real(60);
        
        a_m = a_n + b_n + c_n;
        b_m = b_n + c_n;
        c_m = c_n;
    }
    else if (d_order == 8)
    {
        a_n =  Real(4)/Real(5);
        b_n = -Real(1)/Real(5);
        c_n =  Real(4)/Real(105);
        d_n = -Real(1)/Real(280);
        
        a_m = a_n + b_n + c_n + d_n;
        b_m = b_n + c_n + d_n;
        c_m = c_n + d_n;
        d_m = d_n;
    }
    else if (d_order == 10)
    {
        a_n =  Real(5)/Real(6);
        b_n = -Real(5)/Real(21);
        c_n =  Real(5)/Real(84);
        d_n = -Real(5)/Real(504);
        e_n =  Real(1)/Real(1260);
        
        a_m = a_n + b_n + c_n + d_n + e_n;
        b_m = b_n + c_n + d_n + e_n;
        c_m = c_n + d_n + e_n;
        d_m = d_n + e_n;
        e_m = e_n;
    }
    else if (d_order == 12)
    {
        a_n =  Real(6)/Real(7);
        b_n = -Real(15)/Real(56);
        c_n =  Real(5)/Real(63);
        d_n = -Real(1)/Real(56);
        e_n =  Real(3)/Real(1155);
        f_n = -Real(1)/Real(5544);
        
        a_m = a_n + b_n + c_n + d_n + e_n + f_n;
        b_m = b_n + c_n + d_n + e_n + f_n;
        c_m = c_n + d_n + e_n + f_n;
        d_m = d_n + e_n + f_n;
        e_m = e_n + f_n;
        f_m = f_n;
    }
    else
    {
        TBOX_ERROR("ConvectiveFluxReconstructorCentral::computeConvectiveFluxAndSourceOnPatch:"
            " Only 8th, 10th, 12th order central schemes are implemented!");
    }
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // convective ghost cells.
    hier::Box conv_ghost_box = interior_box;
    conv_ghost_box.grow(d_num_conv_ghosts);
    const hier::IntVector conv_ghostcell_dims = conv_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    // Get the side data of convective flux.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux(
        HAMERS_SHARED_PTR_CAST<pdat::SideData<Real>, hier::PatchData>(
            patch.getPatchData(variable_convective_flux, data_context)));
    
    // Get the cell data of source.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > source(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            patch.getPatchData(variable_source, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    
    TBOX_ASSERT(source);
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        
        if (d_has_advective_eqn_form)
        {
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        }
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(2);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        
        std::vector<Real*> F_node_x;
        F_node_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
        }
        
        t_reconstruct_flux->start();
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        if (d_order == 2)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    
                    const int idx_node_L = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R = i     + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = Real(dt)*(
                        a_m*(F_node_x[ei][idx_node_L] + F_node_x[ei][idx_node_R])
                        );
                }
            }
        }
        else if (d_order == 4)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    
                    const int idx_node_LL = i - 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_L  = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R  = i     + num_subghosts_0_convective_flux_x;
                    const int idx_node_RR = i + 1 + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = Real(dt)*(
                        a_m*(F_node_x[ei][idx_node_L]  + F_node_x[ei][idx_node_R]) +
                        b_m*(F_node_x[ei][idx_node_LL] + F_node_x[ei][idx_node_RR])
                        );
                }
            }
        }
        else if (d_order == 6)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    
                    const int idx_node_LLL = i - 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LL  = i - 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_L   = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R   = i     + num_subghosts_0_convective_flux_x;
                    const int idx_node_RR  = i + 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRR = i + 2 + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = Real(dt)*(
                        a_m*(F_node_x[ei][idx_node_L]   + F_node_x[ei][idx_node_R]) +
                        b_m*(F_node_x[ei][idx_node_LL]  + F_node_x[ei][idx_node_RR]) +
                        c_m*(F_node_x[ei][idx_node_LLL] + F_node_x[ei][idx_node_RRR])
                        );
                }
            }
        }
        else if (d_order == 8)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    
                    const int idx_node_LLLL = i - 4 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLL  = i - 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LL   = i - 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_L    = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R    = i     + num_subghosts_0_convective_flux_x;
                    const int idx_node_RR   = i + 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRR  = i + 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRR = i + 3 + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = Real(dt)*(
                        a_m*(F_node_x[ei][idx_node_L]    + F_node_x[ei][idx_node_R]) +
                        b_m*(F_node_x[ei][idx_node_LL]   + F_node_x[ei][idx_node_RR]) +
                        c_m*(F_node_x[ei][idx_node_LLL]  + F_node_x[ei][idx_node_RRR]) +
                        d_m*(F_node_x[ei][idx_node_LLLL] + F_node_x[ei][idx_node_RRRR])
                        );
                }
            }
        }
        else if (d_order == 10)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    
                    const int idx_node_LLLLL = i - 5 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLLL  = i - 4 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLL   = i - 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LL    = i - 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_L     = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R     = i     + num_subghosts_0_convective_flux_x;
                    const int idx_node_RR    = i + 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRR   = i + 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRR  = i + 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRRR = i + 4 + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = Real(dt)*(
                        a_m*(F_node_x[ei][idx_node_L]     + F_node_x[ei][idx_node_R]) +
                        b_m*(F_node_x[ei][idx_node_LL]    + F_node_x[ei][idx_node_RR]) +
                        c_m*(F_node_x[ei][idx_node_LLL]   + F_node_x[ei][idx_node_RRR]) +
                        d_m*(F_node_x[ei][idx_node_LLLL]  + F_node_x[ei][idx_node_RRRR]) +
                        e_m*(F_node_x[ei][idx_node_LLLLL] + F_node_x[ei][idx_node_RRRRR])
                        );
                }
            }
        }
        else if (d_order == 12)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    
                    const int idx_node_LLLLLL = i - 6 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLLLL  = i - 5 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLLL   = i - 4 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLL    = i - 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LL     = i - 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_L      = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R      = i     + num_subghosts_0_convective_flux_x;
                    const int idx_node_RR     = i + 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRR    = i + 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRR   = i + 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRRR  = i + 4 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRRRR = i + 5 + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = Real(dt)*(
                        a_m*(F_node_x[ei][idx_node_L]      + F_node_x[ei][idx_node_R]) +
                        b_m*(F_node_x[ei][idx_node_LL]     + F_node_x[ei][idx_node_RR]) +
                        c_m*(F_node_x[ei][idx_node_LLL]    + F_node_x[ei][idx_node_RRR]) +
                        d_m*(F_node_x[ei][idx_node_LLLL]   + F_node_x[ei][idx_node_RRRR]) +
                        e_m*(F_node_x[ei][idx_node_LLLLL]  + F_node_x[ei][idx_node_RRRRR]) +
                        f_m*(F_node_x[ei][idx_node_LLLLLL] + F_node_x[ei][idx_node_RRRRRR])
                        );
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Compute the source.
         */
        
        t_compute_source->start();
        
        if (d_has_advective_eqn_form)
        {
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            
            Real* u = velocity->getPointer(0);
            
            std::vector<hier::IntVector> num_subghosts_conservative_var;
            num_subghosts_conservative_var.reserve(d_num_eqn);
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
                d_flow_model->getCellDataOfConservativeVariables();
            
            std::vector<Real*> Q;
            Q.reserve(d_num_eqn);
            
            int count_eqn = 0;
            
            for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
            {
                int depth = conservative_variables[vi]->getDepth();
                
                for (int di = 0; di < depth; di++)
                {
                    // If the last element of the conservative variable vector is not in the system of equations,
                    // ignore it.
                    if (count_eqn >= d_num_eqn)
                        break;
                    
                    Q.push_back(conservative_variables[vi]->getPointer(di));
                    num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                    
                    count_eqn++;
                }
            }
            
            for (int ei = 0; ei < d_num_eqn; ei ++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                    
                    if (d_order == 2)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                            
                            const int idx_cell_wghost_x_L = i - 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_R = i + 1 + num_subghosts_0_velocity;
                            
                            const int idx_cell_nghost = i;
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                (
                                a_n*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L])
                                )/Real(dx[0]));
                        }
                    }
                    else if (d_order == 4)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                            
                            const int idx_cell_wghost_x_LL = i - 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_L  = i - 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_R  = i + 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RR = i + 2 + num_subghosts_0_velocity;
                            
                            const int idx_cell_nghost = i;
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                (
                                a_n*(u[idx_cell_wghost_x_R]  - u[idx_cell_wghost_x_L]) +
                                b_n*(u[idx_cell_wghost_x_RR] - u[idx_cell_wghost_x_LL])
                                )/Real(dx[0]));
                        }
                    }
                    else if (d_order == 6)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                            
                            const int idx_cell_wghost_x_LLL = i - 3 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LL  = i - 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_L   = i - 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_R   = i + 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RR  = i + 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRR = i + 3 + num_subghosts_0_velocity;
                            
                            const int idx_cell_nghost = i;
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                (
                                a_n*(u[idx_cell_wghost_x_R]   - u[idx_cell_wghost_x_L]) +
                                b_n*(u[idx_cell_wghost_x_RR]  - u[idx_cell_wghost_x_LL]) +
                                c_n*(u[idx_cell_wghost_x_RRR] - u[idx_cell_wghost_x_LLL])
                                )/Real(dx[0]));
                        }
                    }
                    else if (d_order == 8)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                            
                            const int idx_cell_wghost_x_LLLL = i - 4 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LLL  = i - 3 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LL   = i - 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_L    = i - 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_R    = i + 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RR   = i + 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRR  = i + 3 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRRR = i + 4 + num_subghosts_0_velocity;
                            
                            const int idx_cell_nghost = i;
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                (
                                a_n*(u[idx_cell_wghost_x_R]    - u[idx_cell_wghost_x_L]) +
                                b_n*(u[idx_cell_wghost_x_RR]   - u[idx_cell_wghost_x_LL]) +
                                c_n*(u[idx_cell_wghost_x_RRR]  - u[idx_cell_wghost_x_LLL]) +
                                d_n*(u[idx_cell_wghost_x_RRRR] - u[idx_cell_wghost_x_LLLL])
                                )/Real(dx[0]));
                        }
                    }
                    else if (d_order == 10)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                            
                            const int idx_cell_wghost_x_LLLLL = i - 5 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LLLL  = i - 4 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LLL   = i - 3 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LL    = i - 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_L     = i - 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_R     = i + 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RR    = i + 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRR   = i + 3 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRRR  = i + 4 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRRRR = i + 5 + num_subghosts_0_velocity;
                            
                            const int idx_cell_nghost = i;
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                (
                                a_n*(u[idx_cell_wghost_x_R]     - u[idx_cell_wghost_x_L]) +
                                b_n*(u[idx_cell_wghost_x_RR]    - u[idx_cell_wghost_x_LL]) +
                                c_n*(u[idx_cell_wghost_x_RRR]   - u[idx_cell_wghost_x_LLL]) +
                                d_n*(u[idx_cell_wghost_x_RRRR]  - u[idx_cell_wghost_x_LLLL]) +
                                e_n*(u[idx_cell_wghost_x_RRRRR] - u[idx_cell_wghost_x_LLLLL])
                                )/Real(dx[0]));
                        }
                    }
                    else if (d_order == 12)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                            
                            const int idx_cell_wghost_x_LLLLLL = i - 6 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LLLLL  = i - 5 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LLLL   = i - 4 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LLL    = i - 3 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_LL     = i - 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_L      = i - 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_R      = i + 1 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RR     = i + 2 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRR    = i + 3 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRRR   = i + 4 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRRRR  = i + 5 + num_subghosts_0_velocity;
                            const int idx_cell_wghost_x_RRRRRR = i + 6 + num_subghosts_0_velocity;
                            
                            const int idx_cell_nghost = i;
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                (
                                a_n*(u[idx_cell_wghost_x_R]      - u[idx_cell_wghost_x_L]) +
                                b_n*(u[idx_cell_wghost_x_RR]     - u[idx_cell_wghost_x_LL]) +
                                c_n*(u[idx_cell_wghost_x_RRR]    - u[idx_cell_wghost_x_LLL]) +
                                d_n*(u[idx_cell_wghost_x_RRRR]   - u[idx_cell_wghost_x_LLLL]) +
                                e_n*(u[idx_cell_wghost_x_RRRRR]  - u[idx_cell_wghost_x_LLLLL]) +
                                f_n*(u[idx_cell_wghost_x_RRRRRR] - u[idx_cell_wghost_x_LLLLLL])
                                )/Real(dx[0]));
                        }
                    }
                }
            }
        }
        
        t_compute_source->stop();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        
        if (d_has_advective_eqn_form)
        {
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        }
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(2);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
        const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
        
        const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
        const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
        const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
        
        std::vector<Real*> F_node_x;
        std::vector<Real*> F_node_y;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
        }
        
        t_reconstruct_flux->start();
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        if (d_order == 2)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            a_m*(F_node_x[ei][idx_node_L] + F_node_x[ei][idx_node_R])
                            );
                    }
                }
            }
        }
        else if (d_order == 4)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        const int idx_node_LL = (i - 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_L  = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R  = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RR = (i + 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            a_m*(F_node_x[ei][idx_node_L]  + F_node_x[ei][idx_node_R]) +
                            b_m*(F_node_x[ei][idx_node_LL] + F_node_x[ei][idx_node_RR])
                            );
                    }
                }
            }
        }
        else if (d_order == 6)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        const int idx_node_LLL = (i - 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LL  = (i - 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_L   = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R   = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RR  = (i + 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRR = (i + 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            a_m*(F_node_x[ei][idx_node_L]   + F_node_x[ei][idx_node_R]) +
                            b_m*(F_node_x[ei][idx_node_LL]  + F_node_x[ei][idx_node_RR]) +
                            c_m*(F_node_x[ei][idx_node_LLL] + F_node_x[ei][idx_node_RRR])
                            );
                    }
                }
            }
        }
        else if (d_order == 8)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        const int idx_node_LLLL = (i - 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLL  = (i - 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LL   = (i - 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_L    = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R    = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RR   = (i + 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRR  = (i + 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRR = (i + 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            a_m*(F_node_x[ei][idx_node_L]    + F_node_x[ei][idx_node_R]) +
                            b_m*(F_node_x[ei][idx_node_LL]   + F_node_x[ei][idx_node_RR]) +
                            c_m*(F_node_x[ei][idx_node_LLL]  + F_node_x[ei][idx_node_RRR]) +
                            d_m*(F_node_x[ei][idx_node_LLLL] + F_node_x[ei][idx_node_RRRR])
                            );
                    }
                }
            }
        }
        else if (d_order == 10)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        const int idx_node_LLLLL = (i - 5 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLLL  = (i - 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLL   = (i - 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LL    = (i - 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_L     = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R     = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RR    = (i + 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRR   = (i + 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRR  = (i + 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRRR = (i + 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            a_m*(F_node_x[ei][idx_node_L]     + F_node_x[ei][idx_node_R]) +
                            b_m*(F_node_x[ei][idx_node_LL]    + F_node_x[ei][idx_node_RR]) +
                            c_m*(F_node_x[ei][idx_node_LLL]   + F_node_x[ei][idx_node_RRR]) +
                            d_m*(F_node_x[ei][idx_node_LLLL]  + F_node_x[ei][idx_node_RRRR]) +
                            e_m*(F_node_x[ei][idx_node_LLLLL] + F_node_x[ei][idx_node_RRRRR])
                            );
                    }
                }
            }
        }
        else if (d_order == 12)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        const int idx_node_LLLLLL = (i - 6 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLLLL  = (i - 5 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLLL   = (i - 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLL    = (i - 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LL     = (i - 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_L      = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R      = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RR     = (i + 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRR    = (i + 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRR   = (i + 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRRR  = (i + 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRRRR = (i + 5 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            a_m*(F_node_x[ei][idx_node_L]      + F_node_x[ei][idx_node_R]) +
                            b_m*(F_node_x[ei][idx_node_LL]     + F_node_x[ei][idx_node_RR]) +
                            c_m*(F_node_x[ei][idx_node_LLL]    + F_node_x[ei][idx_node_RRR]) +
                            d_m*(F_node_x[ei][idx_node_LLLL]   + F_node_x[ei][idx_node_RRRR]) +
                            e_m*(F_node_x[ei][idx_node_LLLLL]  + F_node_x[ei][idx_node_RRRRR]) +
                            f_m*(F_node_x[ei][idx_node_LLLLLL] + F_node_x[ei][idx_node_RRRRRR])
                            );
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        if (d_order == 2)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = Real(dt)*(
                            a_m*(F_node_y[ei][idx_node_B] + F_node_y[ei][idx_node_T])
                            );
                    }
                }
            }
        }
        else if (d_order == 4)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_BB = (i + num_subghosts_0_convective_flux_y) +
                            (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_B  = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T  = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TT = (i + num_subghosts_0_convective_flux_y) +
                            (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = Real(dt)*(
                            a_m*(F_node_y[ei][idx_node_B]  + F_node_y[ei][idx_node_T]) +
                            b_m*(F_node_y[ei][idx_node_BB] + F_node_y[ei][idx_node_TT])
                            );
                    }
                }
            }
        }
        else if (d_order == 6)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_BBB = (i + num_subghosts_0_convective_flux_y) +
                            (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BB  = (i + num_subghosts_0_convective_flux_y) +
                            (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_B   = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T   = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TT  = (i + num_subghosts_0_convective_flux_y) +
                            (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTT = (i + num_subghosts_0_convective_flux_y) +
                            (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = Real(dt)*(
                            a_m*(F_node_y[ei][idx_node_B]   + F_node_y[ei][idx_node_T]) +
                            b_m*(F_node_y[ei][idx_node_BB]  + F_node_y[ei][idx_node_TT]) +
                            c_m*(F_node_y[ei][idx_node_BBB] + F_node_y[ei][idx_node_TTT])
                            );
                    }
                }
            }
        }
        else if (d_order == 8)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_BBBB = (i + num_subghosts_0_convective_flux_y) +
                            (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBB  = (i + num_subghosts_0_convective_flux_y) +
                            (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BB   = (i + num_subghosts_0_convective_flux_y) +
                            (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_B    = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T    = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TT   = (i + num_subghosts_0_convective_flux_y) +
                            (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTT  = (i + num_subghosts_0_convective_flux_y) +
                            (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTT = (i + num_subghosts_0_convective_flux_y) +
                            (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = Real(dt)*(
                            a_m*(F_node_y[ei][idx_node_B]    + F_node_y[ei][idx_node_T]) +
                            b_m*(F_node_y[ei][idx_node_BB]   + F_node_y[ei][idx_node_TT]) +
                            c_m*(F_node_y[ei][idx_node_BBB]  + F_node_y[ei][idx_node_TTT]) +
                            d_m*(F_node_y[ei][idx_node_BBBB] + F_node_y[ei][idx_node_TTTT])
                            );
                    }
                }
            }
        }
        else if (d_order == 10)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_BBBBB = (i + num_subghosts_0_convective_flux_y) +
                            (j - 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBBB  = (i + num_subghosts_0_convective_flux_y) +
                            (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBB   = (i + num_subghosts_0_convective_flux_y) +
                            (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BB    = (i + num_subghosts_0_convective_flux_y) +
                            (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_B     = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T     = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TT    = (i + num_subghosts_0_convective_flux_y) +
                            (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTT   = (i + num_subghosts_0_convective_flux_y) +
                            (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTT  = (i + num_subghosts_0_convective_flux_y) +
                            (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTTT = (i + num_subghosts_0_convective_flux_y) +
                            (j + 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = Real(dt)*(
                            a_m*(F_node_y[ei][idx_node_B]     + F_node_y[ei][idx_node_T]) +
                            b_m*(F_node_y[ei][idx_node_BB]    + F_node_y[ei][idx_node_TT]) +
                            c_m*(F_node_y[ei][idx_node_BBB]   + F_node_y[ei][idx_node_TTT]) +
                            d_m*(F_node_y[ei][idx_node_BBBB]  + F_node_y[ei][idx_node_TTTT]) +
                            e_m*(F_node_y[ei][idx_node_BBBBB] + F_node_y[ei][idx_node_TTTTT])
                            );
                    }
                }
            }
        }
        else if (d_order == 12)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_BBBBBB = (i + num_subghosts_0_convective_flux_y) +
                            (j - 6 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBBBB  = (i + num_subghosts_0_convective_flux_y) +
                            (j - 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBBB   = (i + num_subghosts_0_convective_flux_y) +
                            (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBB    = (i + num_subghosts_0_convective_flux_y) +
                            (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BB     = (i + num_subghosts_0_convective_flux_y) +
                            (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_B      = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T      = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TT     = (i + num_subghosts_0_convective_flux_y) +
                            (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTT    = (i + num_subghosts_0_convective_flux_y) +
                            (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTT   = (i + num_subghosts_0_convective_flux_y) +
                            (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTTT  = (i + num_subghosts_0_convective_flux_y) +
                            (j + 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTTTT = (i + num_subghosts_0_convective_flux_y) +
                            (j + 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = Real(dt)*(
                            a_m*(F_node_y[ei][idx_node_B]      + F_node_y[ei][idx_node_T]) +
                            b_m*(F_node_y[ei][idx_node_BB]     + F_node_y[ei][idx_node_TT]) +
                            c_m*(F_node_y[ei][idx_node_BBB]    + F_node_y[ei][idx_node_TTT]) +
                            d_m*(F_node_y[ei][idx_node_BBBB]   + F_node_y[ei][idx_node_TTTT]) +
                            e_m*(F_node_y[ei][idx_node_BBBBB]  + F_node_y[ei][idx_node_TTTTT]) +
                            f_m*(F_node_y[ei][idx_node_BBBBBB] + F_node_y[ei][idx_node_TTTTTT])
                            );
                    }
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Compute the source.
         */
        
        t_compute_source->start();
        
        if (d_has_advective_eqn_form)
        {
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            const int num_subghosts_1_velocity = num_subghosts_velocity[1];
            const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
            
            Real* u = velocity->getPointer(0);
            Real* v = velocity->getPointer(1);
            
            std::vector<hier::IntVector> num_subghosts_conservative_var;
            num_subghosts_conservative_var.reserve(d_num_eqn);
            
            std::vector<hier::IntVector> subghostcell_dims_conservative_var;
            subghostcell_dims_conservative_var.reserve(d_num_eqn);
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
                d_flow_model->getCellDataOfConservativeVariables();
            
            std::vector<Real*> Q;
            Q.reserve(d_num_eqn);
            
            int count_eqn = 0;
            
            for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
            {
                int depth = conservative_variables[vi]->getDepth();
                
                for (int di = 0; di < depth; di++)
                {
                    // If the last element of the conservative variable vector is not in the system of equations,
                    // ignore it.
                    if (count_eqn >= d_num_eqn)
                        break;
                    
                    Q.push_back(conservative_variables[vi]->getPointer(di));
                    num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_variables[vi]->getGhostBox().numberCells());
                    
                    count_eqn++;
                }
            }
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                    const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                    const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                    
                    if (d_order == 2)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                    (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                                
                                const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                                    (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                                    (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_nghost = i + j*interior_dim_0;
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (
                                    a_n*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L])
                                    )/Real(dx[0]) +
                                    (
                                    a_n*(v[idx_cell_wghost_y_T] - v[idx_cell_wghost_y_B])
                                    )/Real(dx[1])
                                    );
                            }
                        }
                    }
                    else if (d_order == 4)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                    (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                                
                                const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_L  = (i - 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_R  = (i + 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BB = (i + num_subghosts_0_velocity) +
                                    (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_B  = (i + num_subghosts_0_velocity) +
                                    (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_T  = (i + num_subghosts_0_velocity) +
                                    (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TT = (i + num_subghosts_0_velocity) +
                                    (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_nghost = i + j*interior_dim_0;
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (
                                    a_n*(u[idx_cell_wghost_x_R]  - u[idx_cell_wghost_x_L]) +
                                    b_n*(u[idx_cell_wghost_x_RR] - u[idx_cell_wghost_x_LL])
                                    )/Real(dx[0]) +
                                    (
                                    a_n*(v[idx_cell_wghost_y_T]  - v[idx_cell_wghost_y_B]) +
                                    b_n*(v[idx_cell_wghost_y_TT] - v[idx_cell_wghost_y_BB])
                                    )/Real(dx[1])
                                    );
                            }
                        }
                    }
                    else if (d_order == 6)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                    (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                                
                                const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LL  = (i - 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_L   = (i - 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_R   = (i + 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RR  = (i + 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_velocity) +
                                    (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BB  = (i + num_subghosts_0_velocity) +
                                    (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_B   = (i + num_subghosts_0_velocity) +
                                    (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_T   = (i + num_subghosts_0_velocity) +
                                    (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TT  = (i + num_subghosts_0_velocity) +
                                    (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_velocity) +
                                    (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_nghost = i + j*interior_dim_0;
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (
                                    a_n*(u[idx_cell_wghost_x_R]   - u[idx_cell_wghost_x_L]) +
                                    b_n*(u[idx_cell_wghost_x_RR]  - u[idx_cell_wghost_x_LL]) +
                                    c_n*(u[idx_cell_wghost_x_RRR] - u[idx_cell_wghost_x_LLL])
                                    )/Real(dx[0]) +
                                    (
                                    a_n*(v[idx_cell_wghost_y_T]   - v[idx_cell_wghost_y_B]) +
                                    b_n*(v[idx_cell_wghost_y_TT]  - v[idx_cell_wghost_y_BB]) +
                                    c_n*(v[idx_cell_wghost_y_TTT] - v[idx_cell_wghost_y_BBB])
                                    )/Real(dx[1])
                                    );
                            }
                        }
                    }
                    else if (d_order == 8)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                    (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                                
                                const int idx_cell_wghost_x_LLLL = (i - 4 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LLL  = (i - 3 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LL   = (i - 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_L    = (i - 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_R    = (i + 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RR   = (i + 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRR  = (i + 3 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRRR = (i + 4 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBBB = (i + num_subghosts_0_velocity) +
                                    (j - 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBB  = (i + num_subghosts_0_velocity) +
                                    (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BB   = (i + num_subghosts_0_velocity) +
                                    (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_B    = (i + num_subghosts_0_velocity) +
                                    (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_T    = (i + num_subghosts_0_velocity) +
                                    (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TT   = (i + num_subghosts_0_velocity) +
                                    (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTT  = (i + num_subghosts_0_velocity) +
                                    (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTTT = (i + num_subghosts_0_velocity) +
                                    (j + 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_nghost = i + j*interior_dim_0;
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (
                                    a_n*(u[idx_cell_wghost_x_R]    - u[idx_cell_wghost_x_L]) +
                                    b_n*(u[idx_cell_wghost_x_RR]   - u[idx_cell_wghost_x_LL]) +
                                    c_n*(u[idx_cell_wghost_x_RRR]  - u[idx_cell_wghost_x_LLL]) +
                                    d_n*(u[idx_cell_wghost_x_RRRR] - u[idx_cell_wghost_x_LLLL])
                                    )/Real(dx[0]) +
                                    (
                                    a_n*(v[idx_cell_wghost_y_T]    - v[idx_cell_wghost_y_B]) +
                                    b_n*(v[idx_cell_wghost_y_TT]   - v[idx_cell_wghost_y_BB]) +
                                    c_n*(v[idx_cell_wghost_y_TTT]  - v[idx_cell_wghost_y_BBB]) +
                                    d_n*(v[idx_cell_wghost_y_TTTT] - v[idx_cell_wghost_y_BBBB])
                                    )/Real(dx[1])
                                    );
                            }
                        }
                    }
                    else if (d_order == 10)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                    (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                                
                                const int idx_cell_wghost_x_LLLLL = (i - 5 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LLLL  = (i - 4 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LLL   = (i - 3 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LL    = (i - 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_L     = (i - 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_R     = (i + 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RR    = (i + 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRR   = (i + 3 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRRR  = (i + 4 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRRRR = (i + 5 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBBBB = (i + num_subghosts_0_velocity) +
                                    (j - 5 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBBB  = (i + num_subghosts_0_velocity) +
                                    (j - 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBB   = (i + num_subghosts_0_velocity) +
                                    (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BB    = (i + num_subghosts_0_velocity) +
                                    (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_B     = (i + num_subghosts_0_velocity) +
                                    (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_T     = (i + num_subghosts_0_velocity) +
                                    (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TT    = (i + num_subghosts_0_velocity) +
                                    (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTT   = (i + num_subghosts_0_velocity) +
                                    (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTTT  = (i + num_subghosts_0_velocity) +
                                    (j + 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTTTT = (i + num_subghosts_0_velocity) +
                                    (j + 5 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_nghost = i + j*interior_dim_0;
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (
                                    a_n*(u[idx_cell_wghost_x_R]     - u[idx_cell_wghost_x_L]) +
                                    b_n*(u[idx_cell_wghost_x_RR]    - u[idx_cell_wghost_x_LL]) +
                                    c_n*(u[idx_cell_wghost_x_RRR]   - u[idx_cell_wghost_x_LLL]) +
                                    d_n*(u[idx_cell_wghost_x_RRRR]  - u[idx_cell_wghost_x_LLLL]) +
                                    e_n*(u[idx_cell_wghost_x_RRRRR] - u[idx_cell_wghost_x_LLLLL])
                                    )/Real(dx[0]) +
                                    (
                                    a_n*(v[idx_cell_wghost_y_T]     - v[idx_cell_wghost_y_B]) +
                                    b_n*(v[idx_cell_wghost_y_TT]    - v[idx_cell_wghost_y_BB]) +
                                    c_n*(v[idx_cell_wghost_y_TTT]   - v[idx_cell_wghost_y_BBB]) +
                                    d_n*(v[idx_cell_wghost_y_TTTT]  - v[idx_cell_wghost_y_BBBB]) +
                                    e_n*(v[idx_cell_wghost_y_TTTTT] - v[idx_cell_wghost_y_BBBBB])
                                    )/Real(dx[1])
                                    );
                            }
                        }
                    }
                    else if (d_order == 12)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                    (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                                
                                const int idx_cell_wghost_x_LLLLLL = (i - 6 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LLLLL  = (i - 5 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LLLL   = (i - 4 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LLL    = (i - 3 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_LL     = (i - 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_L      = (i - 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_R      = (i + 1 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RR     = (i + 2 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRR    = (i + 3 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRRR   = (i + 4 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRRRR  = (i + 5 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_x_RRRRRR = (i + 6 + num_subghosts_0_velocity) +
                                    (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBBBBB = (i + num_subghosts_0_velocity) +
                                    (j - 6 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBBBB  = (i + num_subghosts_0_velocity) +
                                    (j - 5 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBBB   = (i + num_subghosts_0_velocity) +
                                    (j - 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BBB    = (i + num_subghosts_0_velocity) +
                                    (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_BB     = (i + num_subghosts_0_velocity) +
                                    (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_B      = (i + num_subghosts_0_velocity) +
                                    (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_T      = (i + num_subghosts_0_velocity) +
                                    (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TT     = (i + num_subghosts_0_velocity) +
                                    (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTT    = (i + num_subghosts_0_velocity) +
                                    (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTTT   = (i + num_subghosts_0_velocity) +
                                    (j + 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTTTT  = (i + num_subghosts_0_velocity) +
                                    (j + 5 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_wghost_y_TTTTTT = (i + num_subghosts_0_velocity) +
                                    (j + 6 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                                
                                const int idx_cell_nghost = i + j*interior_dim_0;
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (
                                    a_n*(u[idx_cell_wghost_x_R]      - u[idx_cell_wghost_x_L]) +
                                    b_n*(u[idx_cell_wghost_x_RR]     - u[idx_cell_wghost_x_LL]) +
                                    c_n*(u[idx_cell_wghost_x_RRR]    - u[idx_cell_wghost_x_LLL]) +
                                    d_n*(u[idx_cell_wghost_x_RRRR]   - u[idx_cell_wghost_x_LLLL]) +
                                    e_n*(u[idx_cell_wghost_x_RRRRR]  - u[idx_cell_wghost_x_LLLLL]) +
                                    f_n*(u[idx_cell_wghost_x_RRRRRR] - u[idx_cell_wghost_x_LLLLLL])
                                    )/Real(dx[0]) +
                                    (
                                    a_n*(v[idx_cell_wghost_y_T]      - v[idx_cell_wghost_y_B]) +
                                    b_n*(v[idx_cell_wghost_y_TT]     - v[idx_cell_wghost_y_BB]) +
                                    c_n*(v[idx_cell_wghost_y_TTT]    - v[idx_cell_wghost_y_BBB]) +
                                    d_n*(v[idx_cell_wghost_y_TTTT]   - v[idx_cell_wghost_y_BBBB]) +
                                    e_n*(v[idx_cell_wghost_y_TTTTT]  - v[idx_cell_wghost_y_BBBBB]) +
                                    f_n*(v[idx_cell_wghost_y_TTTTTT] - v[idx_cell_wghost_y_BBBBBB])
                                    )/Real(dx[1])
                                    );
                            }
                        }
                    }
                }
            }
        }
        
        t_compute_source->stop();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Z", d_num_conv_ghosts));
        
        if (d_has_advective_eqn_form)
        {
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        }
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(3);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_z = convective_flux_node[2]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_z = convective_flux_node[2]->getGhostBox().numberCells();
        
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
        const int num_subghosts_2_convective_flux_x = num_subghosts_convective_flux_x[2];
        const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
        const int subghostcell_dim_1_convective_flux_x = subghostcell_dims_convective_flux_x[1];
        
        const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
        const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
        const int num_subghosts_2_convective_flux_y = num_subghosts_convective_flux_y[2];
        const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
        const int subghostcell_dim_1_convective_flux_y = subghostcell_dims_convective_flux_y[1];
        
        const int num_subghosts_0_convective_flux_z = num_subghosts_convective_flux_z[0];
        const int num_subghosts_1_convective_flux_z = num_subghosts_convective_flux_z[1];
        const int num_subghosts_2_convective_flux_z = num_subghosts_convective_flux_z[2];
        const int subghostcell_dim_0_convective_flux_z = subghostcell_dims_convective_flux_z[0];
        const int subghostcell_dim_1_convective_flux_z = subghostcell_dims_convective_flux_z[1];
        
        std::vector<Real*> F_node_x;
        std::vector<Real*> F_node_y;
        std::vector<Real*> F_node_z;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        F_node_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
            F_node_z.push_back(convective_flux_node[2]->getPointer(ei));
        }
        
        t_reconstruct_flux->start();
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        if (d_order == 2)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = Real(dt)*(
                                a_m*(F_node_x[ei][idx_node_L] + F_node_x[ei][idx_node_R])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 4)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_LL = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_L  = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R  = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RR = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = Real(dt)*(
                                a_m*(F_node_x[ei][idx_node_L]  + F_node_x[ei][idx_node_R]) +
                                b_m*(F_node_x[ei][idx_node_LL] + F_node_x[ei][idx_node_RR])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 6)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_LLL = (i - 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LL  = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_L   = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R   = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RR  = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRR = (i + 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = Real(dt)*(
                                a_m*(F_node_x[ei][idx_node_L]   + F_node_x[ei][idx_node_R]) +
                                b_m*(F_node_x[ei][idx_node_LL]  + F_node_x[ei][idx_node_RR]) +
                                c_m*(F_node_x[ei][idx_node_LLL] + F_node_x[ei][idx_node_RRR])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 8)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_LLLL = (i - 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLL  = (i - 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LL   = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_L    = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R    = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RR   = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRR  = (i + 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRR = (i + 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = Real(dt)*(
                                a_m*(F_node_x[ei][idx_node_L]    + F_node_x[ei][idx_node_R]) +
                                b_m*(F_node_x[ei][idx_node_LL]   + F_node_x[ei][idx_node_RR]) +
                                c_m*(F_node_x[ei][idx_node_LLL]  + F_node_x[ei][idx_node_RRR]) +
                                d_m*(F_node_x[ei][idx_node_LLLL] + F_node_x[ei][idx_node_RRRR])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 10)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_LLLLL = (i - 5 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLLL  = (i - 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLL   = (i - 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LL    = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_L     = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R     = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RR    = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRR   = (i + 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRR  = (i + 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRRR = (i + 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = Real(dt)*(
                                a_m*(F_node_x[ei][idx_node_L]     + F_node_x[ei][idx_node_R]) +
                                b_m*(F_node_x[ei][idx_node_LL]    + F_node_x[ei][idx_node_RR]) +
                                c_m*(F_node_x[ei][idx_node_LLL]   + F_node_x[ei][idx_node_RRR]) +
                                d_m*(F_node_x[ei][idx_node_LLLL]  + F_node_x[ei][idx_node_RRRR]) +
                                e_m*(F_node_x[ei][idx_node_LLLLL] + F_node_x[ei][idx_node_RRRRR])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 12)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_x = convective_flux->getPointer(0, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_LLLLLL = (i - 6 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLLLL  = (i - 5 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLLL   = (i - 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLL    = (i - 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LL     = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_L      = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R      = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RR     = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRR    = (i + 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRR   = (i + 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRRR  = (i + 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRRRR = (i + 5 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = Real(dt)*(
                                a_m*(F_node_x[ei][idx_node_L]      + F_node_x[ei][idx_node_R]) +
                                b_m*(F_node_x[ei][idx_node_LL]     + F_node_x[ei][idx_node_RR]) +
                                c_m*(F_node_x[ei][idx_node_LLL]    + F_node_x[ei][idx_node_RRR]) +
                                d_m*(F_node_x[ei][idx_node_LLLL]   + F_node_x[ei][idx_node_RRRR]) +
                                e_m*(F_node_x[ei][idx_node_LLLLL]  + F_node_x[ei][idx_node_RRRRR]) +
                                f_m*(F_node_x[ei][idx_node_LLLLLL] + F_node_x[ei][idx_node_RRRRRR])
                                );
                        }
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        if (d_order == 2)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = Real(dt)*(
                                a_m*(F_node_y[ei][idx_node_B] + F_node_y[ei][idx_node_T])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 4)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_BB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_B  = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T  = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = Real(dt)*(
                                a_m*(F_node_y[ei][idx_node_B]  + F_node_y[ei][idx_node_T]) +
                                b_m*(F_node_y[ei][idx_node_BB] + F_node_y[ei][idx_node_TT])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 6)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_BBB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BB  = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_B   = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T   = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TT  = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = Real(dt)*(
                                a_m*(F_node_y[ei][idx_node_B]   + F_node_y[ei][idx_node_T]) +
                                b_m*(F_node_y[ei][idx_node_BB]  + F_node_y[ei][idx_node_TT]) +
                                c_m*(F_node_y[ei][idx_node_BBB] + F_node_y[ei][idx_node_TTT])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 8)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_BBBB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBB  = (i + num_subghosts_0_convective_flux_y) +
                                (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BB   = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_B    = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T    = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TT   = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTT  = (i + num_subghosts_0_convective_flux_y) +
                                (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = Real(dt)*(
                                a_m*(F_node_y[ei][idx_node_B]    + F_node_y[ei][idx_node_T]) +
                                b_m*(F_node_y[ei][idx_node_BB]   + F_node_y[ei][idx_node_TT]) +
                                c_m*(F_node_y[ei][idx_node_BBB]  + F_node_y[ei][idx_node_TTT]) +
                                d_m*(F_node_y[ei][idx_node_BBBB] + F_node_y[ei][idx_node_TTTT])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 10)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_BBBBB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBBB  = (i + num_subghosts_0_convective_flux_y) +
                                (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBB   = (i + num_subghosts_0_convective_flux_y) +
                                (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BB    = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_B     = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T     = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TT    = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTT   = (i + num_subghosts_0_convective_flux_y) +
                                (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTT  = (i + num_subghosts_0_convective_flux_y) +
                                (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTTT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = Real(dt)*(
                                a_m*(F_node_y[ei][idx_node_B]     + F_node_y[ei][idx_node_T]) +
                                b_m*(F_node_y[ei][idx_node_BB]    + F_node_y[ei][idx_node_TT]) +
                                c_m*(F_node_y[ei][idx_node_BBB]   + F_node_y[ei][idx_node_TTT]) +
                                d_m*(F_node_y[ei][idx_node_BBBB]  + F_node_y[ei][idx_node_TTTT]) +
                                e_m*(F_node_y[ei][idx_node_BBBBB] + F_node_y[ei][idx_node_TTTTT])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 12)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_y = convective_flux->getPointer(1, ei);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_BBBBBB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 6 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBBBB  = (i + num_subghosts_0_convective_flux_y) +
                                (j - 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBBB   = (i + num_subghosts_0_convective_flux_y) +
                                (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBB    = (i + num_subghosts_0_convective_flux_y) +
                                (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BB     = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_B      = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T      = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TT     = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTT    = (i + num_subghosts_0_convective_flux_y) +
                                (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTT   = (i + num_subghosts_0_convective_flux_y) +
                                (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTTT  = (i + num_subghosts_0_convective_flux_y) +
                                (j + 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTTTT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = Real(dt)*(
                                a_m*(F_node_y[ei][idx_node_B]      + F_node_y[ei][idx_node_T]) +
                                b_m*(F_node_y[ei][idx_node_BB]     + F_node_y[ei][idx_node_TT]) +
                                c_m*(F_node_y[ei][idx_node_BBB]    + F_node_y[ei][idx_node_TTT]) +
                                d_m*(F_node_y[ei][idx_node_BBBB]   + F_node_y[ei][idx_node_TTTT]) +
                                e_m*(F_node_y[ei][idx_node_BBBBB]  + F_node_y[ei][idx_node_TTTTT]) +
                                f_m*(F_node_y[ei][idx_node_BBBBBB] + F_node_y[ei][idx_node_TTTTTT])
                                );
                        }
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the z-direction.
         */
        
        if (d_order == 2)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_z = convective_flux->getPointer(2, ei);
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_B = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = Real(dt)*(
                                a_m*(F_node_z[ei][idx_node_B] + F_node_z[ei][idx_node_F])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 4)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_z = convective_flux->getPointer(2, ei);
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_BB = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_B  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FF = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = Real(dt)*(
                                a_m*(F_node_z[ei][idx_node_B]  + F_node_z[ei][idx_node_F]) +
                                b_m*(F_node_z[ei][idx_node_BB] + F_node_z[ei][idx_node_FF])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 6)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_z = convective_flux->getPointer(2, ei);
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_BBB = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BB  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_B   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FF  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFF = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = Real(dt)*(
                                a_m*(F_node_z[ei][idx_node_B]   + F_node_z[ei][idx_node_F]) +
                                b_m*(F_node_z[ei][idx_node_BB]  + F_node_z[ei][idx_node_FF]) +
                                c_m*(F_node_z[ei][idx_node_BBB] + F_node_z[ei][idx_node_FFF])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 8)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_z = convective_flux->getPointer(2, ei);
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_BBBB = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBB  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BB   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_B    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FF   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFF  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFF = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = Real(dt)*(
                                a_m*(F_node_z[ei][idx_node_B]    + F_node_z[ei][idx_node_F]) +
                                b_m*(F_node_z[ei][idx_node_BB]   + F_node_z[ei][idx_node_FF]) +
                                c_m*(F_node_z[ei][idx_node_BBB]  + F_node_z[ei][idx_node_FFF]) +
                                d_m*(F_node_z[ei][idx_node_BBBB] + F_node_z[ei][idx_node_FFFF])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 10)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_z = convective_flux->getPointer(2, ei);
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_BBBBB = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 5 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBBB  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBB   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BB    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_B     = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F     = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FF    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFF   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFF  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFFF = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = Real(dt)*(
                                a_m*(F_node_z[ei][idx_node_B]     + F_node_z[ei][idx_node_F]) +
                                b_m*(F_node_z[ei][idx_node_BB]    + F_node_z[ei][idx_node_FF]) +
                                c_m*(F_node_z[ei][idx_node_BBB]   + F_node_z[ei][idx_node_FFF]) +
                                d_m*(F_node_z[ei][idx_node_BBBB]  + F_node_z[ei][idx_node_FFFF]) +
                                e_m*(F_node_z[ei][idx_node_BBBBB] + F_node_z[ei][idx_node_FFFFF])
                                );
                        }
                    }
                }
            }
        }
        else if (d_order == 12)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                Real* F_face_z = convective_flux->getPointer(2, ei);
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_BBBBBB = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 6 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBBBB  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 5 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBBB   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBB    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BB     = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_B      = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F      = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FF     = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFF    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFF   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFFF  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFFFF = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 5 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = Real(dt)*(
                                a_m*(F_node_z[ei][idx_node_B]      + F_node_z[ei][idx_node_F]) +
                                b_m*(F_node_z[ei][idx_node_BB]     + F_node_z[ei][idx_node_FF]) +
                                c_m*(F_node_z[ei][idx_node_BBB]    + F_node_z[ei][idx_node_FFF]) +
                                d_m*(F_node_z[ei][idx_node_BBBB]   + F_node_z[ei][idx_node_FFFF]) +
                                e_m*(F_node_z[ei][idx_node_BBBBB]  + F_node_z[ei][idx_node_FFFFF]) +
                                f_m*(F_node_z[ei][idx_node_BBBBBB] + F_node_z[ei][idx_node_FFFFFF])
                                );
                        }
                    }
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Compute the source.
         */
        
        t_compute_source->start();
        
        if (d_has_advective_eqn_form)
        {
            HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
            
            hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
            hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
            
            const int num_subghosts_0_velocity = num_subghosts_velocity[0];
            const int num_subghosts_1_velocity = num_subghosts_velocity[1];
            const int num_subghosts_2_velocity = num_subghosts_velocity[2];
            const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
            const int subghostcell_dim_1_velocity = subghostcell_dims_velocity[1];
            
            Real* u = velocity->getPointer(0);
            Real* v = velocity->getPointer(1);
            Real* w = velocity->getPointer(2);
            
            std::vector<hier::IntVector> num_subghosts_conservative_var;
            num_subghosts_conservative_var.reserve(d_num_eqn);
            
            std::vector<hier::IntVector> subghostcell_dims_conservative_var;
            subghostcell_dims_conservative_var.reserve(d_num_eqn);
            
            std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
                d_flow_model->getCellDataOfConservativeVariables();
            
            std::vector<Real*> Q;
            Q.reserve(d_num_eqn);
            
            int count_eqn = 0;
            
            for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
            {
                int depth = conservative_variables[vi]->getDepth();
                
                for (int di = 0; di < depth; di++)
                {
                    // If the last element of the conservative variable vector is not in the system of equations,
                    // ignore it.
                    if (count_eqn >= d_num_eqn)
                        break;
                    
                    Q.push_back(conservative_variables[vi]->getPointer(di));
                    num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_variables[vi]->getGhostBox().numberCells());
                    
                    count_eqn++;
                }
            }
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                    const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                    const int num_subghosts_2_conservative_var = num_subghosts_conservative_var[ei][2];
                    const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                    const int subghostcell_dim_1_conservative_var = subghostcell_dims_conservative_var[ei][1];
                    
                    if (d_order == 2)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                        (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                        (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                            subghostcell_dim_1_conservative_var;
                                    
                                    const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                                        (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                                        (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_B = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_F = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_nghost = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*
                                            interior_dim_1;
                                    
                                    S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                        (
                                        a_n*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L])
                                        )/Real(dx[0]) +
                                        (
                                        a_n*(v[idx_cell_wghost_y_T] - v[idx_cell_wghost_y_B])
                                        )/Real(dx[1]) +
                                        (
                                        a_n*(w[idx_cell_wghost_z_F] - w[idx_cell_wghost_z_B])
                                        )/Real(dx[2])
                                        );
                                }
                            }
                        }
                    }
                    else if (d_order == 4)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                        (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                        (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                            subghostcell_dim_1_conservative_var;
                                    
                                    const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_L  = (i - 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_R  = (i + 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BB = (i + num_subghosts_0_velocity) +
                                        (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_B  = (i + num_subghosts_0_velocity) +
                                        (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_T  = (i + num_subghosts_0_velocity) +
                                        (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TT = (i + num_subghosts_0_velocity) +
                                        (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BB = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_B  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_F  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FF = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_nghost = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*
                                            interior_dim_1;
                                    
                                    S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                        (
                                        a_n*(u[idx_cell_wghost_x_R]  - u[idx_cell_wghost_x_L]) +
                                        b_n*(u[idx_cell_wghost_x_RR] - u[idx_cell_wghost_x_LL])
                                        )/Real(dx[0]) +
                                        (
                                        a_n*(v[idx_cell_wghost_y_T]  - v[idx_cell_wghost_y_B]) +
                                        b_n*(v[idx_cell_wghost_y_TT] - v[idx_cell_wghost_y_BB])
                                        )/Real(dx[1]) +
                                        (
                                        a_n*(w[idx_cell_wghost_z_F]  - w[idx_cell_wghost_z_B]) +
                                        b_n*(w[idx_cell_wghost_z_FF] - w[idx_cell_wghost_z_BB])
                                        )/Real(dx[2])
                                        );
                                }
                            }
                        }
                    }
                    else if (d_order == 6)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                        (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                        (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                            subghostcell_dim_1_conservative_var;
                                    
                                    const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LL  = (i - 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_L   = (i - 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_R   = (i + 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RR  = (i + 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_velocity) +
                                        (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BB  = (i + num_subghosts_0_velocity) +
                                        (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_B   = (i + num_subghosts_0_velocity) +
                                        (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_T   = (i + num_subghosts_0_velocity) +
                                        (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TT  = (i + num_subghosts_0_velocity) +
                                        (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_velocity) +
                                        (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBB = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BB  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_B   = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_F   = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FF  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFF = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_nghost = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*
                                            interior_dim_1;
                                    
                                    S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                        (
                                        a_n*(u[idx_cell_wghost_x_R]   - u[idx_cell_wghost_x_L]) +
                                        b_n*(u[idx_cell_wghost_x_RR]  - u[idx_cell_wghost_x_LL]) +
                                        c_n*(u[idx_cell_wghost_x_RRR] - u[idx_cell_wghost_x_LLL])
                                        )/Real(dx[0]) +
                                        (
                                        a_n*(v[idx_cell_wghost_y_T]   - v[idx_cell_wghost_y_B]) +
                                        b_n*(v[idx_cell_wghost_y_TT]  - v[idx_cell_wghost_y_BB]) +
                                        c_n*(v[idx_cell_wghost_y_TTT] - v[idx_cell_wghost_y_BBB])
                                        )/Real(dx[1]) +
                                        (
                                        a_n*(w[idx_cell_wghost_z_F]   - w[idx_cell_wghost_z_B]) +
                                        b_n*(w[idx_cell_wghost_z_FF]  - w[idx_cell_wghost_z_BB]) +
                                        c_n*(w[idx_cell_wghost_z_FFF] - w[idx_cell_wghost_z_BBB])
                                        )/Real(dx[2])
                                        );
                                }
                            }
                        }
                    }
                    else if (d_order == 8)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                        (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                        (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                            subghostcell_dim_1_conservative_var;
                                    
                                    const int idx_cell_wghost_x_LLLL = (i - 4 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LLL  = (i - 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LL   = (i - 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_L    = (i - 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_R    = (i + 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RR   = (i + 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRR  = (i + 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRRR = (i + 4 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBBB = (i + num_subghosts_0_velocity) +
                                        (j - 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBB  = (i + num_subghosts_0_velocity) +
                                        (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BB   = (i + num_subghosts_0_velocity) +
                                        (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_B    = (i + num_subghosts_0_velocity) +
                                        (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_T    = (i + num_subghosts_0_velocity) +
                                        (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TT   = (i + num_subghosts_0_velocity) +
                                        (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTT  = (i + num_subghosts_0_velocity) +
                                        (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTTT = (i + num_subghosts_0_velocity) +
                                        (j + 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBBB = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 4 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBB  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BB   = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_B    = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_F    = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FF   = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFF  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFFF = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 4 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_nghost = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*
                                            interior_dim_1;
                                    
                                    S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                        (
                                        a_n*(u[idx_cell_wghost_x_R]    - u[idx_cell_wghost_x_L]) +
                                        b_n*(u[idx_cell_wghost_x_RR]   - u[idx_cell_wghost_x_LL]) +
                                        c_n*(u[idx_cell_wghost_x_RRR]  - u[idx_cell_wghost_x_LLL]) +
                                        d_n*(u[idx_cell_wghost_x_RRRR] - u[idx_cell_wghost_x_LLLL])
                                        )/Real(dx[0]) +
                                        (
                                        a_n*(v[idx_cell_wghost_y_T]    - v[idx_cell_wghost_y_B]) +
                                        b_n*(v[idx_cell_wghost_y_TT]   - v[idx_cell_wghost_y_BB]) +
                                        c_n*(v[idx_cell_wghost_y_TTT]  - v[idx_cell_wghost_y_BBB]) +
                                        d_n*(v[idx_cell_wghost_y_TTTT] - v[idx_cell_wghost_y_BBBB])
                                        )/Real(dx[1]) +
                                        (
                                        a_n*(w[idx_cell_wghost_z_F]    - w[idx_cell_wghost_z_B]) +
                                        b_n*(w[idx_cell_wghost_z_FF]   - w[idx_cell_wghost_z_BB]) +
                                        c_n*(w[idx_cell_wghost_z_FFF]  - w[idx_cell_wghost_z_BBB]) +
                                        d_n*(w[idx_cell_wghost_z_FFFF] - w[idx_cell_wghost_z_BBBB])
                                        )/Real(dx[2])
                                        );
                                }
                            }
                        }
                    }
                    else if (d_order == 10)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                        (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                        (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                            subghostcell_dim_1_conservative_var;
                                    
                                    const int idx_cell_wghost_x_LLLLL = (i - 5 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LLLL  = (i - 4 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LLL   = (i - 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LL    = (i - 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_L     = (i - 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_R     = (i + 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RR    = (i + 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRR   = (i + 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRRR  = (i + 4 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRRRR = (i + 5 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBBBB = (i + num_subghosts_0_velocity) +
                                        (j - 5 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBBB  = (i + num_subghosts_0_velocity) +
                                        (j - 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBB   = (i + num_subghosts_0_velocity) +
                                        (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BB    = (i + num_subghosts_0_velocity) +
                                        (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_B     = (i + num_subghosts_0_velocity) +
                                        (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_T     = (i + num_subghosts_0_velocity) +
                                        (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TT    = (i + num_subghosts_0_velocity) +
                                        (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTT   = (i + num_subghosts_0_velocity) +
                                        (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTTT  = (i + num_subghosts_0_velocity) +
                                        (j + 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTTTT = (i + num_subghosts_0_velocity) +
                                        (j + 5 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBBBB = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 5 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBBB  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 4 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBB   = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BB    = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_B     = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_F     = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FF    = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFF   = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFFF  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 4 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFFFF = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 5 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_nghost = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*
                                            interior_dim_1;
                                    
                                    S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                        (
                                        a_n*(u[idx_cell_wghost_x_R]     - u[idx_cell_wghost_x_L]) +
                                        b_n*(u[idx_cell_wghost_x_RR]    - u[idx_cell_wghost_x_LL]) +
                                        c_n*(u[idx_cell_wghost_x_RRR]   - u[idx_cell_wghost_x_LLL]) +
                                        d_n*(u[idx_cell_wghost_x_RRRR]  - u[idx_cell_wghost_x_LLLL]) +
                                        e_n*(u[idx_cell_wghost_x_RRRRR] - u[idx_cell_wghost_x_LLLLL])
                                        )/Real(dx[0]) +
                                        (
                                        a_n*(v[idx_cell_wghost_y_T]     - v[idx_cell_wghost_y_B]) +
                                        b_n*(v[idx_cell_wghost_y_TT]    - v[idx_cell_wghost_y_BB]) +
                                        c_n*(v[idx_cell_wghost_y_TTT]   - v[idx_cell_wghost_y_BBB]) +
                                        d_n*(v[idx_cell_wghost_y_TTTT]  - v[idx_cell_wghost_y_BBBB]) +
                                        e_n*(v[idx_cell_wghost_y_TTTTT] - v[idx_cell_wghost_y_BBBBB])
                                        )/Real(dx[1]) +
                                        (
                                        a_n*(w[idx_cell_wghost_z_F]     - w[idx_cell_wghost_z_B]) +
                                        b_n*(w[idx_cell_wghost_z_FF]    - w[idx_cell_wghost_z_BB]) +
                                        c_n*(w[idx_cell_wghost_z_FFF]   - w[idx_cell_wghost_z_BBB]) +
                                        d_n*(w[idx_cell_wghost_z_FFFF]  - w[idx_cell_wghost_z_BBBB]) +
                                        e_n*(w[idx_cell_wghost_z_FFFFF] - w[idx_cell_wghost_z_BBBBB])
                                        )/Real(dx[2])
                                        );
                                }
                            }
                        }
                    }
                    else if (d_order == 12)
                    {
                       for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                HAMERS_PRAGMA_SIMD
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute the linear indices.
                                    const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                        (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                        (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                            subghostcell_dim_1_conservative_var;
                                    
                                    const int idx_cell_wghost_x_LLLLLL = (i - 6 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LLLLL  = (i - 5 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LLLL   = (i - 4 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LLL    = (i - 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_LL     = (i - 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_L      = (i - 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_R      = (i + 1 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RR     = (i + 2 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRR    = (i + 3 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRRR   = (i + 4 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRRRR  = (i + 5 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_x_RRRRRR = (i + 6 + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBBBBB = (i + num_subghosts_0_velocity) +
                                        (j - 6 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBBBB  = (i + num_subghosts_0_velocity) +
                                        (j - 5 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBBB   = (i + num_subghosts_0_velocity) +
                                        (j - 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BBB    = (i + num_subghosts_0_velocity) +
                                        (j - 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_BB     = (i + num_subghosts_0_velocity) +
                                        (j - 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_B      = (i + num_subghosts_0_velocity) +
                                        (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_T      = (i + num_subghosts_0_velocity) +
                                        (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TT     = (i + num_subghosts_0_velocity) +
                                        (j + 2 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTT    = (i + num_subghosts_0_velocity) +
                                        (j + 3 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTTT   = (i + num_subghosts_0_velocity) +
                                        (j + 4 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTTTT  = (i + num_subghosts_0_velocity) +
                                        (j + 5 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_y_TTTTTT = (i + num_subghosts_0_velocity) +
                                        (j + 6 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBBBBB = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 6 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBBBB  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 5 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBBB   = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 4 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BBB    = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_BB     = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_B      = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_F      = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FF     = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 2 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFF    = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 3 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFFF   = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 4 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFFFF  = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 5 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_wghost_z_FFFFFF = (i + num_subghosts_0_velocity) +
                                        (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                        (k + 6 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                            subghostcell_dim_1_velocity;
                                    
                                    const int idx_cell_nghost = i +
                                        j*interior_dim_0 +
                                        k*interior_dim_0*
                                            interior_dim_1;
                                    
                                    S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                        (
                                        a_n*(u[idx_cell_wghost_x_R]      - u[idx_cell_wghost_x_L]) +
                                        b_n*(u[idx_cell_wghost_x_RR]     - u[idx_cell_wghost_x_LL]) +
                                        c_n*(u[idx_cell_wghost_x_RRR]    - u[idx_cell_wghost_x_LLL]) +
                                        d_n*(u[idx_cell_wghost_x_RRRR]   - u[idx_cell_wghost_x_LLLL]) +
                                        e_n*(u[idx_cell_wghost_x_RRRRR]  - u[idx_cell_wghost_x_LLLLL]) +
                                        f_n*(u[idx_cell_wghost_x_RRRRRR] - u[idx_cell_wghost_x_LLLLLL])
                                        )/Real(dx[0]) +
                                        (
                                        a_n*(v[idx_cell_wghost_y_T]      - v[idx_cell_wghost_y_B]) +
                                        b_n*(v[idx_cell_wghost_y_TT]     - v[idx_cell_wghost_y_BB]) +
                                        c_n*(v[idx_cell_wghost_y_TTT]    - v[idx_cell_wghost_y_BBB]) +
                                        d_n*(v[idx_cell_wghost_y_TTTT]   - v[idx_cell_wghost_y_BBBB]) +
                                        e_n*(v[idx_cell_wghost_y_TTTTT]  - v[idx_cell_wghost_y_BBBBB]) +
                                        f_n*(v[idx_cell_wghost_y_TTTTTT] - v[idx_cell_wghost_y_BBBBBB])
                                        )/Real(dx[1]) +
                                        (
                                        a_n*(w[idx_cell_wghost_z_F]      - w[idx_cell_wghost_z_B]) +
                                        b_n*(w[idx_cell_wghost_z_FF]     - w[idx_cell_wghost_z_BB]) +
                                        c_n*(w[idx_cell_wghost_z_FFF]    - w[idx_cell_wghost_z_BBB]) +
                                        d_n*(w[idx_cell_wghost_z_FFFF]   - w[idx_cell_wghost_z_BBBB]) +
                                        e_n*(w[idx_cell_wghost_z_FFFFF]  - w[idx_cell_wghost_z_BBBBB]) +
                                        f_n*(w[idx_cell_wghost_z_FFFFFF] - w[idx_cell_wghost_z_BBBBBB])
                                        )/Real(dx[2])
                                        );
                                }
                            }
                        }
                    }
                }
            }
        }
        
        t_compute_source->stop();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(3))
}


/*
 * Compute the convective flux and source due to splitting using shock-capturing scheme.
 */
void
ConvectiveFluxReconstructorCentral::computeConvectiveFluxAndSourceOnPatchShockCapturing(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const hier::Box& domain,
    const Real coeff_blending,
    const double dt) const
{
    d_flow_model->setupRiemannSolver();
    d_flow_model->setupBasicUtilities();
    
    HAMERS_SHARED_PTR<FlowModelRiemannSolver> riemann_solver = d_flow_model->getFlowModelRiemannSolver();
    HAMERS_SHARED_PTR<FlowModelBasicUtilities> basic_utilities = d_flow_model->getFlowModelBasicUtilities();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    domain_lo = domain.lower() - interior_box.lower();
    domain_dims = domain.numberCells();
    
    // Allocate temporary patch data.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity_midpoint;
    
    if (d_has_advective_eqn_form)
    {
        velocity_midpoint.reset(new pdat::SideData<Real>(
            interior_box, d_dim.getValue(), hier::IntVector::getOne(d_dim)));
    }
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux_midpoint(
        new pdat::SideData<Real>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux_midpoint_HLLC(
        new pdat::SideData<Real>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
    
    /*
     * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
    if (d_dim > tbox::Dimension(1))
    {
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
    }
    if (d_dim > tbox::Dimension(2))
    {
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Z", d_num_conv_ghosts));
    }
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
    
    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
    
    d_flow_model->allocateMemoryForDerivedCellData();
    
    d_flow_model->computeDerivedCellData();
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(d_dim.getValue());
    convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
    if (d_dim > tbox::Dimension(1))
    {
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
    }
    if (d_dim > tbox::Dimension(2))
    {
        convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
    }
    
    
    /*
     * Get the pointers to the conservative variables and primitive variables.
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
        d_flow_model->getCellDataOfConservativeVariables();
    
    std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > primitive_variables =
        d_flow_model->getCellDataOfPrimitiveVariables();
    
    std::vector<hier::IntVector> num_subghosts_conservative_var;
    num_subghosts_conservative_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> subghostcell_dims_conservative_var;
    subghostcell_dims_conservative_var.reserve(d_num_eqn);
    
    std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
    
    std::vector<Real*> Q;
    Q.reserve(d_num_eqn);
    
    std::vector<Real*> V;
    V.reserve(d_num_eqn);
    
    int count_eqn = 0;
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        int depth = conservative_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            // If the last element of the conservative variable vector is not in the system of equations,
            // ignore it.
            if (count_eqn >= d_num_eqn)
                break;
            
            Q.push_back(conservative_variables[vi]->getPointer(di));
            num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
            subghostcell_dims_conservative_var.push_back(
                conservative_variables[vi]->getGhostBox().numberCells());
            
            count_eqn++;
        }
    }
    
    count_eqn = 0;
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        int depth = primitive_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            // If the last element of the primitive variable vector is not in the system of equations,
            // ignore it.
            if (count_eqn >= d_num_eqn)
                break;
            
            V.push_back(primitive_variables[vi]->getPointer(di));
            num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
            subghostcell_dims_primitive_var.push_back(
                primitive_variables[vi]->getGhostBox().numberCells());
            
            count_eqn++;
        }
    }
    
    /*
     * Compute the convective flux and source using shock-capturing scheme.
     */
    if (d_dim == tbox::Dimension(1))
    {
        // NOT IMPLEMENTED.
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices and the number of cells in each dimension.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        /*
         * Get the interior dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Get the pointers to the velocity and convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(2);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        
        hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
        hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        const int num_subghosts_0_velocity = num_subghosts_velocity[0];
        const int num_subghosts_1_velocity = num_subghosts_velocity[1];
        const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
        
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
        const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
        
        const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
        const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
        const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
        
        Real* u = velocity->getPointer(0);
        Real* v = velocity->getPointer(1);
        
        std::vector<Real*> F_node_x;
        std::vector<Real*> F_node_y;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
        }
        
        std::vector<Real*> F_midpoint_x;
        std::vector<Real*> F_midpoint_y;
        std::vector<Real*> F_midpoint_HLLC_x;
        std::vector<Real*> F_midpoint_HLLC_y;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        F_midpoint_HLLC_x.reserve(d_num_eqn);
        F_midpoint_HLLC_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
        }
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_HLLC_x.push_back(convective_flux_midpoint_HLLC->getPointer(0, ei));
            F_midpoint_HLLC_y.push_back(convective_flux_midpoint_HLLC->getPointer(1, ei));
        }
        
        /*
         * Declare temporary data containers for WENO interpolation.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > projection_variables;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > characteristic_variables;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_plus;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_plus;
        
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_minus;
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
         */
        
        const int num_projection_var = basic_utilities->getNumberOfProjectionVariablesForPrimitiveVariables();
        projection_variables.reserve(num_projection_var);
        
        for (int vi = 0; vi < num_projection_var; vi++)
        {
            projection_variables.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        characteristic_variables.resize(6);
        
        for (int m = 0; m < 6; m++)
        {
            characteristic_variables[m].reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                characteristic_variables[m].push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                    interior_box, 1, hier::IntVector::getOne(d_dim)));
            }
        }
        
        characteristic_variables_minus.reserve(d_num_eqn);
        characteristic_variables_plus.reserve(d_num_eqn);
        primitive_variables_minus.reserve(d_num_eqn);
        primitive_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            characteristic_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            characteristic_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        
        basic_utilities->computeSideDataOfProjectionVariablesForPrimitiveVariables(
            projection_variables);
        
        /*
         * Transform primitive variables to characteristic variables.
         */
        
        for (int m = 0; m < 6; m++)
        {
            basic_utilities->computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
                characteristic_variables[m],
                primitive_variables,
                projection_variables,
                m - 3);
        }
        
        /*
         * Peform WENO interpolation.
         */
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables,
            domain);
        
        /*
         * Transform characteristic variables back to primitive variables.
         */
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_minus,
            characteristic_variables_minus,
            projection_variables);
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_plus,
            characteristic_variables_plus,
            projection_variables);
        
        /*
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<Real*> V_minus;
        std::vector<Real*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus);
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_plus,
            primitive_variables_plus);
        
        /*
         * Use first order interpolation if interpolated side primitive variables in x-direction
         * are out of bounds.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
        flag_minus = bounded_flag_minus->getPointer(0);
        flag_plus = bounded_flag_plus->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                        (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                    
                    const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                        (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                    
                    if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                    {
                        V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                        V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                    }
                }
            }
        }
        
        /*
         * Use first order interpolation if interpolated side primitive variables in y-direction
         * are out of bounds.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
        }
        
        flag_minus = bounded_flag_minus->getPointer(1);
        flag_plus = bounded_flag_plus->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            
            for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                        (j - 1 + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                    
                    const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                        (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                    
                    if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0)
                    {
                        V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                        V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                    }
                }
            }
        }
        
        /*
         * Compute mid-point flux in the x-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                velocity_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    F_midpoint_x[ei][idx_midpoint_x] = F_midpoint_HLLC_x[ei][idx_midpoint_x];
                }
            }
        }
        
        /*
         * Compute mid-point flux in the y-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                velocity_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    F_midpoint_y[ei][idx_midpoint_y] = F_midpoint_HLLC_y[ei][idx_midpoint_y];
                }
            }
        }
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* F_face_x = convective_flux->getPointer(0, ei);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    const int idx_midpoint_x_L = i +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    const int idx_midpoint_x_R = (i + 2) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                        (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                    
                    const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                        (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = Real(dt)*(
                        Real(1)/Real(30)*(F_midpoint_x[ei][idx_midpoint_x_R] +
                            F_midpoint_x[ei][idx_midpoint_x_L]) -
                        Real(3)/Real(10)*(F_node_x[ei][idx_node_R] +
                            F_node_x[ei][idx_node_L]) +
                        Real(23)/Real(15)*F_midpoint_x[ei][idx_midpoint_x]);
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* F_face_y = convective_flux->getPointer(1, ei);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const int idx_midpoint_y_B = (i + 1) +
                        j*(interior_dim_0 + 2);
                    
                    const int idx_midpoint_y_T = (i + 1) +
                        (j + 2)*(interior_dim_0 + 2);
                    
                    const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                        (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                    
                    const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                        (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                    
                    F_face_y[idx_face_y] = Real(dt)*(
                        Real(1)/Real(30)*(F_midpoint_y[ei][idx_midpoint_y_T] +
                            F_midpoint_y[ei][idx_midpoint_y_B]) -
                        Real(3)/Real(10)*(F_node_y[ei][idx_node_T] +
                            F_node_y[ei][idx_node_B]) +
                        Real(23)/Real(15)*F_midpoint_y[ei][idx_midpoint_y]);
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Compute the source.
         */
        
        t_compute_source->start();
        
        if (d_has_advective_eqn_form)
        {
            Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
            Real* v_midpoint_y = velocity_midpoint->getPointer(1, 1);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                    const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                    const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                    
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                            
                            const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                                (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                                (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                            
                            const int idx_cell_nghost = i + j*interior_dim_0;
                            
                            const int idx_midpoint_x_LL = i + 
                                (j + 1)*(interior_dim_0 + 3);
                            
                            const int idx_midpoint_x_L = (i + 1) +
                                (j + 1)*(interior_dim_0 + 3);
                            
                            const int idx_midpoint_x_R = (i + 2) +
                                (j + 1)*(interior_dim_0 + 3);
                            
                            const int idx_midpoint_x_RR = (i + 3) +
                                (j + 1)*(interior_dim_0 + 3);
                            
                            const int idx_midpoint_y_BB = (i + 1) +
                                j*(interior_dim_0 + 2);
                            
                            const int idx_midpoint_y_B = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2);
                            
                            const int idx_midpoint_y_T = (i + 1) +
                                (j + 2)*(interior_dim_0 + 2);
                            
                            const int idx_midpoint_y_TT = (i + 1) +
                                (j + 3)*(interior_dim_0 + 2);
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                (Real(3)/Real(2)*(u_midpoint_x[idx_midpoint_x_R] -
                                     u_midpoint_x[idx_midpoint_x_L]) -
                                 Real(3)/Real(10)*(u[idx_cell_wghost_x_R] -
                                     u[idx_cell_wghost_x_L]) +
                                 Real(1)/Real(30)*(u_midpoint_x[idx_midpoint_x_RR] -
                                     u_midpoint_x[idx_midpoint_x_LL]))/Real(dx[0]) +
                                (Real(3)/Real(2)*(v_midpoint_y[idx_midpoint_y_T] -
                                     v_midpoint_y[idx_midpoint_y_B]) -
                                 Real(3)/Real(10)*(v[idx_cell_wghost_y_T] -
                                     v[idx_cell_wghost_y_B]) +
                                 Real(1)/Real(30)*(v_midpoint_y[idx_midpoint_y_TT] -
                                     v_midpoint_y[idx_midpoint_y_BB]))/Real(dx[1]));
                        }
                    }
                }
            }
        }
        
        t_compute_source->stop();
        

    }
    else if (d_dim == tbox::Dimension(3))
    {
    }
    
    /*
     * Unregister the patch and data of all registered derived cell variables in the flow model.
     */
    
    d_flow_model->unregisterPatch();
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorCentral::performWENOInterpolation(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_minus,
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_plus,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& variables,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(variables_plus.size()) == d_num_eqn);
    
    TBOX_ASSERT(static_cast<int>(variables.size()) == 6);
#endif
    
    // Default scheme constants.
    const int constant_p = 2;
    const int constant_q = 4;
    const Real constant_C = Real(1.0e9);
    const Real constant_alpha_tau = Real(35);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = variables_minus[0]->getBox();
    const hier::IntVector interior_dims = variables_minus[0]->getBox().numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(variables_minus[ei]->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(variables_plus[ei]->getBox().numberCells() == interior_dims);
        
        TBOX_ASSERT(variables_minus[ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
        TBOX_ASSERT(variables_plus[ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
    }
    
    TBOX_ASSERT(static_cast<int>(variables.size()) == 6);
    
    for (int m = 0; m < 6; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(variables[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ASSERT(variables[m][ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
        }
    }
#endif
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    domain_lo = domain.lower() - interior_box.lower();
    domain_dims = domain.numberCells();
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index and the number of cells in each dimension.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                performLocalWENOInterpolationMinus(
                    U_L,
                    U_array.data(),
                    idx_midpoint_x,
                    constant_p,
                    constant_q,
                    constant_C,
                    constant_alpha_tau);
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_R = variables_plus[ei]->getPointer(0);
            
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                performLocalWENOInterpolationPlus(
                    U_R,
                    U_array.data(),
                    idx_midpoint_x,
                    constant_p,
                    constant_q,
                    constant_C,
                    constant_alpha_tau);
            }
        }
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices and the number of cells in each dimension.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        /*
         * Get the interior dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    performLocalWENOInterpolationMinus(
                        U_L,
                        U_array.data(),
                        idx_midpoint_x,
                        constant_p,
                        constant_q,
                        constant_C,
                        constant_alpha_tau);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_R = variables_plus[ei]->getPointer(0);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    performLocalWENOInterpolationPlus(
                        U_R,
                        U_array.data(),
                        idx_midpoint_x,
                        constant_p,
                        constant_q,
                        constant_C,
                        constant_alpha_tau);
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(1);
            
            for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    performLocalWENOInterpolationMinus(
                        U_B,
                        U_array.data(),
                        idx_midpoint_y,
                        constant_p,
                        constant_q,
                        constant_C,
                        constant_alpha_tau);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_T = variables_plus[ei]->getPointer(1);
            
            for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    performLocalWENOInterpolationPlus(
                        U_T,
                        U_array.data(),
                        idx_midpoint_y,
                        constant_p,
                        constant_q,
                        constant_C,
                        constant_alpha_tau);
                }
            }
        }
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices and the number of cells in each dimension.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        /*
         * Get the interior dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_L = variables_minus[ei]->getPointer(0);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolationMinus(
                            U_L,
                            U_array.data(),
                            idx_midpoint_x,
                            constant_p,
                            constant_q,
                            constant_C,
                            constant_alpha_tau);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            Real* U_R = variables_plus[ei]->getPointer(0);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0 - 1; i < domain_lo_0 + domain_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolationPlus(
                            U_R,
                            U_array.data(),
                            idx_midpoint_x,
                            constant_p,
                            constant_q,
                            constant_C,
                            constant_alpha_tau);
                    }
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(1);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 2; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        performLocalWENOInterpolationMinus(
                            U_B,
                            U_array.data(),
                            idx_midpoint_y,
                            constant_p,
                            constant_q,
                            constant_C,
                            constant_alpha_tau);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            Real* U_T = variables_plus[ei]->getPointer(1);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1 - 1; j < domain_lo_1 + domain_dim_1 + 2; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        performLocalWENOInterpolationPlus(
                            U_T,
                            U_array.data(),
                            idx_midpoint_y,
                            constant_p,
                            constant_q,
                            constant_C,
                            constant_alpha_tau);
                    }
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(2));
            }
            
            Real* U_B = variables_minus[ei]->getPointer(2);
            
            for (int k = domain_lo_2 - 1; k < domain_lo_2 + domain_dim_2 + 2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolationMinus(
                            U_B,
                            U_array.data(),
                            idx_midpoint_z,
                            constant_p,
                            constant_q,
                            constant_C,
                            constant_alpha_tau);
                    }
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<Real*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(2));
            }
            
            Real* U_F = variables_plus[ei]->getPointer(2);
            
            for (int k = domain_lo_2 - 1; k < domain_lo_2 + domain_dim_2 + 2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolationPlus(
                            U_F,
                            U_array.data(),
                            idx_midpoint_z,
                            constant_p,
                            constant_q,
                            constant_C,
                            constant_alpha_tau);
                    }
                }
            }
        }
        
    } // if (d_dim == tbox::Dimension(3))
}