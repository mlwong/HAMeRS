#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
NavierStokesInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const double data_time,
    const bool initial_time)
{
    // Follow Reckinger, Scott J., Daniel Livescu, and Oleg V. Vasilyev.
    // "Comprehensive numerical methodology for direct numerical simulations of compressible Rayleigh–Taylor instability."
    // Journal of Computational Physics 313 (2016): 181-208.
    // Note that the sign of gravity in the paper is flipped.
    NULL_USE(data_time);
    
    if ((d_project_name != "2D discontinuous Rayleigh-Taylor instability") &&
        (d_project_name != "2D smooth Rayleigh-Taylor instability") &&
        (d_project_name != "2D smooth multi-mode Rayleigh-Taylor instability")
       ) 
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D discontinuous Rayleigh-Taylor instability' or "
            << "'2D smooth Rayleigh-Taylor instability' or "
            << "'2D smooth multi-mode Rayleigh-Taylor instability'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if (d_dim != tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Dimension of problem should be 2!"
            << std::endl);
    }
    
    if (d_flow_model_type != FLOW_MODEL::FOUR_EQN_CONSERVATIVE)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be conservative four-equation models!"
            << std::endl);
    }
    
    if (d_flow_model->getNumberOfSpecies() != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species should be 2!"
            << std::endl);
    }
    
    if (initial_time)
    {
        const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
            HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif
        
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        // Get the dimensions of box that covers the interior of Patch.
        hier::Box patch_box = patch.getBox();
        const hier::IntVector patch_dims = patch_box.numberCells();
        
        /*
         * Initialize data for a 2D Rayleigh-Taylor instability problem (At = 0.04, M = 0.3).
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<double> > partial_density = conservative_variables[0];
        HAMERS_SHARED_PTR<pdat::CellData<double> > momentum        = conservative_variables[1];
        HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy    = conservative_variables[2];
        
        double* rho_Y_0 = partial_density->getPointer(0);
        double* rho_Y_1 = partial_density->getPointer(1);
        double* rho_u   = momentum->getPointer(0);
        double* rho_v   = momentum->getPointer(1);
        double* E       = total_energy->getPointer(0);
        
        const double gamma = double(7)/double(5); // assume both gases have the same ratio of specific heat ratios
        // const double gamma_0 = double(7)/double(5);
        // const double gamma_1 = double(7)/double(5);
        
              double lambda = 701.53278340668; // wavelength of single-mode perturbation
              double eta_0  = 0.01*lambda;      // 1% perturbation
        // const double eta_0  = 0.0*lambda;      // no perturbation
        
        const double W_1 = 0.03328; // molecular weight of heavier gas
        const double W_2 = 0.03072; // molecular weight of lighter gas
        
        const double p_i = 100000.0; // interface pressure
        const double T_0 = 300.0;    // background temperature
        
        TBOX_ASSERT(d_initial_conditions_db != nullptr);
        TBOX_ASSERT(d_initial_conditions_db->keyExists("gravity"));
        
        std::vector<double> gravity_vector = d_initial_conditions_db->getDoubleVector("gravity");
        const double g = gravity_vector[0]; // gravity
        
        const double R_u = 8.31446261815324; // universal gas constant
        const double R_1 = R_u/W_1;          // gas constant of heavier gas
        const double R_2 = R_u/W_2;          // gas constant of lighter gas
        
        // const double rho_i = p_i/(R_u*T_0)*(W_1 + W_2)/2.0;
        
        if (d_project_name == "2D discontinuous Rayleigh-Taylor instability")
        {
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
                    
                    const double eta = eta_0*cos(2.0*M_PI/lambda*x[1]);
                    
                    if (x[0] < eta) // heavier fluid
                    {
                        const double rho = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                        rho_Y_0[idx_cell] = rho;
                        rho_Y_1[idx_cell] = 0.0;
                        
                        const double p = p_i*exp((g*x[0])/(R_1*T_0));
                        
                        const double u = 0.0;
                        const double v = 0.0;
                        
                        rho_u[idx_cell] = rho*u;
                        rho_v[idx_cell] = rho*v;
                        E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                    }
                    else // lighter fluid
                    {
                        const double rho = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                        rho_Y_0[idx_cell] = 0.0;
                        rho_Y_1[idx_cell] = rho;
                        
                        const double p = p_i*exp((g*x[0])/(R_2*T_0));
                        
                        const double u = 0.0;
                        const double v = 0.0;
                        
                        rho_u[idx_cell] = rho*u;
                        rho_v[idx_cell] = rho*v;
                        E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                    }
                }
            }
        }
        else if (d_project_name == "2D smooth Rayleigh-Taylor instability")
        {
            const double delta = 0.01*lambda; // characteristic length of interface.
            
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
                    
                    const double eta = eta_0*cos(2.0*M_PI/lambda*x[1]);
                    
                    const double X_2_H = 0.5*(1.0 + erf((x[0] - eta)/delta)); // mass fraction of second species (Y_2)
                    const double R_H   = R_1*(1.0 - X_2_H) + X_2_H*R_2;
                    
                    const int N_int = 10000; // number of numerical quadrature points
                    const double dx_p = x[0]/(N_int - 1.0);
                    
                    double integral = 0.0;
                    for (int ii = 0; ii < N_int; ii++)
                    {
                        const double x_p = x[0] + ii*dx_p;
                        integral += 1.0/(0.5*(R_2 - R_1)*erf((x_p - eta)/delta) + 0.5*(R_1 + R_2))*dx_p;
                    }
                    
                    const double p_H = p_i*exp(g/T_0*integral);
                    const double rho_H = p_H/(R_H*T_0);
                    
                    // Scott's implementation
                    // const double dX_2_H_dx = 1.0/(delta*sqrt(M_PI))*exp(-(x[0]/delta)*(x[0]/delta));
                    // const double dlnR_H_dx = (R_2 - R_1)*dX_2_H_dx;
                    // const double p_H = p_i*exp(g/(R_H*T_0)*(x[0] - 0.5*delta*delta*dlnR_H_dx));
                    // const double rho_H = p_H/(R_H*T_0);
                    
                    // const double X_2 = 0.5*(1.0 + erf((x[0] - eta)/delta)); // mass fraction of second species (Y_2)
                    
                    double rho, p;
                    
                    rho = rho_H;
                    p   = p_H;
                    
                    // if (x[0] < eta)
                    // {
                    //     const double p_1_H   = p_i*exp((g*x[0])/(R_1*T_0));
                    //     const double rho_1_H = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                    //     
                    //     const double p_1   = p_i*exp((g*x[0])/(R_1*T_0));
                    //     const double rho_1 = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                    //     
                    //     p   = p_1*p_H/p_1_H;
                    //     rho = rho_1*rho_H/rho_1_H;
                    //     
                    //     
                    // }
                    // else
                    // {
                    //     const double p_2_H   = p_i*exp((g*x[0])/(R_2*T_0));
                    //     const double rho_2_H = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                    //     
                    //     const double p_2   = p_i*exp((g*x[0])/(R_2*T_0));
                    //     const double rho_2 = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                    //     
                    //     p   = p_2*p_H/p_2_H;
                    //     rho = rho_2*rho_H/rho_2_H;
                    // }
                    
                    rho_Y_0[idx_cell] = rho*(1.0 - X_2_H);
                    rho_Y_1[idx_cell] = rho*X_2_H;
                    
                    const double u = 0.0;
                    const double v = 0.0;
                    
                    rho_u[idx_cell] = rho*u;
                    rho_v[idx_cell] = rho*v;
                    E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
                }
            }
        }
        else if (d_project_name == "2D smooth multi-mode Rayleigh-Taylor instability")
        {
            lambda  = lambda/4.0;
            eta_0   = lambda*0.04;
            const double delta   = 0.04*lambda; // characteristic length of interface
            const int    waven   = 16;          // dominant wave number
            const double width   = 16.0*lambda; // domain size in y direction
            
            // Seed for "random" phase shifts.
            int random_seed = 0;
            if (d_initial_conditions_db->keyExists("random_seed"))
            {
                random_seed = d_initial_conditions_db->getInteger("random_seed");
                if (random_seed < 0 || random_seed > 15)
                {
                    TBOX_ERROR(d_object_name << ": "
                        << "'random_seed' should be in between 0 - 15."
                        << std::endl);
                }
            }
            else
            {
                TBOX_WARNING(d_object_name << ": "
                    << "'random_seed' not given for '2D smooth multi-mode Rayleigh-Taylor instability'!"
                    << "  'random_seed' = 0 is assumed."
                    << std::endl);
            }
            
            double rmod[9]; // random seed
            switch (random_seed)
            {
                case 0:
                {
                    rmod[0] = 6.031966614958411e+000;
                    rmod[1] = 1.273017034173460e+000;
                    rmod[2] = 5.934447177754063e+000;
                    rmod[3] = 3.101658133166612e+000;
                    rmod[4] = 2.294026034817427e+000;
                    rmod[5] = 4.916046917518752e+000;
                    rmod[6] = 0.571212135466553e+000;
                    rmod[7] = 4.966766749458944e+000;
                    rmod[8] = 5.027899324302027e+000;
                    
                    break;
                }
                case 1:
                {
                    rmod[0] = 2.620226532717789200e+000;
                    rmod[1] = 4.525932273597345700e+000;
                    rmod[2] = 7.186381718527406600e-004;
                    rmod[3] = 1.899611578242180700e+000;
                    rmod[4] = 9.220944569241362700e-001;
                    rmod[5] = 5.801805019369201700e-001;
                    rmod[6] = 1.170307423440345900e+000;
                    rmod[7] = 2.171222082895173200e+000;
                    rmod[8] = 2.492963564452900500e+000;
                    
                    break;
                }
                case 2:
                {
                    rmod[0] = 2.739436763143839700e+00;
                    rmod[1] = 1.628993188915385800e-01;
                    rmod[2] = 3.453631204915429600e+00;
                    rmod[3] = 2.735211261185420500e+00;
                    rmod[4] = 2.641248797687487200e+00;
                    rmod[5] = 2.075554893781340400e+00;
                    rmod[6] = 1.285845290520944100e+00;
                    rmod[7] = 3.890994236937394200e+00;
                    rmod[8] = 1.882785842859457300e+00;
                    
                    break;
                }
                case 3:
                {
                    rmod[0] = 3.460765288681905800e+00;
                    rmod[1] = 4.449423994385291800e+00;
                    rmod[2] = 1.827808381326725400e+00;
                    rmod[3] = 3.209624503479690600e+00;
                    rmod[4] = 5.610551183647944900e+00;
                    rmod[5] = 5.631575567313184600e+00;
                    rmod[6] = 7.890757775039627400e-01;
                    rmod[7] = 1.302145406935464500e+00;
                    rmod[8] = 3.233779755813989700e-01;
                    
                    break;
                }
                case 4:
                {
                    rmod[0] = 6.076027676094973600e+00;
                    rmod[1] = 3.438361627635736700e+00;
                    rmod[2] = 6.111556079054740700e+00;
                    rmod[3] = 4.491321348791744100e+00;
                    rmod[4] = 4.383959499105254800e+00;
                    rmod[5] = 1.357730343666468900e+00;
                    rmod[6] = 6.134113310024843300e+00;
                    rmod[7] = 3.914584796145817400e-02;
                    rmod[8] = 1.589535062303236700e+00;
                    
                    break;
                }
                case 5:
                {
                    rmod[0] = 1.394824230885255200e+00;
                    rmod[1] = 5.470972432660288700e+00;
                    rmod[2] = 1.298854759541258500e+00;
                    rmod[3] = 5.771802559770447900e+00;
                    rmod[4] = 3.068778005297785300e+00;
                    rmod[5] = 3.843700051147186600e+00;
                    rmod[6] = 4.812340990490530300e+00;
                    rmod[7] = 3.257316284380881800e+00;
                    rmod[8] = 1.864852550667249300e+00;
                    
                    break;
                }
                case 6:
                {
                    rmod[0] = 5.610005784868826100e+00;
                    rmod[1] = 2.085890634948696300e+00;
                    rmod[2] = 5.159934759824945000e+00;
                    rmod[3] = 2.619876261158565800e-01;
                    rmod[4] = 6.764268695934091400e-01;
                    rmod[5] = 3.738822386827532500e+00;
                    rmod[6] = 3.328940665618581800e+00;
                    rmod[7] = 2.631444681644834500e+00;
                    rmod[8] = 2.107429670466887600e+00;
                    
                    break;
                }
                case 7:
                {
                    rmod[0] = 4.794591226104558700e-01;
                    rmod[1] = 4.900374296196336100e+00;
                    rmod[2] = 2.754606441521316700e+00;
                    rmod[3] = 4.545665775603436200e+00;
                    rmod[4] = 6.144889332352788000e+00;
                    rmod[5] = 3.383469340939719400e+00;
                    rmod[6] = 3.148632734395143500e+00;
                    rmod[7] = 4.527106224916906400e-01;
                    rmod[8] = 1.686651855650350300e+00;
                    
                    break;
                }
                case 8:
                {
                    rmod[0] = 5.487918790480180500e+00;
                    rmod[1] = 6.085520462042457400e+00;
                    rmod[2] = 5.461310364152818200e+00;
                    rmod[3] = 3.335464681414797900e+00;
                    rmod[4] = 1.462275210911384800e+00;
                    rmod[5] = 7.162079955536335100e-02;
                    rmod[6] = 2.704715354294335400e+00;
                    rmod[7] = 2.528048153710494200e+00;
                    rmod[8] = 3.284061815477384600e+00;
                    
                    break;
                }
                case 9:
                {
                    rmod[0] = 6.518273126904997100e-02;
                    rmod[1] = 3.153371063435702800e+00;
                    rmod[2] = 3.115035471112504800e+00;
                    rmod[3] = 8.408757300236915400e-01;
                    rmod[4] = 8.929102841152991600e-01;
                    rmod[5] = 1.373244659450403700e+00;
                    rmod[6] = 2.629564450717737100e+00;
                    rmod[7] = 1.558865616070162600e+00;
                    rmod[8] = 5.281623651287690200e-01;
                    
                    break;
                }
                case 10:
                {
                    rmod[0] = 4.846350532897925100e+00;
                    rmod[1] = 1.303883433103263400e-01;
                    rmod[2] = 3.981329279609052500e+00;
                    rmod[3] = 4.704873552725635100e+00;
                    rmod[4] = 3.132211935225629200e+00;
                    rmod[5] = 1.412438980302679600e+00;
                    rmod[6] = 1.244465681755566800e+00;
                    rmod[7] = 4.778555396547323800e+00;
                    rmod[8] = 1.062554723574571100e+00;
                    
                    break;
                }
                case 11:
                {
                    rmod[0] = 1.132667860480351500e+00;
                    rmod[1] = 1.223665511688170900e-01;
                    rmod[2] = 2.910487839707776500e+00;
                    rmod[3] = 4.554894212576069600e+00;
                    rmod[4] = 2.640217114369509700e+00;
                    rmod[5] = 3.050028410914633200e+00;
                    rmod[6] = 8.030422644949872200e-02;
                    rmod[7] = 3.062246122248716100e+00;
                    rmod[8] = 5.917545720207830800e+00;
                    
                    break;
                }
                case 12:
                {
                    rmod[0] = 9.686337061529999300e-01;
                    rmod[1] = 4.649869379728302800e+00;
                    rmod[2] = 1.654457034571007900e+00;
                    rmod[3] = 3.353583514350031900e+00;
                    rmod[4] = 9.157719014108255100e-02;
                    rmod[5] = 5.772657702308401400e+00;
                    rmod[6] = 5.659358337346415800e+00;
                    rmod[7] = 2.099930230068142400e-01;
                    rmod[8] = 6.012690009399070900e+00;
                    
                    break;
                }
                case 13:
                {
                    rmod[0] = 4.886448359475572500e+00;
                    rmod[1] = 1.492515503572874100e+00;
                    rmod[2] = 5.179094765441459600e+00;
                    rmod[3] = 6.067981171564244200e+00;
                    rmod[4] = 6.111033028633724700e+00;
                    rmod[5] = 2.849105648924096900e+00;
                    rmod[6] = 3.826726653470131600e+00;
                    rmod[7] = 4.872776801893367700e+00;
                    rmod[8] = 4.031375540680533800e+00;
                    
                    break;
                }
                case 14:
                {
                    rmod[0] = 3.229201266315314500e+00;
                    rmod[1] = 4.857939295249377000e+00;
                    rmod[2] = 5.469058445908473200e+00;
                    rmod[3] = 5.056046877018226200e-02;
                    rmod[4] = 1.946128216239969300e+00;
                    rmod[5] = 6.016801745283995500e+00;
                    rmod[6] = 3.224007387389883600e+00;
                    rmod[7] = 1.999840021860475000e+00;
                    rmod[8] = 3.387893124464984600e+00;
                    
                    break;
                }
                case 15:
                {
                    rmod[0] = 5.333278883951943600e+00;
                    rmod[1] = 1.124036246977920200e+00;
                    rmod[2] = 3.415741493812250500e-01;
                    rmod[3] = 2.271613052442061200e+00;
                    rmod[4] = 1.730395068203421900e+00;
                    rmod[5] = 3.330089625864812500e+00;
                    rmod[6] = 1.922145236540217400e+00;
                    rmod[7] = 1.913068819850232400e+00;
                    rmod[8] = 7.020911446719929600e-01;
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name << ": "
                        << "'random_seed' should be in between 0 - 15."
                        << std::endl);
                }
            }
            
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
                    
                    double eta = 0.0;        
                    for (int m = waven - 4; m <= waven + 4; m++)
                    {
                        eta += eta_0/3.0*cos(2.0*M_PI*m/width*x[1] + rmod[m-waven+4]);
                    }
                    
                    const double X_2_H = 0.5*(1.0 + erf((x[0] - eta)/delta)); // mass fraction of second species (Y_2)
                    const double R_H   = R_1*(1.0 - X_2_H) + X_2_H*R_2;
                    
                    const int N_int = 10000; // number of numerical quadrature points
                    const double dx_p = x[0]/(N_int - 1.0);
                    
                    double integral = 0.0;
                    for (int ii = 0; ii < N_int; ii++)
                    {
                        const double x_p = x[0] + ii*dx_p;
                        integral += 1.0/(0.5*(R_2 - R_1)*erf((x_p - eta)/delta) + 0.5*(R_1 + R_2))*dx_p;
                    }
                    
                    const double p_H = p_i*exp(g/T_0*integral);
                    const double rho_H = p_H/(R_H*T_0);
                    
                    double rho, p;
                    rho = rho_H;
                    p   = p_H;
                    
                    // if (x[0] < eta)
                    // {
                    //     const double p_1_H   = p_i*exp((g*x[0])/(R_1*T_0));
                    //     const double rho_1_H = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                    //
                    //     const double p_1   = p_i*exp((g*x[0])/(R_1*T_0));
                    //     const double rho_1 = p_i/(R_1*T_0)*exp((g*x[0])/(R_1*T_0));
                    //
                    //     p   = p_1*p_H/p_1_H;
                    //     rho = rho_1*rho_H/rho_1_H;
                    //
                    //
                    // }
                    // else
                    // {
                    //     const double p_2_H   = p_i*exp((g*x[0])/(R_2*T_0));
                    //     const double rho_2_H = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                    //
                    //     const double p_2   = p_i*exp((g*x[0])/(R_2*T_0));
                    //     const double rho_2 = p_i/(R_2*T_0)*exp((g*x[0])/(R_2*T_0));
                    //
                    //     p   = p_2*p_H/p_2_H;
                    //     rho = rho_2*rho_H/rho_2_H;
                    // }
                    
                    rho_Y_0[idx_cell] = rho*(1.0 - X_2_H);
                    rho_Y_1[idx_cell] = rho*X_2_H;
                    
                    const double u = 0.0;
                    const double v = 0.0;
                    
                    rho_u[idx_cell] = rho*u;
                    rho_v[idx_cell] = rho*v;
                    E[idx_cell]     = p/(gamma - double(1)) + double(1)/double(2)*rho*(u*u + v*v);
               }
            }
        }
    }
}
