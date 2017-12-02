#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"

#include <sstream>

/*
 * Set the data on the patch interior to some initial values.
 */
void
NavierStokesInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& conservative_variables,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(data_time);
    
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
            << "Flow model should be conservative four-equation model!"
            << std::endl);
    }
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of species should be 2!"
            << std::endl);
    }
    
    if (initial_time)
    {
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
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
         * Initialize data for 2D binary mass diffusion problem.
         */
        
        boost::shared_ptr<pdat::CellVariable<double> > var_partial_density = conservative_variables[0];
        boost::shared_ptr<pdat::CellVariable<double> > var_momentum        = conservative_variables[1];
        boost::shared_ptr<pdat::CellVariable<double> > var_total_energy    = conservative_variables[2];
        
        boost::shared_ptr<pdat::CellData<double> > partial_density(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(var_partial_density, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > momentum(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(var_momentum, data_context)));
        
        boost::shared_ptr<pdat::CellData<double> > total_energy(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                patch.getPatchData(var_total_energy, data_context)));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(partial_density);
        TBOX_ASSERT(momentum);
        TBOX_ASSERT(total_energy);
#endif
        
        double* rho_Y_1 = partial_density->getPointer(0);
        double* rho_Y_2 = partial_density->getPointer(1);
        double* rho_u   = momentum->getPointer(0);
        double* rho_v   = momentum->getPointer(1);
        double* E       = total_energy->getPointer(0);
        
        if (d_project_name.find("2D Poggi's RMI 1mm smooth interface reduced domain size") != std::string::npos)
        {
            /*
             * Get the settings.
             */
            
            std::string settings = d_project_name.substr(56);
            
            std::stringstream ss(settings);
            
            std::string A_idx_str;
            std::string m_idx_str;
            std::string random_seed_str;
            
            std::getline(ss, A_idx_str, '-');
            std::getline(ss, m_idx_str, '-');
            std::getline(ss, random_seed_str);
            
            double A_candidates[] = {0.01e-3, sqrt(2.0)*0.01e-3, 0.02e-3};
            int m_min_candidates[] = {25, 20, 15, 10};
            int m_max_candidates[] = {30, 30, 30, 30};
            
            // Seed for "random" phase shifts.
            int random_seed = std::stoi(random_seed_str) - 1;
            
            double phase_shifts[21];
            if (random_seed == 0)
            {
                phase_shifts[0] = 5.119059895756811000e+000;
                phase_shifts[1] = 5.691258590395267300e+000;
                phase_shifts[2] = 7.978816983408705300e-001;
                phase_shifts[3] = 5.738909759225261800e+000;
                phase_shifts[4] = 3.973230324742651500e+000;
                phase_shifts[5] = 6.128644395486362300e-001;
                phase_shifts[6] = 1.749855916861123200e+000;
                phase_shifts[7] = 3.436157926236805200e+000;
                phase_shifts[8] = 6.016192879924800800e+000;
                phase_shifts[9] = 6.062573467430127000e+000;
                phase_shifts[10] = 9.903121990156674700e-001;
                phase_shifts[11] = 6.098414305612863000e+000;
                phase_shifts[12] = 6.014057305717998700e+000;
                phase_shifts[13] = 3.049705144518116000e+000;
                phase_shifts[14] = 5.028310483744898600e+000;
                phase_shifts[15] = 8.914981581520268200e-001;
                phase_shifts[16] = 2.650004294134627800e+000;
                phase_shifts[17] = 5.753735997130328400e+000;
                phase_shifts[18] = 4.977585453328568800e+000;
                phase_shifts[19] = 6.028668715861979200e+000;
                phase_shifts[20] = 4.120140326260335300e+000;
            }
            else if (random_seed == 1)
            {
                phase_shifts[0] = 2.620226532717789200e+000;
                phase_shifts[1] = 4.525932273597345700e+000;
                phase_shifts[2] = 7.186381718527406600e-004;
                phase_shifts[3] = 1.899611578242180700e+000;
                phase_shifts[4] = 9.220944569241362700e-001;
                phase_shifts[5] = 5.801805019369201700e-001;
                phase_shifts[6] = 1.170307423440345900e+000;
                phase_shifts[7] = 2.171222082895173200e+000;
                phase_shifts[8] = 2.492963564452900500e+000;
                phase_shifts[9] = 3.385485386352383500e+000;
                phase_shifts[10] = 2.633876813749063600e+000;
                phase_shifts[11] = 4.305361097085856200e+000;
                phase_shifts[12] = 1.284611371532881700e+000;
                phase_shifts[13] = 5.517374574309792800e+000;
                phase_shifts[14] = 1.720813231802212400e-001;
                phase_shifts[15] = 4.212671608894216200e+000;
                phase_shifts[16] = 2.622003402848613000e+000;
                phase_shifts[17] = 3.510351721361030500e+000;
                phase_shifts[18] = 8.820771499014955500e-001;
                phase_shifts[19] = 1.244708365548507600e+000;
                phase_shifts[20] = 5.031226508705986900e+000;
            }
            else if (random_seed == 2)
            {
                phase_shifts[0] =  2.739436763143839700e+00;
                phase_shifts[1] =  1.628993188915385800e-01;
                phase_shifts[2] =  3.453631204915429600e+00;
                phase_shifts[3] =  2.735211261185420500e+00;
                phase_shifts[4] =  2.641248797687487200e+00;
                phase_shifts[5] =  2.075554893781340400e+00;
                phase_shifts[6] =  1.285845290520944100e+00;
                phase_shifts[7] =  3.890994236937394200e+00;
                phase_shifts[8] =  1.882785842859457300e+00;
                phase_shifts[9] =  1.676525214481097100e+00;
                phase_shifts[10] =  3.902698971848176200e+00;
                phase_shifts[11] =  3.324697832171727100e+00;
                phase_shifts[12] =  8.455907352323183100e-01;
                phase_shifts[13] =  3.226906505625833700e+00;
                phase_shifts[14] =  1.158869853890875000e+00;
                phase_shifts[15] =  4.934406261973431500e+00;
                phase_shifts[16] =  5.365685011406823100e+00;
                phase_shifts[17] =  3.105381634905035200e+00;
                phase_shifts[18] =  5.319102686422169800e+00;
                phase_shifts[19] =  5.004272909266416200e-01;
                phase_shifts[20] =  3.174554809962623300e+00;
            }
            else if (random_seed == 3)
            {
                phase_shifts[0] =  3.460765288681905800e+00;
                phase_shifts[1] =  4.449423994385291800e+00;
                phase_shifts[2] =  1.827808381326725400e+00;
                phase_shifts[3] =  3.209624503479690600e+00;
                phase_shifts[4] =  5.610551183647944900e+00;
                phase_shifts[5] =  5.631575567313184600e+00;
                phase_shifts[6] =  7.890757775039627400e-01;
                phase_shifts[7] =  1.302145406935464500e+00;
                phase_shifts[8] =  3.233779755813989700e-01;
                phase_shifts[9] =  2.769689932885809600e+00;
                phase_shifts[10] =  1.877177692264108100e-01;
                phase_shifts[11] =  2.870367803348323800e+00;
                phase_shifts[12] =  4.078692342216150700e+00;
                phase_shifts[13] =  1.749787202570126900e+00;
                phase_shifts[14] =  4.249034864029732200e+00;
                phase_shifts[15] =  3.712500572949150300e+00;
                phase_shifts[16] =  1.506826109907144800e-01;
                phase_shifts[17] =  3.511383794521557400e+00;
                phase_shifts[18] =  1.628931165259342300e+00;
                phase_shifts[19] =  2.608157742046328200e+00;
                phase_shifts[20] =  1.781440628002440600e+00;
            }
            else if (random_seed == 4)
            {
                phase_shifts[0] =  6.076027676094973600e+00;
                phase_shifts[1] =  3.438361627635736700e+00;
                phase_shifts[2] =  6.111556079054740700e+00;
                phase_shifts[3] =  4.491321348791744100e+00;
                phase_shifts[4] =  4.383959499105254800e+00;
                phase_shifts[5] =  1.357730343666468900e+00;
                phase_shifts[6] =  6.134113310024843300e+00;
                phase_shifts[7] =  3.914584796145817400e-02;
                phase_shifts[8] =  1.589535062303236700e+00;
                phase_shifts[9] =  2.731875768089710600e+00;
                phase_shifts[10] =  4.897007322881202100e+00;
                phase_shifts[11] =  1.242091956177010500e+00;
                phase_shifts[12] =  5.422346418112404400e+00;
                phase_shifts[13] =  6.178888685898380500e+00;
                phase_shifts[14] =  1.029451163889373700e+00;
                phase_shifts[15] =  3.753159859998575600e+00;
                phase_shifts[16] =  5.646131683366221300e-02;
                phase_shifts[17] =  2.428899003284019200e+00;
                phase_shifts[18] =  2.774658271593974900e-01;
                phase_shifts[19] =  6.010827870811834100e+00;
                phase_shifts[20] =  2.740390203450787900e+00;
            }
            else if (random_seed == 5)
            {
                phase_shifts[0] =  1.394824230885255200e+00;
                phase_shifts[1] =  5.470972432660288700e+00;
                phase_shifts[2] =  1.298854759541258500e+00;
                phase_shifts[3] =  5.771802559770447900e+00;
                phase_shifts[4] =  3.068778005297785300e+00;
                phase_shifts[5] =  3.843700051147186600e+00;
                phase_shifts[6] =  4.812340990490530300e+00;
                phase_shifts[7] =  3.257316284380881800e+00;
                phase_shifts[8] =  1.864852550667249300e+00;
                phase_shifts[9] =  1.179487265770075700e+00;
                phase_shifts[10] =  5.073123535864998400e-01;
                phase_shifts[11] =  4.639757219306710000e+00;
                phase_shifts[12] =  2.772827625222693500e+00;
                phase_shifts[13] =  9.946902347936740200e-01;
                phase_shifts[14] =  5.528807425687100300e+00;
                phase_shifts[15] =  1.722136030886382000e+00;
                phase_shifts[16] =  2.602715385609317300e+00;
                phase_shifts[17] =  1.860325083102777100e+00;
                phase_shifts[18] =  3.950790950403745900e+00;
                phase_shifts[19] =  3.643228409530135300e+00;
                phase_shifts[20] =  3.769466313582174900e+00;
            }
            else if (random_seed == 6)
            {
                phase_shifts[0] =  5.610005784868826100e+00;
                phase_shifts[1] =  2.085890634948696300e+00;
                phase_shifts[2] =  5.159934759824945000e+00;
                phase_shifts[3] =  2.619876261158565800e-01;
                phase_shifts[4] =  6.764268695934091400e-01;
                phase_shifts[5] =  3.738822386827532500e+00;
                phase_shifts[6] =  3.328940665618581800e+00;
                phase_shifts[7] =  2.631444681644834500e+00;
                phase_shifts[8] =  2.107429670466887600e+00;
                phase_shifts[9] =  3.911404949808252600e+00;
                phase_shifts[10] =  2.752923770993987800e+00;
                phase_shifts[11] =  4.623683638426456400e+00;
                phase_shifts[12] =  3.254918772457836300e+00;
                phase_shifts[13] =  3.637075851875532200e+00;
                phase_shifts[14] =  4.054885656444668000e+00;
                phase_shifts[15] =  6.221762592706061100e+00;
                phase_shifts[16] =  5.151320976870061400e+00;
                phase_shifts[17] =  2.596218042549897300e+00;
                phase_shifts[18] =  5.505752056350989000e+00;
                phase_shifts[19] =  5.175833163377692600e+00;
                phase_shifts[20] =  3.422734271652580500e-01;
            }
            else if (random_seed == 7)
            {
                phase_shifts[0] =  4.794591226104558700e-01;
                phase_shifts[1] =  4.900374296196336100e+00;
                phase_shifts[2] =  2.754606441521316700e+00;
                phase_shifts[3] =  4.545665775603436200e+00;
                phase_shifts[4] =  6.144889332352788000e+00;
                phase_shifts[5] =  3.383469340939719400e+00;
                phase_shifts[6] =  3.148632734395143500e+00;
                phase_shifts[7] =  4.527106224916906400e-01;
                phase_shifts[8] =  1.686651855650350300e+00;
                phase_shifts[9] =  3.140854384503345600e+00;
                phase_shifts[10] =  4.267727931822740600e+00;
                phase_shifts[11] =  5.050041302457694700e+00;
                phase_shifts[12] =  2.393523730699239000e+00;
                phase_shifts[13] =  4.142902860882791700e-01;
                phase_shifts[14] =  1.810472195900441500e+00;
                phase_shifts[15] =  5.715144688873524000e+00;
                phase_shifts[16] =  1.340739718380646000e+00;
                phase_shifts[17] =  2.840778633916690900e+00;
                phase_shifts[18] =  5.850939980867245500e+00;
                phase_shifts[19] =  1.564464607044678000e-01;
                phase_shifts[20] =  3.773360134453180900e+00;
            }
            else if (random_seed == 8)
            {
                phase_shifts[0] =  5.487918790480180500e+00;
                phase_shifts[1] =  6.085520462042457400e+00;
                phase_shifts[2] =  5.461310364152818200e+00;
                phase_shifts[3] =  3.335464681414797900e+00;
                phase_shifts[4] =  1.462275210911384800e+00;
                phase_shifts[5] =  7.162079955536335100e-02;
                phase_shifts[6] =  2.704715354294335400e+00;
                phase_shifts[7] =  2.528048153710494200e+00;
                phase_shifts[8] =  3.284061815477384600e+00;
                phase_shifts[9] =  3.005824302460911900e+00;
                phase_shifts[10] =  3.489407636567531100e+00;
                phase_shifts[11] =  3.414195041550118300e+00;
                phase_shifts[12] =  4.780847900926277200e+00;
                phase_shifts[13] =  4.475981457080911800e+00;
                phase_shifts[14] =  3.893577441359380200e+00;
                phase_shifts[15] =  2.677213551766147500e+00;
                phase_shifts[16] =  1.816311968715775800e+00;
                phase_shifts[17] =  6.118912942830475800e+00;
                phase_shifts[18] =  2.097164178670534600e+00;
                phase_shifts[19] =  1.374767610538450100e+00;
                phase_shifts[20] =  4.134862952337833900e-01;
            }
            else if (random_seed == 9)
            {
                phase_shifts[0] =  6.518273126904997100e-02;
                phase_shifts[1] =  3.153371063435702800e+00;
                phase_shifts[2] =  3.115035471112504800e+00;
                phase_shifts[3] =  8.408757300236915400e-01;
                phase_shifts[4] =  8.929102841152991600e-01;
                phase_shifts[5] =  1.373244659450403700e+00;
                phase_shifts[6] =  2.629564450717737100e+00;
                phase_shifts[7] =  1.558865616070162600e+00;
                phase_shifts[8] =  5.281623651287690200e-01;
                phase_shifts[9] =  2.170831978815202400e+00;
                phase_shifts[10] =  1.047886690116973200e+00;
                phase_shifts[11] =  5.520149537707936800e+00;
                phase_shifts[12] =  5.975083231173248200e+00;
                phase_shifts[13] =  2.434632257067045100e-01;
                phase_shifts[14] =  4.392621289777343600e+00;
                phase_shifts[15] =  3.598756057728601700e+00;
                phase_shifts[16] =  5.642345130974941700e+00;
                phase_shifts[17] =  4.190249828021817900e+00;
                phase_shifts[18] =  3.442166309155490800e+00;
                phase_shifts[19] =  4.413481670264409300e+00;
                phase_shifts[20] =  2.428307673683218300e+00;
            }
            else if (random_seed == 10)
            {
                phase_shifts[0] =  4.846350532897925100e+00;
                phase_shifts[1] =  1.303883433103263400e-01;
                phase_shifts[2] =  3.981329279609052500e+00;
                phase_shifts[3] =  4.704873552725635100e+00;
                phase_shifts[4] =  3.132211935225629200e+00;
                phase_shifts[5] =  1.412438980302679600e+00;
                phase_shifts[6] =  1.244465681755566800e+00;
                phase_shifts[7] =  4.778555396547323800e+00;
                phase_shifts[8] =  1.062554723574571100e+00;
                phase_shifts[9] =  5.550554224571162500e-01;
                phase_shifts[10] =  4.306242740899813600e+00;
                phase_shifts[11] =  5.990347064774805800e+00;
                phase_shifts[12] =  2.480768898038398000e-02;
                phase_shifts[13] =  3.218198903756568400e+00;
                phase_shifts[14] =  5.105848086558705900e+00;
                phase_shifts[15] =  3.848614783366913100e+00;
                phase_shifts[16] =  4.534922405866221400e+00;
                phase_shifts[17] =  1.833911423047069700e+00;
                phase_shifts[18] =  5.766544881882963700e+00;
                phase_shifts[19] =  4.489812063110711900e+00;
                phase_shifts[20] =  3.408906801581391000e+00;
            }
            else if (random_seed == 11)
            {
                phase_shifts[0] =  1.132667860480351500e+00;
                phase_shifts[1] =  1.223665511688170900e-01;
                phase_shifts[2] =  2.910487839707776500e+00;
                phase_shifts[3] =  4.554894212576069600e+00;
                phase_shifts[4] =  2.640217114369509700e+00;
                phase_shifts[5] =  3.050028410914633200e+00;
                phase_shifts[6] =  8.030422644949872200e-02;
                phase_shifts[7] =  3.062246122248716100e+00;
                phase_shifts[8] =  5.917545720207830800e+00;
                phase_shifts[9] =  5.345703204992719100e+00;
                phase_shifts[10] =  4.586502034054578100e+00;
                phase_shifts[11] =  6.832088891275668300e-01;
                phase_shifts[12] =  5.616565548761781400e+00;
                phase_shifts[13] =  5.385658971194583700e+00;
                phase_shifts[14] =  1.037269810086306000e+00;
                phase_shifts[15] =  3.973071784885684000e+00;
                phase_shifts[16] =  1.287023349278289800e-01;
                phase_shifts[17] =  7.334818919988721500e-01;
                phase_shifts[18] =  1.987794443995819500e+00;
                phase_shifts[19] =  9.921922851228006700e-01;
                phase_shifts[20] =  4.768809396779301000e+00;
            }
            else if (random_seed == 12)
            {
                phase_shifts[0] =  9.686337061529999300e-01;
                phase_shifts[1] =  4.649869379728302800e+00;
                phase_shifts[2] =  1.654457034571007900e+00;
                phase_shifts[3] =  3.353583514350031900e+00;
                phase_shifts[4] =  9.157719014108255100e-02;
                phase_shifts[5] =  5.772657702308401400e+00;
                phase_shifts[6] =  5.659358337346415800e+00;
                phase_shifts[7] =  2.099930230068142400e-01;
                phase_shifts[8] =  6.012690009399070900e+00;
                phase_shifts[9] =  8.621115919525816900e-01;
                phase_shifts[10] =  1.783346137066347000e+00;
                phase_shifts[11] =  3.808132958892010300e+00;
                phase_shifts[12] =  5.932741501518120400e+00;
                phase_shifts[13] =  5.357895422807692900e+00;
                phase_shifts[14] =  1.419518284901189600e-02;
                phase_shifts[15] =  3.274959715950131100e+00;
                phase_shifts[16] =  3.468554746338117200e+00;
                phase_shifts[17] =  3.049716233962365600e+00;
                phase_shifts[18] =  4.826329230760398700e+00;
                phase_shifts[19] =  1.009813141856190900e+00;
                phase_shifts[20] =  4.803874988019850400e+00;
            }
            else if (random_seed == 13)
            {
                phase_shifts[0] =  4.886448359475572500e+00;
                phase_shifts[1] =  1.492515503572874100e+00;
                phase_shifts[2] =  5.179094765441459600e+00;
                phase_shifts[3] =  6.067981171564244200e+00;
                phase_shifts[4] =  6.111033028633724700e+00;
                phase_shifts[5] =  2.849105648924096900e+00;
                phase_shifts[6] =  3.826726653470131600e+00;
                phase_shifts[7] =  4.872776801893367700e+00;
                phase_shifts[8] =  4.031375540680533800e+00;
                phase_shifts[9] =  4.536574331216701100e+00;
                phase_shifts[10] =  2.201409734487944100e-01;
                phase_shifts[11] =  1.875213330426413600e+00;
                phase_shifts[12] =  3.676448292799168700e-01;
                phase_shifts[13] =  5.385072721821375200e+00;
                phase_shifts[14] =  2.342710949665718400e+00;
                phase_shifts[15] =  4.271610660471637300e+00;
                phase_shifts[16] =  1.610254412133811200e+00;
                phase_shifts[17] =  2.183917184097764300e+00;
                phase_shifts[18] =  5.914217867260662000e-02;
                phase_shifts[19] =  2.251477558557948300e+00;
                phase_shifts[20] =  5.963334617450164500e+00;
            }
            else if (random_seed == 14)
            {
                phase_shifts[0] =  3.229201266315314500e+00;
                phase_shifts[1] =  4.857939295249377000e+00;
                phase_shifts[2] =  5.469058445908473200e+00;
                phase_shifts[3] =  5.056046877018226200e-02;
                phase_shifts[4] =  1.946128216239969300e+00;
                phase_shifts[5] =  6.016801745283995500e+00;
                phase_shifts[6] =  3.224007387389883600e+00;
                phase_shifts[7] =  1.999840021860475000e+00;
                phase_shifts[8] =  3.387893124464984600e+00;
                phase_shifts[9] =  1.390185803402966200e+00;
                phase_shifts[10] =  5.067271818624014700e+00;
                phase_shifts[11] =  2.150449234976434800e+00;
                phase_shifts[12] =  3.385938498661600800e+00;
                phase_shifts[13] =  3.690608307193345200e-02;
                phase_shifts[14] =  4.229541760395403700e+00;
                phase_shifts[15] =  1.319621358272680800e+00;
                phase_shifts[16] =  5.859432165674284900e+00;
                phase_shifts[17] =  2.351449111231216400e+00;
                phase_shifts[18] =  4.727587516094319900e+00;
                phase_shifts[19] =  4.794943763731478900e+00;
                phase_shifts[20] =  5.469509387271293700e+00;
            }
            else if (random_seed == 15)
            {
                phase_shifts[0] =  5.333278883951943600e+00;
                phase_shifts[1] =  1.124036246977920200e+00;
                phase_shifts[2] =  3.415741493812250500e-01;
                phase_shifts[3] =  2.271613052442061200e+00;
                phase_shifts[4] =  1.730395068203421900e+00;
                phase_shifts[5] =  3.330089625864812500e+00;
                phase_shifts[6] =  1.922145236540217400e+00;
                phase_shifts[7] =  1.913068819850232400e+00;
                phase_shifts[8] =  7.020911446719929600e-01;
                phase_shifts[9] =  1.570161812290202800e+00;
                phase_shifts[10] =  5.765638691515394300e+00;
                phase_shifts[11] =  1.659683627299319400e+00;
                phase_shifts[12] =  4.509905087013959400e+00;
                phase_shifts[13] =  5.439447983028048700e+00;
                phase_shifts[14] =  5.071029944570252500e+00;
                phase_shifts[15] =  1.322928327380049000e+00;
                phase_shifts[16] =  1.050818957975311600e+00;
                phase_shifts[17] =  2.934649139681677100e-01;
                phase_shifts[18] =  2.476976900097651600e-01;
                phase_shifts[19] =  1.258087282025922100e+00;
                phase_shifts[20] =  6.274033236596475900e+00;
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Random seed number > 16 is not supported!"
                    << std::endl);
            }
            
            // Amplitude.
            double A = A_candidates[std::stoi(A_idx_str)];
            
            // Bounds of wavenumbers for initial perturbations.
            int m_min = m_min_candidates[std::stoi(m_idx_str)];
            int m_max = m_max_candidates[std::stoi(m_idx_str)];
            
            // Characteristic length of the initial interface thickness.
            const double epsilon_i = 0.001;
            
            // species 0: SF6.
            // species 1: air.
            const double gamma_0 = 1.09312;
            const double gamma_1 = 1.39909;
            
            const double c_p_SF6 = 668.286;
            const double c_p_air = 1040.50;
            
            const double c_v_SF6 = 611.359;
            const double c_v_air = 743.697;
            
            NULL_USE(gamma_1);
            
            // Unshocked SF6.
            const double rho_unshocked = 5.97286552525647;
            const double u_unshocked   = 0.0;
            const double v_unshocked   = 0.0;
            const double p_unshocked   = 101325.0;
            
            // Shocked SF6.
            const double rho_shocked = 11.9708247309869;
            const double u_shocked   = 98.9344103891513;
            const double v_shocked   = 0.0;
            const double p_shocked   = 218005.430874;
            
            // Air.
            const double rho_air = 1.14560096494103;
            const double u_air   = 0.0;
            const double v_air   = 0.0;
            const double p_air   = 101325.0;
            
            // Shock hits the interface after 0.05 ms.
            const double L_x_shock = 0.190127254739019;
            const double L_x_interface = 0.2;
            
            for (int j = 0; j < patch_dims[1]; j++)
            {
                for (int i = 0; i < patch_dims[0]; i++)
                {
                    // Compute the linear index.
                    int idx_cell = i + j*patch_dims[0];
                    
                    // Compute the coordinates.
                    double x[2];
                    x[0] = patch_xlo[0] + (i + 0.5)*dx[0];
                    x[1] = patch_xlo[1] + (j + 0.5)*dx[1];
                    
                    double S = 0.0;
                    for (int m = m_min; m <= m_max; m++)
                    {
                        S += A*cos(2.0*M_PI*m/0.025*x[1] + phase_shifts[30 - m]);
                    }
                    
                    if (x[0] < L_x_shock)
                    {
                        rho_Y_1[idx_cell] = rho_shocked;
                        rho_Y_2[idx_cell] = 0.0;
                        rho_u[idx_cell]   = rho_shocked*u_shocked;
                        rho_v[idx_cell]   = rho_shocked*v_shocked;
                        E[idx_cell]       = p_shocked/(gamma_0 - 1.0) +
                            0.5*rho_shocked*(u_shocked*u_shocked + v_shocked*v_shocked);
                    }
                    else
                    {
                        const double f_sm = 0.5*(1.0 + erf((x[0] - (L_x_interface + S))/epsilon_i));
                        
                        // Smooth the primitive variables.
                        const double rho_Y_1_i = rho_unshocked*(1.0 - f_sm);
                        const double rho_Y_2_i = rho_air*f_sm;
                        const double u_i = u_unshocked*(1.0 - f_sm) + u_air*f_sm;
                        const double v_i = v_unshocked*(1.0 - f_sm) + v_air*f_sm;
                        const double p_i = p_unshocked*(1.0 - f_sm) + p_air*f_sm;
                        
                        const double rho_i = rho_Y_1_i + rho_Y_2_i;
                        const double Y_0_i = rho_Y_1_i/rho_i;
                        const double Y_1_i = 1.0 - Y_0_i;
                        
                        const double gamma_i = (Y_0_i*c_p_SF6 + Y_1_i*c_p_air)/
                            (Y_0_i*c_v_SF6 + Y_1_i*c_v_air);
                        
                        rho_Y_1[idx_cell] = rho_Y_1_i;
                        rho_Y_2[idx_cell] = rho_Y_2_i;
                        rho_u[idx_cell]   = rho_i*u_i;
                        rho_v[idx_cell]   = rho_i*v_i;
                        E[idx_cell]       = p_i/(gamma_i - 1.0) + 0.5*rho_i*(u_i*u_i + v_i*v_i);
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Cannot initialize data for unknown problem with name = '"
                << d_project_name
                << "'."
                << std::endl);
        }
    }
}
