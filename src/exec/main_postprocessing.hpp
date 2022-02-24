void runPostProcessing(
    HAMERS_SHARED_PTR<tbox::InputDatabase> input_db,
    const bool& is_from_restart,
    const std::string& restart_read_dirname,
    const int& restore_num)
{
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    /*
     * Retrieve "Main" section of the input database.  First, read dump
     * information, which is used for writing plot files.  Second, if
     * proper restart information was given on command line, and the
     * restart interval is non-zero, create a restart database.
     */
    
    HAMERS_SHARED_PTR<tbox::Database> main_db(input_db->getDatabase("Main"));
    
    const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));
    
    const std::string base_name = main_db->getStringWithDefault("base_name", "unnamed");
    
    const std::string log_filename = main_db->getStringWithDefault("log_filename", base_name + ".log");
    
    bool log_all_nodes = false;
    if (main_db->keyExists("log_all_nodes"))
    {
        log_all_nodes = main_db->getBool("log_all_nodes");
    }
    if (log_all_nodes)
    {
        tbox::PIO::logAllNodes(log_filename);
    }
    else
    {
        tbox::PIO::logOnlyNodeZero(log_filename);
    }
    
#ifdef _OPENMP
    tbox::plog << "Compiled with OpenMP version "
               << _OPENMP
               << ".  Running with "
               << omp_get_max_threads()
               << " threads."
               << std::endl;
#else
    tbox::plog << "Compiled without OpenMP.\n";
#endif
    
    std::string stat_dump_filename = "";
    
    if (main_db->keyExists("stat_dump_filename"))
    {
        stat_dump_filename = main_db->getString("stat_dump_filename");
    }
    else
    {
        TBOX_ERROR("Key data 'stat_dump_filename' not found in input."
            << std::endl);
    }
    
    /*
     * Get restart manager and root restart database.  If run is from
     * restart, open the restart file.
     */
    
    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
    
#ifdef HAVE_SILO
    /*
     * If SILO is present then use SILO as the file storage format
     * for this example, otherwise it will default to HDF5.
     */
    HAMERS_SHARED_PTR<tbox::SiloDatabaseFactory> silo_database_factory(
        new tbox::SiloDatabaseFactory());
    restart_manager->setDatabaseFactory(silo_database_factory);
#endif
    
    if (is_from_restart)
    {
        restart_manager->openRestartFile(
            restart_read_dirname,
            restore_num,
            mpi.getSize());
    }
    
    /*
     * Setup the timer manager to trace timing statistics during execution
     * of the code.  The list of timers is given in the tbox::TimerManager
     * section of the input file.  Timing information is stored in the
     * restart file.  Timers will automatically be initialized to their
     * previous state if the run is restarted, unless they are explicitly
     * reset using the tbox::TimerManager::resetAllTimers() routine.
     */
    
    tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
    
    /*
     * Create major algorithm and data objects which comprise application.
     * Each object is initialized either from input data or restart
     * files, or a combination of both.  Refer to each class constructor
     * for details.  For more information on the composition of objects
     * and the roles they play in this application, see comments at top of file.
     */
       
    HAMERS_SHARED_PTR<geom::CartesianGridGeometry> grid_geometry(
        new geom::CartesianGridGeometry(
            dim,
            "Cartesian grid geometry",
            input_db->getDatabase("CartesianGeometry")));
    
    HAMERS_SHARED_PTR<hier::PatchHierarchy> patch_hierarchy(
        new hier::PatchHierarchy(
            "PatchHierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy")));
    
    APPLICATION_LABEL app_label = EULER;
    
    std::string app_string = input_db->getString("Application");
    
    if (app_string == "Euler")
    {
        app_label = EULER;
    }
    else if (app_string == "Navier-Stokes")
    {
        app_label = NAVIER_STOKES;
    }
    else
    {
       TBOX_ERROR("Unkonwn application string = '"
            << app_string
            << "' found in input database."
            << std::endl); 
    }
    
    HAMERS_SHARED_PTR<Euler> Euler_app;
    HAMERS_SHARED_PTR<NavierStokes> Navier_Stokes_app;
    
    /*
     * Set bogus database for Runge Kutta level integrator if it is not in simulation mode.
     */
    
    HAMERS_SHARED_PTR<tbox::Database> RK_level_integrator_db(
        new tbox::InputDatabase("RungeKuttaLevelIntegrator"));
    
    RK_level_integrator_db->putDouble("cfl", double(0));
    RK_level_integrator_db->putDouble("cfl_init", double(0));
    RK_level_integrator_db->putBool("lag_dt_computation", false);
    RK_level_integrator_db->putBool("use_ghosts_to_compute_dt", false);
    
    HAMERS_SHARED_PTR<RungeKuttaLevelIntegrator> RK_level_integrator;
    
    const bool use_refined_timestepping = true;
    
    switch (app_label)
    {
        case EULER:
        {
            Euler_app = HAMERS_MAKE_SHARED<Euler>(
                "Euler_app",
                dim,
                input_db->getDatabase("Euler"),
                grid_geometry,
                stat_dump_filename);
            
                RK_level_integrator.reset(new RungeKuttaLevelIntegrator(
                    "Runge-Kutta level integrator",
                    RK_level_integrator_db,
                    Euler_app.get(),
                    use_refined_timestepping));
            
            break;
        }
        case NAVIER_STOKES:
        {
            Navier_Stokes_app = HAMERS_MAKE_SHARED<NavierStokes>(
                "Navier_Stokes_app",
                dim,
                input_db->getDatabase("NavierStokes"),
                grid_geometry,
                stat_dump_filename);
            
            RK_level_integrator.reset(
                new RungeKuttaLevelIntegrator(
                    "Runge-Kutta level integrator",
                    RK_level_integrator_db,
                    Navier_Stokes_app.get(),
                    use_refined_timestepping));
            
            break;
        }
    }
    
    HAMERS_SHARED_PTR<ExtendedTagAndInitialize> error_detector(
        new ExtendedTagAndInitialize(
            "ExtendedTagAndInitialize",
            RK_level_integrator.get(),
            input_db->getDatabase("ExtendedTagAndInitialize")));
    
    /*
     * Set up the clustering.
     */
    
    const std::string clustering_type =
        main_db->getStringWithDefault("clustering_type", "BergerRigoutsos");
    
    HAMERS_SHARED_PTR<mesh::BoxGeneratorStrategy> box_generator;
    
    if (clustering_type == "BergerRigoutsos")
    {
        HAMERS_SHARED_PTR<tbox::Database> abr_db(
            input_db->getDatabaseWithDefault("BergerRigoutsos", HAMERS_SHARED_PTR<tbox::Database>()));
        
        HAMERS_SHARED_PTR<mesh::BoxGeneratorStrategy> berger_rigoutsos(new mesh::BergerRigoutsos(dim, abr_db));
        box_generator = berger_rigoutsos;
    }
    else if (clustering_type == "TileClustering")
    {
        HAMERS_SHARED_PTR<tbox::Database> tc_db(
            input_db->getDatabaseWithDefault("TileClustering", HAMERS_SHARED_PTR<tbox::Database>()));
        
        HAMERS_SHARED_PTR<mesh::BoxGeneratorStrategy> tile_clustering(new mesh::TileClustering(dim, tc_db));
        box_generator = tile_clustering;
    }
    
    /*
     * Set up the load balancers.
     */
    
    HAMERS_SHARED_PTR<mesh::LoadBalanceStrategy> load_balancer;
    HAMERS_SHARED_PTR<mesh::LoadBalanceStrategy> load_balancer0;
    
    const std::string partitioner_type =
        main_db->getStringWithDefault("partitioner_type", "CascadePartitioner");
    
    if (partitioner_type == "TreeLoadBalancer")
    {
        HAMERS_SHARED_PTR<mesh::TreeLoadBalancer> tree_load_balancer(
            new mesh::TreeLoadBalancer(
                dim,
                "mesh::TreeLoadBalancer",
                input_db->getDatabaseWithDefault("TreeLoadBalancer", HAMERS_SHARED_PTR<tbox::Database>()),
                HAMERS_SHARED_PTR<tbox::RankTreeStrategy>(new tbox::BalancedDepthFirstTree)));
        
        tree_load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
        
        HAMERS_SHARED_PTR<mesh::TreeLoadBalancer> tree_load_balancer0(
            new mesh::TreeLoadBalancer(
                dim,
                "mesh::TreeLoadBalancer0",
                input_db->getDatabaseWithDefault("TreeLoadBalancer", HAMERS_SHARED_PTR<tbox::Database>()),
                HAMERS_SHARED_PTR<tbox::RankTreeStrategy>(new tbox::BalancedDepthFirstTree)));
        
        tree_load_balancer0->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
        
        load_balancer = tree_load_balancer;
        load_balancer0 = tree_load_balancer0;
    }
    else if (partitioner_type == "CascadePartitioner")
    {
        HAMERS_SHARED_PTR<mesh::CascadePartitioner> cascade_partitioner(
            new mesh::CascadePartitioner(
                dim,
                "mesh::CascadePartitioner",
                input_db->getDatabaseWithDefault("CascadePartitioner", HAMERS_SHARED_PTR<tbox::Database>())));
        
        cascade_partitioner->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
        
        HAMERS_SHARED_PTR<mesh::CascadePartitioner> cascade_partitioner0(
            new mesh::CascadePartitioner(
                dim,
                "mesh::CascadePartitioner0",
                input_db->getDatabaseWithDefault("CascadePartitioner", HAMERS_SHARED_PTR<tbox::Database>())));
        
        cascade_partitioner0->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
        
        load_balancer = cascade_partitioner;
        load_balancer0 = cascade_partitioner0;
    }
    else if (partitioner_type == "ChopAndPackLoadBalancer")
    {
        HAMERS_SHARED_PTR<mesh::ChopAndPackLoadBalancer> cap_load_balancer(
            new mesh::ChopAndPackLoadBalancer(
                dim,
                "mesh::ChopAndPackLoadBalancer",
                input_db->getDatabaseWithDefault("ChopAndPackLoadBalancer", HAMERS_SHARED_PTR<tbox::Database>())));
        
        load_balancer = cap_load_balancer;
        
        /*
         * ChopAndPackLoadBalancer has trouble on L0 for some reason.
         * Work around by using the CascadePartitioner for L0.
         */
        
        HAMERS_SHARED_PTR<mesh::CascadePartitioner> cascade_partitioner0(
            new mesh::CascadePartitioner(
                dim,
                "mesh::CascadePartitioner0",
                input_db->getDatabaseWithDefault("CascadePartitioner", HAMERS_SHARED_PTR<tbox::Database>())));
        
        cascade_partitioner0->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
        
        load_balancer0 = cascade_partitioner0;
    }
    
    HAMERS_SHARED_PTR<mesh::GriddingAlgorithm> gridding_algorithm(
        new mesh::GriddingAlgorithm(
            patch_hierarchy,
            "GriddingAlgorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            error_detector,
            box_generator,
            load_balancer,
            load_balancer0));
    
    HAMERS_SHARED_PTR<tbox::Database> time_integrator_db(
        new tbox::InputDatabase("TimeRefinementIntegrator"));
    
    time_integrator_db->putDouble("start_time", tbox::MathUtilities<double>::getMin());
    time_integrator_db->putDouble("end_time", tbox::MathUtilities<double>::getMax());
    time_integrator_db->putInteger("max_integrator_steps", tbox::MathUtilities<int>::getMax());
    
    HAMERS_SHARED_PTR<algs::TimeRefinementIntegrator> time_integrator(
        new algs::TimeRefinementIntegrator(
            "TimeRefinementIntegrator",
            time_integrator_db,
            patch_hierarchy,
            RK_level_integrator,
            gridding_algorithm));
    
    /*
     * Initialize hierarchy configuration and data on all patches.
     * Then, close restart file and write initial state for visualization.
     */
    
    double dt_now = time_integrator->initializeHierarchy();
    
    restart_manager->closeRestartFile();
    
    /*
     * After creating all objects and initializing their state, we
     * print the input database and variable database contents
     * to the log file.
     */
    
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    switch (app_label)
    {
        case EULER:
        {
            Euler_app->printClassData(tbox::pout);
            break;
        }
        case NAVIER_STOKES:
        {
            Navier_Stokes_app->printClassData(tbox::pout);
            break;
        }
    }
    
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    RK_level_integrator->printClassData(tbox::pout);
    
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    tbox::pout << std::endl;
    
    /*
     * Create timers for measuring I/O.
     */
    
    HAMERS_SHARED_PTR<tbox::Timer> t_write_stat(
        tbox::TimerManager::getManager()->getTimer("apps::main::write_stat"));
    
    tbox::pout << "Output statistics..." << std::endl;
    tbox::pout << std::endl;
    
    /*
     * Output the statistics.
     */
    
    t_write_stat->start();
    RK_level_integrator->outputDataStatistics(
        patch_hierarchy,
        time_integrator->getIntegratorTime());
    t_write_stat->stop();
    
    // /*
    //  * Output timer results.
    //  */
    
    // tbox::TimerManager::getManager()->print(tbox::plog);
    
    /*
     * At conclusion of post-processing, deallocate objects.
     */
    
    patch_hierarchy.reset();
    grid_geometry.reset();
    
    box_generator.reset();
    load_balancer.reset();
    load_balancer0.reset();
    RK_level_integrator.reset();
    error_detector.reset();
    gridding_algorithm.reset();
    time_integrator.reset();
    
    Euler_app.reset();
    Navier_Stokes_app.reset();
}