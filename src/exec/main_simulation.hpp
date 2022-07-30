void runSimulation(
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
    
    bool is_viz_dumping = false;
    std::string viz_dump_setting;
    int viz_dump_timestep_interval = 0;
    double viz_dump_time_interval = double(0);
    std::string visit_dump_dirname = "";
    int visit_dump_directory_name_zero_padding_length = 5;
    int visit_number_procs_per_file = 1;
    
    if (main_db->keyExists("viz_dump_setting"))
    {
        viz_dump_setting = main_db->getString("viz_dump_setting");
        
        if ((viz_dump_setting != "CONSTANT_TIME_INTERVAL") &&
            (viz_dump_setting != "CONSTANT_TIMESTEP_INTERVAL"))
        {
            TBOX_ERROR("Unknown viz_dump_setting string = "
                << viz_dump_setting
                << " found in input."
                << std::endl);
        }
        
        is_viz_dumping = true;
    }
    
    if (is_viz_dumping)
    {
        if (main_db->keyExists("viz_dump_interval"))
        {
            if (viz_dump_setting == "CONSTANT_TIME_INTERVAL")
            {
                viz_dump_time_interval = main_db->getDouble("viz_dump_interval");
            }
            else if (viz_dump_setting == "CONSTANT_TIMESTEP_INTERVAL")
            {
                viz_dump_timestep_interval = main_db->getInteger("viz_dump_interval");
            }
            else
            {
                TBOX_ERROR("Unknown viz_dump_setting = "
                    << viz_dump_setting
                    << "."
                    << std::endl);
            }
        }
        else
        {
            TBOX_ERROR("Key data 'viz_dump_interval' not found in input."
                << std::endl);
        }
        
        visit_dump_dirname =
            main_db->getStringWithDefault("viz_dump_dirname", base_name + ".visit");
        
        if (main_db->keyExists("visit_dump_directory_name_zero_padding_length"))
        {
            visit_dump_directory_name_zero_padding_length =
                main_db->getInteger("visit_dump_directory_name_zero_padding_length");
        }
        
        if (main_db->keyExists("visit_number_procs_per_file"))
        {
            visit_number_procs_per_file = main_db->getInteger("visit_number_procs_per_file");
        }
    }
    
    bool is_stat_dumping = false;
    std::string stat_dump_setting;
    int stat_dump_timestep_interval = 0;
    double stat_dump_time_interval = double(0);
    std::string stat_dump_filename = "";
    
    if (main_db->keyExists("stat_dump_setting"))
    {
        stat_dump_setting = main_db->getString("stat_dump_setting");
        
        if ((stat_dump_setting != "CONSTANT_TIME_INTERVAL") &&
            (stat_dump_setting != "CONSTANT_TIMESTEP_INTERVAL"))
        {
            TBOX_ERROR("Unknown stat_dump_setting string = "
                << stat_dump_setting
                << " found in input."
                << std::endl);
        }
        
        is_stat_dumping = true;
    }
    
    if (is_stat_dumping)
    {
        if (main_db->keyExists("stat_dump_interval"))
        {
            if (stat_dump_setting == "CONSTANT_TIME_INTERVAL")
            {
                stat_dump_time_interval = main_db->getDouble("stat_dump_interval");
            }
            else if (stat_dump_setting == "CONSTANT_TIMESTEP_INTERVAL")
            {
                stat_dump_timestep_interval = main_db->getInteger("stat_dump_interval");
            }
            else
            {
                TBOX_ERROR("Unknown stat_dump_setting = "
                    << stat_dump_setting
                    << "."
                    << std::endl);
            }
        }
        else
        {
            TBOX_ERROR("Key data 'stat_dump_interval' not found in input."
                << std::endl);
        }
        
        if (main_db->keyExists("stat_dump_filename"))
        {
            stat_dump_filename = main_db->getString("stat_dump_filename");
        }
        else
        {
            TBOX_ERROR("Key data 'stat_dump_filename' not found in input."
                << std::endl);
        }
    }
    
    int restart_interval = 0;
    std::string restart_write_dirname;
    bool write_restart = false;
    bool use_refined_timestepping = true;
    
    if (main_db->keyExists("restart_interval"))
    {
        restart_interval = main_db->getInteger("restart_interval");
    }
    
    restart_write_dirname =
        main_db->getStringWithDefault("restart_write_dirname",
                                      base_name + ".restart");
    
    write_restart = (restart_interval > 0) && !(restart_write_dirname.empty());
    
    if (main_db->keyExists("timestepping"))
    {
        std::string timestepping_method = main_db->getString("timestepping");
        if (timestepping_method == "SYNCHRONIZED")
        {
            use_refined_timestepping = false;
        }
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
    
    if (app_string == "EULER" || app_string == "Euler")
    {
        app_label = EULER;
    }
    else if (app_string == "NAVIER_STOKES" || app_string == "Navier-Stokes")
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
    
    HAMERS_SHARED_PTR<tbox::Database> RK_level_integrator_db = input_db->getDatabase("RungeKuttaLevelIntegrator");
    
    HAMERS_SHARED_PTR<RungeKuttaLevelIntegrator> RK_level_integrator;
    
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
    
    HAMERS_SHARED_PTR<tbox::Database> time_integrator_db;
    
    time_integrator_db = input_db->getDatabase("TimeRefinementIntegrator");
    
    HAMERS_SHARED_PTR<algs::TimeRefinementIntegrator> time_integrator(
        new algs::TimeRefinementIntegrator(
            "TimeRefinementIntegrator",
            time_integrator_db,
            patch_hierarchy,
            RK_level_integrator,
            gridding_algorithm));
    
    /*
     * Set up visualization writer(s).  Note that the Euler application
     * creates some derived data quantities so we register the Euler model
     * as a derived data writer.  If no derived data is written, this step
     * is not necessary.
     */
#ifdef HAVE_HDF5
    HAMERS_SHARED_PTR<ExtendedVisItDataWriter> visit_data_writer;
    
    switch (app_label)
    {
        case EULER:
        {
            visit_data_writer.reset(
                new ExtendedVisItDataWriter(
                    dim,
                    "Euler VisIt Writer",
                    visit_dump_dirname,
                    visit_dump_directory_name_zero_padding_length,
                    visit_number_procs_per_file));
            
            Euler_app->registerVisItDataWriter(visit_data_writer);
            
            break;
        }
        case NAVIER_STOKES:
        {
            visit_data_writer.reset(
                new ExtendedVisItDataWriter(
                    dim,
                    "Navier-Stokes VisIt Writer",
                    visit_dump_dirname,
                    visit_dump_directory_name_zero_padding_length,
                    visit_number_procs_per_file));
            
            Navier_Stokes_app->registerVisItDataWriter(visit_data_writer);
            
            break;
        }
    }
#endif
    
    /*
     * Initialize hierarchy configuration and data on all patches.
     * Then, close restart file and write initial state for visualization.
     */
    
    double dt_now = time_integrator->initializeHierarchy();
    
    double dt_const = 0.0;
    if (!(RK_level_integrator->usingRefinedTimestepping()))
    {
        dt_const = dt_now;
    }
    
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
    HAMERS_SHARED_PTR<tbox::Timer> t_write_viz(
        tbox::TimerManager::getManager()->getTimer("apps::main::write_viz"));
    
    HAMERS_SHARED_PTR<tbox::Timer> t_write_stat(
        tbox::TimerManager::getManager()->getTimer("apps::main::write_stat"));
    
    HAMERS_SHARED_PTR<tbox::Timer> t_write_restart(
        tbox::TimerManager::getManager()->getTimer(
            "apps::main::write_restart"));
    
    t_write_viz->start();
#ifdef HAVE_HDF5
    if (is_viz_dumping)
    {
        visit_data_writer->writePlotData(
            patch_hierarchy,
            time_integrator->getIntegratorStep(),
            time_integrator->getIntegratorTime());
    }
#endif
    t_write_viz->stop();
    
    t_write_stat->start();
    if (is_stat_dumping)
    {
        RK_level_integrator->outputHeaderStatistics();
        
        RK_level_integrator->computeAndOutputDataStatistics(
            patch_hierarchy,
            time_integrator->getIntegratorTime(),
            true);
    }
    t_write_stat->stop();
    
    /*
     * Time step loop. Note that the step count and integration
     * time are maintained by algs::TimeRefinementIntegrator.
     */
    
    double loop_time = time_integrator->getIntegratorTime();
    double loop_time_end = time_integrator->getEndTime();
    
    double last_viz_dump_time = floor((time_integrator->getIntegratorTime() + double(10)*HAMERS_REAL_EPSILON)/
        viz_dump_time_interval)*viz_dump_time_interval;
    bool dump_viz = true;
    
    double last_stat_dump_time = floor((time_integrator->getIntegratorTime() + double(10)*HAMERS_REAL_EPSILON)/
        stat_dump_time_interval)*stat_dump_time_interval;
    bool dump_stat = true;
    
    int iteration_num = 0;
    
    tbox::pout << "Start simulation..." << std::endl;
    tbox::pout << std::endl;
    
    while (loop_time < loop_time_end && time_integrator->stepsRemaining())
    {
        dump_viz = false;
        dump_stat = false;
        
        iteration_num = time_integrator->getIntegratorStep() + 1;
        
        // Check whether dt_now is larger than the time interval to next files dumping time.
        if ((viz_dump_setting == "CONSTANT_TIME_INTERVAL") &&
            (stat_dump_setting == "CONSTANT_TIME_INTERVAL"))
        {
            if (viz_dump_time_interval == stat_dump_time_interval)
            {
                if ((loop_time + dt_now) - (last_viz_dump_time + viz_dump_time_interval) >= -double(10)*HAMERS_REAL_EPSILON)
                {
                    dt_now = last_viz_dump_time + viz_dump_time_interval - loop_time;
                    dump_viz = true;
                    dump_stat = true;
                }
            }
            else
            {
                if ((loop_time + dt_now) - (last_viz_dump_time + viz_dump_time_interval) >= -double(10)*HAMERS_REAL_EPSILON)
                {
                    dt_now = last_viz_dump_time + viz_dump_time_interval - loop_time;
                    dump_viz = true;
                    
                    if ((loop_time + dt_now) - (last_stat_dump_time + stat_dump_time_interval) >= -double(10)*HAMERS_REAL_EPSILON)
                    {
                        dt_now = last_stat_dump_time + stat_dump_time_interval - loop_time;
                        dump_viz = false;
                        dump_stat = true;
                    }
                }
            }
        }
        else if (viz_dump_setting == "CONSTANT_TIME_INTERVAL")
        {
            if ((loop_time + dt_now) - (last_viz_dump_time + viz_dump_time_interval) >= -double(10)*HAMERS_REAL_EPSILON)
            {
                dt_now = last_viz_dump_time + viz_dump_time_interval - loop_time;
                dump_viz = true;
            }
        }
        else if (stat_dump_setting == "CONSTANT_TIME_INTERVAL")
        {
            if ((loop_time + dt_now) - (last_stat_dump_time + stat_dump_time_interval) >= -double(10)*HAMERS_REAL_EPSILON)
            {
                dt_now = last_stat_dump_time + stat_dump_time_interval - loop_time;
                dump_stat = true;
            }
        }
        
        tbox::pout << "At begining of timestep # " << iteration_num - 1 << std::endl;
        tbox::pout << "Simulation time is " << loop_time << std::endl;
        tbox::pout << "Current dt is " << dt_now << std::endl;
        
        // Advance the solution.
        double dt_new = time_integrator->advanceHierarchy(dt_now);
        
        loop_time += dt_now;
        
        if (!(RK_level_integrator->usingRefinedTimestepping()))
        {
            dt_now = dt_const;
            
            if (loop_time + dt_now > loop_time_end)
            {
                dt_now = loop_time_end - loop_time;
            }
        }
        else
        {
            dt_now = dt_new;
        }
        
        tbox::pout << "At end of timestep # " << iteration_num - 1 << std::endl;
        tbox::pout << "Simulation time is " << loop_time << std::endl;
        switch (app_label)
        {
            case EULER:
            {
                Euler_app->computeAndPrintDataStatistics(tbox::pout, patch_hierarchy, loop_time);
                break;
            }
            case NAVIER_STOKES:
            {
                Navier_Stokes_app->computeAndPrintDataStatistics(tbox::pout, patch_hierarchy, loop_time);
                break;
            }
        }
        
        /*
         * At specified intervals, write out data files for plotting.
         * If restart_interval = -1, also write restart files when writing out data
         * files for plotting
         */
#ifdef HAVE_HDF5
        if (is_viz_dumping)
        {
            if (viz_dump_setting == "CONSTANT_TIME_INTERVAL")
            {
                if (dump_viz)
                {
                    t_write_viz->start();
                    
                    visit_data_writer->writePlotData(
                        patch_hierarchy,
                        iteration_num,
                        loop_time);
                    
                    t_write_viz->stop();
                    
                    last_viz_dump_time = loop_time;
                    
                    tbox::pout << "Files for plotting are written." << std::endl;
                    
                    if ((restart_interval == -1) && !(restart_write_dirname.empty()))
                    {
                        t_write_restart->start();
                        
                        restart_manager->writeRestartFile(restart_write_dirname,
                            iteration_num);
                        
                        t_write_restart->stop();
                        
                        tbox::pout << "Files for restart are written." << std::endl;
                    }
                }
            }
            else if (viz_dump_setting == "CONSTANT_TIMESTEP_INTERVAL")
            {
                if ((iteration_num % viz_dump_timestep_interval) == 0)
                {
                    t_write_viz->start();
                    
                    visit_data_writer->writePlotData(
                        patch_hierarchy,
                        iteration_num,
                        loop_time);
                    
                    t_write_viz->stop();
                    
                    tbox::pout << "Files for plotting are written." << std::endl;
                    
                    if ((restart_interval == -1) && !(restart_write_dirname.empty()))
                    {
                        t_write_restart->start();
                        
                        restart_manager->writeRestartFile(restart_write_dirname,
                            iteration_num);
                        
                        t_write_restart->stop();
                        
                        tbox::pout << "Files for restart are written." << std::endl;
                    }
                }
            }
            else
            {
                TBOX_ERROR("Unknown viz_dump_setting = "
                    << viz_dump_setting
                    << "."
                    << std::endl); 
            }
        }
#endif
        
        if (is_stat_dumping)
        {
            if (stat_dump_setting == "CONSTANT_TIME_INTERVAL")
            {
                if (dump_stat)
                {
                    t_write_stat->start();
                    
                    RK_level_integrator->computeAndOutputDataStatistics(patch_hierarchy, loop_time, true);
                    
                    t_write_stat->stop();
                    
                    last_stat_dump_time = loop_time;
                    
                    tbox::pout << "File of statistics is updated." << std::endl;
                }
            }
            else if (stat_dump_setting == "CONSTANT_TIMESTEP_INTERVAL")
            {
                if ((iteration_num % stat_dump_timestep_interval) == 0)
                {
                    t_write_stat->start();
                    
                    RK_level_integrator->computeAndOutputDataStatistics(patch_hierarchy, loop_time, true);
                    
                    t_write_stat->stop();
                    
                    tbox::pout << "File of statistics is updated." << std::endl;
                }
            }
            else
            {
                TBOX_ERROR("Unknown stat_dump_setting = "
                    << stat_dump_setting
                    << "."
                    << std::endl); 
            }
        }
        
        /*
         * At specified intervals, write restart files.
         */
        if (write_restart)
        {
            if ((iteration_num % restart_interval) == 0)
            {
                t_write_restart->start();
                
                restart_manager->writeRestartFile(restart_write_dirname,
                    iteration_num);
                
                t_write_restart->stop();
            }
            
            tbox::pout << "Files for restart are written." << std::endl;
        }
        
        tbox::pout << "--------------------------------------------------------------------------------";
        tbox::pout << std::endl;
    }
    
#ifdef HAVE_HDF5
    if (is_viz_dumping)
    {
        if (viz_dump_setting == "CONSTANT_TIME_INTERVAL" && !dump_viz)
        {
            t_write_viz->start();
            
            visit_data_writer->writePlotData(
                patch_hierarchy,
                iteration_num,
                loop_time);
            
            t_write_viz->stop();
            
            last_viz_dump_time = loop_time;
            
            tbox::pout << "Files for plotting at last time step are written." << std::endl;
            
            if ((restart_interval == -1) && !(restart_write_dirname.empty()))
            {
                t_write_restart->start();
                
                restart_manager->writeRestartFile(restart_write_dirname,
                    iteration_num);
                
                t_write_restart->stop();
                
                tbox::pout << "Files for restart at last time step are written." << std::endl;
            }
        }
    }
#endif
    
    if (is_stat_dumping)
    {
        if (stat_dump_setting == "CONSTANT_TIME_INTERVAL" && !dump_stat)
        {
            t_write_stat->start();
            
            RK_level_integrator->computeAndOutputDataStatistics(patch_hierarchy, loop_time, true);
            
            t_write_stat->stop();
            
            last_stat_dump_time = loop_time;
            
            tbox::pout << "File of statistics at last time step is updated." << std::endl;
        }
    }
    
    tbox::plog << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    tbox::plog << std::endl;
    tbox::plog << "Error statistics:\n";
    switch (app_label)
    {
        case EULER:
        {
            Euler_app->printErrorStatistics(tbox::pout, patch_hierarchy, time_integrator->getIntegratorTime());
            break;
        }
        case NAVIER_STOKES:
        {
            Navier_Stokes_app->printErrorStatistics(tbox::pout, patch_hierarchy, time_integrator->getIntegratorTime());
            break;
        }
    }
    
    tbox::plog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    tbox::plog << std::endl;
    tbox::plog << "GriddingAlgorithm statistics:\n";
    gridding_algorithm->printStatistics();
    
    /*
     * Output timer results.
     */
    
    tbox::TimerManager::getManager()->print(tbox::plog);
    
    /*
     * At conclusion of simulation, deallocate objects.
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
#ifdef HAVE_HDF5
    visit_data_writer.reset();
#endif
    
    Euler_app.reset();
    Navier_Stokes_app.reset();
}