#/usr/bin/env python

import numpy
import subprocess
import sys
import os

convective_flux_schemes = ["WCNS5_JS_HLLC_HLL", "WCNS5_Z_HLLC_HLL", "WCNS6_LD_HLLC_HLL", "DRP4_9pt", "DRP4_11pt", "DRP4_13pt"]
L2_convergence_rates_schemes_expected = [4.8, 4.8, 5.8, 3.8, 3.8, 3.8]
num_grid_levels = 4

N_base  = 8
dx_base = 2.0/N_base
dt_base = 0.005*dx_base
num_steps_base = 16

executable_path = "../../../build_convergence_test_five_eqn_allaire/src/exec/main"

input_file_template = """
Application = "Euler"

Euler
{{
    // Name of project
    project_name = "2D convergence test five-eqn by Allaire"

    // Number of species
    num_species = 2

    // Flow model to use
    flow_model = "FIVE_EQN_ALLAIRE"

    Flow_model
    {{
        // Equation of state to use
        equation_of_state = "IDEAL_GAS"

        Equation_of_state_mixing_rules
        {{
            species_gamma = 1.6, 1.4
            species_R     = 1.0, 1.0
        }}
    }}

    // Convective flux reconstructor to use
    convective_flux_reconstructor = "{:s}"

    Convective_flux_reconstructor
    {{
        stencil_width = {:d}
    }}

    Boundary_data
    {{
    }}
}}

Main
{{
    // Dimension of problem
    dim = 2

    // Base name of log file
    base_name = "2D_convergence_test_five_eqn_allaire"

    // Whether all nodes log to individual files,
    // if false only node 0 will log
    log_all_nodes = TRUE

    // Visualization dump parameters
    // Setting to dump viz output: CONSTANT_TIME_INTERVAL, CONSTANT_TIMESTEP_INTERVAL
    viz_dump_setting = "CONSTANT_TIME_INTERVAL"
    // Frequency at which to dump viz output (0 to turn off)
    viz_dump_interval = 0.02
    // Name of directory in which to place viz output
    viz_dump_dirname = "viz_2D_convergence_test_five_eqn_allaire"
    // Number of processors which write to each viz file
    visit_number_procs_per_file = 1

    // Restart dump parameters
    // Frequency at which to dump restart output
    // (-1 to be as same as viz_dump_interval, 0 to turn off)
    restart_interval = -1
}}

CartesianGeometry
{{
    domain_boxes = [(0, 0), ({:d}, {:d})] // Lower and upper indices of compuational domain
    x_lo         = -1.0, -1.0 // Lower end of computational domain
    x_up         =  1.0,  1.0 // Upper end of computational domain

    // Periodic_dimension. A non-zero value indicates that the direction is periodic
    periodic_dimension = 1, 1
}}

ExtendedTagAndInitialize
{{
}}

PatchHierarchy
{{
    // Maximum number of levels in hierarchy
    max_levels = 1

    ratio_to_coarser
    {{
    }}

    largest_patch_size
    {{
        level_0 = 8, 8
    }}

    smallest_patch_size
    {{
        level_0 = 8, 8
    }}
}}

GriddingAlgorithm
{{
}}

RungeKuttaLevelIntegrator
{{
    use_cfl = FALSE
    dt      = {:24.16e}
}}

TimeRefinementIntegrator
{{
    start_time           = 0.0e0   // Initial simulation time
    end_time             = 2.0e-2  // Final simulation time
    grow_dt              = 1.0e0   // Growth factor for timesteps
    max_integrator_steps = {:d}
}}

TimerManager
{{
    print_threshold = 0.01
    timer_list      = "apps::main::*",
                      "apps::Euler::*",
                      "algs::GriddingAlgorithm::*"
}}
"""
numpy.set_printoptions(formatter={'float': '{:24.16e}'.format})

L1_errors_schemes   = []
L2_errors_schemes   = []
Linf_errors_schemes = []

L1_convergence_rates_schemes   = []
L2_convergence_rates_schemes   = []
Linf_convergence_rates_schemes = []

for scheme in convective_flux_schemes:
    stencil_width = 0
    if scheme.find("DRP4_") != -1:
        stencil_width_str = scheme.split("_")[1]
        stencil_width = int(stencil_width_str.split("pt")[0])
        
        scheme = scheme.split("_")[0]

    L1_errors   = []
    L2_errors   = []
    Linf_errors = []

    L1_convergence_rates   = []
    L2_convergence_rates   = []
    Linf_convergence_rates = []

    for level in range(num_grid_levels):
        factor = 2**level

        N_x_level       = N_base*factor
        N_y_level       = N_base*factor
        dx_level        = dx_base/factor
        dt_level        = dt_base/factor
        num_steps_level = num_steps_base*factor
        
        level_dir = "./" + scheme + "_N_" + str(N_x_level)
        if not os.path.isdir(level_dir):
            os.mkdir(level_dir)

        input_file = open(level_dir + "/input_2D_convergence_single_species.txt", "w")

        input_file.write(input_file_template.format(scheme, stencil_width, (N_x_level - 1), (N_y_level - 1), dt_level, num_steps_level))

        input_file.close()

        os.chdir(level_dir)
        os.system(executable_path + " input_2D_convergence_single_species.txt")
        output_L1 = subprocess.check_output('grep "L1_error" 2D_convergence_test_five_eqn_allaire.log.0000000', shell=True)
        output_L2 = subprocess.check_output('grep "L2_error" 2D_convergence_test_five_eqn_allaire.log.0000000', shell=True)
        output_Linf = subprocess.check_output('grep "Linf_error" 2D_convergence_test_five_eqn_allaire.log.0000000', shell=True)

        L1_error = float(output_L1.split()[1])
        L1_errors.append(L1_error)

        L2_error = float(output_L2.split()[1])
        L2_errors.append(L2_error)
        
        Linf_error = float(output_Linf.split()[1])
        Linf_errors.append(Linf_error)

        if level == 0:
            L1_convergence_rates.append(0.0)
            L2_convergence_rates.append(0.0)
            Linf_convergence_rates.append(0.0)
        else:
            L1_convergence_rate = numpy.log2(L1_errors[level - 1]/L1_errors[level])
            L1_convergence_rates.append(L1_convergence_rate)
            L2_convergence_rate = numpy.log2(L2_errors[level - 1]/L2_errors[level])
            L2_convergence_rates.append(L2_convergence_rate)
            Linf_convergence_rate = numpy.log2(Linf_errors[level - 1]/Linf_errors[level])
            Linf_convergence_rates.append(Linf_convergence_rate)
        
        os.chdir("..")

    L1_errors_schemes.append(L1_errors)
    L2_errors_schemes.append(L2_errors)
    Linf_errors_schemes.append(Linf_errors)

    L1_convergence_rates_schemes.append(L1_convergence_rates)
    L2_convergence_rates_schemes.append(L2_convergence_rates)
    Linf_convergence_rates_schemes.append(Linf_convergence_rates)

print("")
print("--------------------------------------------------------------------------------")
print("                 Summay of convergence tests of different schemes               ")
print("--------------------------------------------------------------------------------")

count = 0
for scheme in convective_flux_schemes:
    print("Convergence test of scheme: " + scheme)
    
    L1_errors   = L1_errors_schemes[count]
    L2_errors   = L2_errors_schemes[count]
    Linf_errors = Linf_errors_schemes[count]
    
    L1_convergence_rates   = L1_convergence_rates_schemes[count]
    L2_convergence_rates   = L2_convergence_rates_schemes[count]
    Linf_convergence_rates = Linf_convergence_rates_schemes[count]
    
    print("  Number of grid points |        L1 error        |        L2 error        |       Linf error       |   L1 convergence rate  |   L2 convergence rate  |  Linf convergence rate ")
    for level in range(num_grid_levels):
        factor = 2**level
        N_x_level = N_base*factor

        if level == 0:
            print("  %21s %24.16e %24.16e %24.16e" % ( str(N_x_level) + "^2        ", L1_errors[level], L2_errors[level], Linf_errors[level] ) )
        else:
            print("  %21s %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e" % ( str(N_x_level) + "^2        ", L1_errors[level], L2_errors[level], Linf_errors[level], L1_convergence_rates[level], L2_convergence_rates[level], Linf_convergence_rates[level] ) )

    assert L2_convergence_rates[num_grid_levels - 1] > L2_convergence_rates_schemes_expected[count]
    count += 1
