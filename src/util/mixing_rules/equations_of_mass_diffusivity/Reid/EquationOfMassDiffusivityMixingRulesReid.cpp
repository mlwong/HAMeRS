#include "util/mixing_rules/equations_of_mass_diffusivity/Reid/EquationOfMassDiffusivityMixingRulesReid.hpp"

#include <cmath>

#define EPSILON HAMERS_EPSILON

EquationOfMassDiffusivityMixingRulesReid::EquationOfMassDiffusivityMixingRulesReid(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const boost::shared_ptr<tbox::Database>& equation_of_mass_diffusivity_mixing_rules_db):
        EquationOfMassDiffusivityMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_mass_diffusivity_mixing_rules_db)
{
    /*
     * Get the Lennard-Jones energy parameter of each species from the database.
     */
    
    if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("species_epsilon_by_k"))
    {
        size_t species_epsilon_by_k_array_size =
            equation_of_mass_diffusivity_mixing_rules_db->getArraySize("species_epsilon_by_k");
        if (static_cast<int>(species_epsilon_by_k_array_size) == d_num_species)
        {
            d_species_epsilon_by_k =
                equation_of_mass_diffusivity_mixing_rules_db->getDoubleVector("species_epsilon_by_k");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_epsilon_by_k' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("d_species_epsilon_by_k"))
    {
        size_t species_epsilon_by_k_array_size =
            equation_of_mass_diffusivity_mixing_rules_db->getArraySize("d_species_epsilon_by_k");
        if (static_cast<int>(species_epsilon_by_k_array_size) == d_num_species)
        {
            d_species_epsilon_by_k =
                equation_of_mass_diffusivity_mixing_rules_db->getDoubleVector("d_species_epsilon_by_k");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_epsilon_by_k' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_epsilon_by_k'/'d_species_epsilon_by_k'"
            << "not found in data for equation of mass diffusivity mixing rules."
            << std::endl);
    }
    
    /*
     * Get the collision diameter of each species from the database.
     */
    
    if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("species_sigma"))
    {
        size_t species_sigma_array_size =
            equation_of_mass_diffusivity_mixing_rules_db->getArraySize("species_sigma");
        if (static_cast<int>(species_sigma_array_size) == d_num_species)
        {
            d_species_sigma =
                equation_of_mass_diffusivity_mixing_rules_db->getDoubleVector("species_sigma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_sigma' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("d_species_sigma"))
    {
        size_t species_sigma_array_size =
            equation_of_mass_diffusivity_mixing_rules_db->getArraySize("d_species_sigma");
        if (static_cast<int>(species_sigma_array_size) == d_num_species)
        {
            d_species_sigma =
                equation_of_mass_diffusivity_mixing_rules_db->getDoubleVector("d_species_sigma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_sigma' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_sigma'/'d_species_sigma'"
            << "not found in data for equation of mass diffusivity mixing rules."
            << std::endl);
    }
    
    /*
     * Get the molecular weight of each species from the database.
     */
    
    if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("species_M"))
    {
        size_t species_M_array_size =
            equation_of_mass_diffusivity_mixing_rules_db->getArraySize("species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M =
                equation_of_mass_diffusivity_mixing_rules_db->getDoubleVector("species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_M' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("d_species_M"))
    {
        size_t species_M_array_size =
            equation_of_mass_diffusivity_mixing_rules_db->getArraySize("d_species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M =
                equation_of_mass_diffusivity_mixing_rules_db->getDoubleVector("d_species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_M' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_M'/'d_species_M'"
            << "not found in data for equation of mass diffusivity mixing rules."
            << std::endl);
    }
}


/*
 * Print all characteristics of the equation of mass diffusivity class.
 */
void
EquationOfMassDiffusivityMixingRulesReid::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfMassDiffusivityMixingRulesReid object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfMassDiffusivityMixingRulesReid: this = "
       << (EquationOfMassDiffusivityMixingRulesReid *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    /*
     * Print the Lennard-Jones energy parameter of each species.
     */
    
    os << "d_species_epsilon_by_k = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_epsilon_by_k[si] << ", ";
    }
    os << d_species_epsilon_by_k[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the collision diameter of each species.
     */
    
    os << "d_species_sigma = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_sigma[si] << ", ";
    }
    os << d_species_sigma[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the molecular weight of each species.
     */
    
    os << "d_species_M = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_M[si] << ", ";
    }
    os << d_species_M[d_num_species - 1];
    os << std::endl;
}


/*
 * Put the characteristics of the equation of mass diffusivity mixing rules class into the restart
 * database.
 */
void
EquationOfMassDiffusivityMixingRulesReid::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_epsilon_by_k", d_species_epsilon_by_k);
    restart_db->putDoubleVector("d_species_sigma", d_species_sigma);
    restart_db->putDoubleVector("d_species_M", d_species_M);
}


/*
 * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfMassDiffusivityMixingRulesReid::getMassDiffusivities(
    std::vector<double*>& mass_diffusivities,
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(static_cast<int>(mass_diffusivities.size()) == d_num_species);
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    if (d_num_species > 1)
    {
        /*
         * Initialize the containers and pointers to the containers for the molecular properties
         * of species i and j.
         */
        
        std::vector<double> species_molecular_properties_i;
        std::vector<double> species_molecular_properties_j;
        std::vector<double*> species_molecular_properties_ptr_i;
        std::vector<double*> species_molecular_properties_ptr_j;
        std::vector<const double*> species_molecular_properties_const_ptr_i;
        std::vector<const double*> species_molecular_properties_const_ptr_j;
        
        const int num_molecular_properties = getNumberOfSpeciesMolecularProperties();
        
        species_molecular_properties_i.resize(num_molecular_properties);
        species_molecular_properties_j.resize(num_molecular_properties);
        species_molecular_properties_ptr_i.reserve(num_molecular_properties);
        species_molecular_properties_ptr_j.reserve(num_molecular_properties);
        species_molecular_properties_const_ptr_i.reserve(num_molecular_properties);
        species_molecular_properties_const_ptr_j.reserve(num_molecular_properties);
        
        for (int mi = 0; mi < num_molecular_properties; mi++)
        {
            species_molecular_properties_ptr_i.push_back(&species_molecular_properties_i[mi]);
            species_molecular_properties_ptr_j.push_back(&species_molecular_properties_j[mi]);
            species_molecular_properties_const_ptr_i.push_back(&species_molecular_properties_i[mi]);
            species_molecular_properties_const_ptr_j.push_back(&species_molecular_properties_j[mi]);
        }
        
        /*
         * Get the binary mass diffusivity for each pair of species.
         */
        
        std::vector<double> D_ij((d_num_species - 1)*d_num_species/2);
        
        for (int i = 0; i < d_num_species; i++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr_i, i);
            
            for (int j = i + 1; j < d_num_species; j++)
            {
                const int idx_ij = (d_num_species - 1)*d_num_species/2 -
                    (d_num_species - 1 - i)*(d_num_species - i)/2 +
                    (j - (i + 1));
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr_j, j);
                
                D_ij[idx_ij] = getMassDiffusivity(
                    pressure,
                    temperature,
                    species_molecular_properties_const_ptr_i,
                    species_molecular_properties_const_ptr_j);
            }
        }
        
        /*
         * Compute the mole fractions.
         */
        
        std::vector<double> X;
        X.reserve(d_num_species);
        
        double sum = double(0);
        
        if (static_cast<int>(mass_fractions.size()) == d_num_species)
        {
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr_i, si);
                X.push_back((*(mass_fractions[si]))/(species_molecular_properties_i[2]));
                X[si] = std::max(double(0), X[si]);
                sum += X[si];
            }
        }
        else if (static_cast<int>(mass_fractions.size()) == d_num_species - 1)
        {
            double Y_last = double(1);
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr_i, si);
                X.push_back((*(mass_fractions[si]))/(species_molecular_properties_i[2]));
                X[si] = std::max(double(0), X[si]);
                sum += X[si];
                
                // Compute the mass fraction of the last species.
                Y_last -= *(mass_fractions[si]);
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr_i, d_num_species - 1);
            X.push_back(Y_last/(species_molecular_properties_i[2]));
            X[d_num_species - 1] = std::max(double(0), X[d_num_species - 1]);
            sum += X[d_num_species - 1];
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of mass fractions provided is not"
                << " equal to the total number of species or (total number of species - 1)."
                << std::endl);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            X[si] = X[si]/sum;
        }
        
        /*
         * Compute the effective binary diffusivity for each species.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            double& D = *(mass_diffusivities[si]);
            D = double(0);
            
            for (int sj = 0; sj < d_num_species; sj++)
            {
                if (si != sj)
                {
                    int idx = 0;
                    if (sj > si)
                    {
                        idx = (d_num_species - 1)*d_num_species/2 -
                            (d_num_species - 1 - si)*(d_num_species - si)/2 +
                            (sj - (si + 1));
                    }
                    else
                    {
                        idx = (d_num_species - 1)*d_num_species/2 -
                            (d_num_species - 1 - sj)*(d_num_species - sj)/2 +
                            (si - (sj + 1));
                    }
                    
                    D += (X[sj] + double(EPSILON))/(D_ij[idx] + double(EPSILON));
                }
            }
            
            D = (double(1) - X[si] + double(EPSILON))/D;
        }
    }
    else
    {
        *(mass_diffusivities[0]) = double(0);
    }
}


/*
 * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfMassDiffusivityMixingRulesReid::computeMassDiffusivities(
    boost::shared_ptr<pdat::CellData<double> >& data_mass_diffusivities,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_mass_diffusivities->getDepth() == d_num_species);
    TBOX_ASSERT(data_mass_fractions->getDepth() == d_num_species ||
                data_mass_fractions->getDepth() == d_num_species - 1);
#endif
    
    if (d_num_species == 1)
    {
        if (domain.empty())
        {
            data_mass_diffusivities->fillAll(double(0));
        }
        else
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(data_mass_diffusivities->getGhostBox().contains(domain));
#endif
            
            data_mass_diffusivities->fillAll(double(0), domain);
        }
        
        return;
    }
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_mass_diffusivities = data_mass_diffusivities->getGhostBox();
    const hier::IntVector ghostcell_dims_mass_diffusivities = ghost_box_mass_diffusivities.numberCells();
    
    const hier::Box ghost_box_pressure = data_pressure->getGhostBox();
    const hier::IntVector ghostcell_dims_pressure = ghost_box_pressure.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
    
    // Delcare data containers for binary mass diffusivities and mole fractions.
    boost::shared_ptr<pdat::CellData<double> > data_binary_mass_diffusivities;
    boost::shared_ptr<pdat::CellData<double> > data_mole_fractions;
    boost::shared_ptr<pdat::CellData<double> > data_sum;
    
    // Declare data container for last volume fraction.
    boost::shared_ptr<pdat::CellData<double> > data_mass_fractions_last;
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     * Also, get the offsets of all data and dimensions of the ghost cell box for binary
     * mass diffusivities and mole fractions and allocate memory.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_mass_diffusivities(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_mass_fractions(d_dim);
    hier::IntVector offset_min(d_dim);
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells.
        const hier::IntVector num_ghosts_mass_diffusivities = data_mass_diffusivities->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
        
        // Get the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = data_mass_diffusivities->getBox();
        const hier::IntVector interior_dims = interior_box.numberCells();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimensions of the ghost cell box for binary
         * mass diffusivities and mole fractions.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mass_diffusivities;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_min = num_ghosts_min;
        offset_mass_diffusivities = num_ghosts_mass_diffusivities;
        offset_pressure = num_ghosts_pressure;
        offset_temperature = num_ghosts_temperature;
        offset_mass_fractions = num_ghosts_mass_fractions;
        
        ghostcell_dims_min = interior_dims + num_ghosts_min*2;
        
        data_binary_mass_diffusivities = boost::make_shared<pdat::CellData<double> >(
            interior_box, (d_num_species - 1)*d_num_species/2, num_ghosts_min);
        
        data_mole_fractions = boost::make_shared<pdat::CellData<double> >(
            interior_box, d_num_species, num_ghosts_min);
        
        data_sum = boost::make_shared<pdat::CellData<double> >(
            interior_box, 1, num_ghosts_min);
        
        if (data_mass_fractions->getDepth() == d_num_species - 1)
        {
            data_mass_fractions_last = boost::make_shared<pdat::CellData<double> >(
                interior_box, 1, num_ghosts_mass_fractions);
        }
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(data_mass_diffusivities->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_min = hier::IntVector::getZero(d_dim);
        offset_mass_diffusivities = domain.lower() - ghost_box_mass_diffusivities.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
        
        ghostcell_dims_min = domain_dims;
        
        data_binary_mass_diffusivities = boost::make_shared<pdat::CellData<double> >(
            domain, (d_num_species - 1)*d_num_species/2, hier::IntVector::getZero(d_dim));
        
        data_mole_fractions = boost::make_shared<pdat::CellData<double> >(
            domain, d_num_species, hier::IntVector::getZero(d_dim));
        
        data_sum = boost::make_shared<pdat::CellData<double> >(
            domain, 1, hier::IntVector::getZero(d_dim));
        
        if (data_mass_fractions->getDepth() == d_num_species - 1)
        {
            data_mass_fractions_last = boost::make_shared<pdat::CellData<double> >(
                domain, 1, hier::IntVector::getZero(d_dim));
        }
    }
    
    if (domain.empty())
    {
        data_sum->fillAll(double(0));
    }
    else
    {
        data_sum->fillAll(double(0), domain);
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    std::vector<double*> D;
    D.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        D.push_back(data_mass_diffusivities->getPointer(si));
    }
    double* p = data_pressure->getPointer(0);
    double* T = data_temperature->getPointer(0);
    std::vector<double*> D_ij;
    D_ij.reserve((d_num_species - 1)*d_num_species/2);
    for (int si = 0; si < (d_num_species - 1)*d_num_species/2; si++)
    {
        D_ij.push_back(data_binary_mass_diffusivities->getPointer(si));
    }
    std::vector<double*> X;
    X.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        X.push_back(data_mole_fractions->getPointer(si));
    }
    double* sum = data_sum->getPointer(0);
    
    /*
     * Compute the binary mass diffusivities.
     */
    
    {
        const double A = double(1.06036);
        const double B = double(-0.1561);
        const double C = double(0.19300);
        const double D = double(-0.47635);
        const double E = double(1.03587);
        const double F = double(-1.52996);
        const double G = double(1.76474);
        const double H = double(-3.89411);
        
        /*
         * Initialize the containers and pointers to the containers for the molecular properties of
         * species i and j.
         */
        
        std::vector<double> species_molecular_properties_i;
        std::vector<double> species_molecular_properties_j;
        std::vector<double*> species_molecular_properties_ptr_i;
        std::vector<double*> species_molecular_properties_ptr_j;
        
        const int num_molecular_properties = getNumberOfSpeciesMolecularProperties();
        
        species_molecular_properties_i.resize(num_molecular_properties);
        species_molecular_properties_j.resize(num_molecular_properties);
        species_molecular_properties_ptr_i.reserve(num_molecular_properties);
        species_molecular_properties_ptr_j.reserve(num_molecular_properties);
        
        for (int mi = 0; mi < num_molecular_properties; mi++)
        {
            species_molecular_properties_ptr_i.push_back(&species_molecular_properties_i[mi]);
            species_molecular_properties_ptr_j.push_back(&species_molecular_properties_j[mi]);
        }
                
        for (int si = 0; si < d_num_species; si++)
        {
            getSpeciesMolecularProperties(species_molecular_properties_ptr_i, si);
            
            const double& epsilon_by_k_i = species_molecular_properties_i[0];
            const double& sigma_i = species_molecular_properties_i[1];
            const double& M_i = species_molecular_properties_i[2];
            
            for (int sj = si + 1; sj < d_num_species; sj++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr_j, sj);
                
                const double& epsilon_by_k_j = species_molecular_properties_j[0];
                const double& sigma_j = species_molecular_properties_j[1];
                const double& M_j = species_molecular_properties_j[2];
                
                const double T_epsilon_ij = sqrt(epsilon_by_k_i*epsilon_by_k_j);
                const double sigma_ij = double(1)/double(2)*(sigma_i + sigma_j);
                const double sigma_ij_sq = sigma_ij*sigma_ij;
                const double M_ij = double(2)/(double(1)/M_i + double(1)/M_j);
                const double M_ij_sqrt = sqrt(M_ij);
                
                const int idx_ij = (d_num_species - 1)*d_num_species/2 -
                    (d_num_species - 1 - si)*(d_num_species - si)/2 +
                    (sj - (si + 1));
                
                if (d_dim == tbox::Dimension(1))
                {
                    /*
                     * Get the local lower index, numbers of cells in each dimension and offsets.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_dim_0 = domain_dims[0];
                    
                    const int offset_0_min = offset_min[0];
                    const int offset_0_pressure = offset_pressure[0];
                    const int offset_0_temperature = offset_temperature[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = i + offset_0_min;
                        const int idx_pressure = i + offset_0_pressure;
                        const int idx_temperature = i + offset_0_temperature;
                        
                        const double T_star_ij = T[idx_temperature]/T_epsilon_ij;
                        const double Omega_D_ij = A*pow(T_star_ij, B) + C*exp(D*T_star_ij) + E*exp(F*T_star_ij) +
                            G*exp(H*T_star_ij);
                        
                        D_ij[idx_ij][idx_min] = double(0.0266)*pow(T[idx_temperature], double(3)/double(2))/
                            (Omega_D_ij*p[idx_pressure]*M_ij_sqrt*sigma_ij_sq);
                    }
                }
                else if (d_dim == tbox::Dimension(2))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and offsets.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    
                    const int offset_0_min = offset_min[0];
                    const int offset_1_min = offset_min[1];
                    const int ghostcell_dim_0_min = ghostcell_dims_min[0];
                    
                    const int offset_0_pressure = offset_pressure[0];
                    const int offset_1_pressure = offset_pressure[1];
                    const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                    
                    const int offset_0_temperature = offset_temperature[0];
                    const int offset_1_temperature = offset_temperature[1];
                    const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                    
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min;
                            
                            const int idx_pressure = (i + offset_0_pressure) +
                                (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                            
                            const int idx_temperature = (i + offset_0_temperature) +
                                (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                            
                            const double T_star_ij = T[idx_temperature]/T_epsilon_ij;
                            const double Omega_D_ij = A*pow(T_star_ij, B) + C*exp(D*T_star_ij) + E*exp(F*T_star_ij) +
                                G*exp(H*T_star_ij);
                            
                            D_ij[idx_ij][idx_min] = double(0.0266)*pow(T[idx_temperature], double(3)/double(2))/
                                (Omega_D_ij*p[idx_pressure]*M_ij_sqrt*sigma_ij_sq);
                        }
                    }
                }
                else if (d_dim == tbox::Dimension(3))
                {
                    /*
                     * Get the local lower indices, numbers of cells in each dimension and offsets.
                     */
                    
                    const int domain_lo_0 = domain_lo[0];
                    const int domain_lo_1 = domain_lo[1];
                    const int domain_lo_2 = domain_lo[2];
                    const int domain_dim_0 = domain_dims[0];
                    const int domain_dim_1 = domain_dims[1];
                    const int domain_dim_2 = domain_dims[2];
                    
                    const int offset_0_min = offset_min[0];
                    const int offset_1_min = offset_min[1];
                    const int offset_2_min = offset_min[2];
                    const int ghostcell_dim_0_min = ghostcell_dims_min[0];
                    const int ghostcell_dim_1_min = ghostcell_dims_min[1];
                    
                    const int offset_0_pressure = offset_pressure[0];
                    const int offset_1_pressure = offset_pressure[1];
                    const int offset_2_pressure = offset_pressure[2];
                    const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                    const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
                    
                    const int offset_0_temperature = offset_temperature[0];
                    const int offset_1_temperature = offset_temperature[1];
                    const int offset_2_temperature = offset_temperature[2];
                    const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                    const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
                    
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_min = (i + offset_0_min) +
                                    (j + offset_1_min)*ghostcell_dim_0_min +
                                    (k + offset_2_min)*ghostcell_dim_0_min*
                                        ghostcell_dim_1_min;
                                
                                const int idx_pressure = (i + offset_0_pressure) +
                                    (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                                    (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                                        ghostcell_dim_1_pressure;
                                
                                const int idx_temperature = (i + offset_0_temperature) +
                                    (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                                    (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                                        ghostcell_dim_1_temperature;
                                
                                const double T_star_ij = T[idx_temperature]/T_epsilon_ij;
                                const double Omega_D_ij = A*pow(T_star_ij, B) + C*exp(D*T_star_ij) +
                                    E*exp(F*T_star_ij) + G*exp(H*T_star_ij);
                                
                                D_ij[idx_ij][idx_min] = double(0.0266)*pow(T[idx_temperature], double(3)/double(2))/
                                    (Omega_D_ij*p[idx_pressure]*M_ij_sqrt*sigma_ij_sq);
                            }
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Compute the mole fractions.
     */
    
    if (data_mass_fractions->getDepth() == d_num_species)
    {
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<double*> Y;
        Y.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        /*
         * Initialize the container and pointer to the container for the molecular properties.
         */
        
        std::vector<double> species_molecular_properties;
        std::vector<double*> species_molecular_properties_ptr;
        
        const int num_molecular_properties = getNumberOfSpeciesMolecularProperties();
        
        species_molecular_properties.resize(num_molecular_properties);
        species_molecular_properties_ptr.reserve(num_molecular_properties);
        
        for (int mi = 0; mi < num_molecular_properties; mi++)
        {
            species_molecular_properties_ptr.push_back(&species_molecular_properties[mi]);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_min = offset_min[0];
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_min = i + offset_0_min;
                    const int idx_mass_fractions = i + offset_0_mass_fractions;
                    
                    X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[2];
                    X[si][idx_min] = std::max(double(0), X[si][idx_min]);
                    sum[idx_min] += X[si][idx_min];
                }
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            const int offset_1_mass_fractions = offset_mass_fractions[1];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                        
                        X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[2];
                        X[si][idx_min] = std::max(double(0), X[si][idx_min]);
                        sum[idx_min] += X[si][idx_min];
                    }
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int offset_2_min = offset_min[2];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            const int ghostcell_dim_1_min = ghostcell_dims_min[1];
            
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            const int offset_1_mass_fractions = offset_mass_fractions[1];
            const int offset_2_mass_fractions = offset_mass_fractions[2];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
            
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min +
                                (k + offset_2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                                (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                    ghostcell_dim_1_mass_fractions;
                            
                            X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[2];
                            X[si][idx_min] = std::max(double(0), X[si][idx_min]);
                            sum[idx_min] += X[si][idx_min];
                        }
                    }
                }
            }
        }
    }
    else if (data_mass_fractions->getDepth() == d_num_species - 1)
    {
        data_mass_fractions_last->fillAll(double(1));
        
        /*
         * Get the pointers to the cell data of mass fractions.
         */
        
        std::vector<double*> Y;
        Y.reserve(d_num_species - 1);
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Y.push_back(data_mass_fractions->getPointer(si));
        }
        
        double* Y_last = data_mass_fractions_last->getPointer(0);
        
        /*
         * Initialize the container and pointer to the container for the molecular properties.
         */
        
        std::vector<double> species_molecular_properties;
        std::vector<double*> species_molecular_properties_ptr;
        
        const int num_molecular_properties = getNumberOfSpeciesMolecularProperties();
        
        species_molecular_properties.resize(num_molecular_properties);
        species_molecular_properties_ptr.reserve(num_molecular_properties);
        
        for (int mi = 0; mi < num_molecular_properties; mi++)
        {
            species_molecular_properties_ptr.push_back(&species_molecular_properties[mi]);
        }
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_min = offset_min[0];
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_min = i + offset_0_min;
                    const int idx_mass_fractions = i + offset_0_mass_fractions;
                    
                    X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[2];
                    X[si][idx_min] = std::max(double(0), X[si][idx_min]);
                    sum[idx_min] += X[si][idx_min];
                    
                    // Compute the mass fraction of the last species.
                    Y_last[idx_min] -= Y[si][idx_mass_fractions];
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_min = i + offset_0_min;
                
                X[d_num_species - 1][idx_min] = Y_last[idx_min]/species_molecular_properties[2];
                X[d_num_species - 1][idx_min] = std::max(double(0), X[d_num_species - 1][idx_min]);
                sum[idx_min] += X[d_num_species - 1][idx_min];
            }
        }
        else if (d_dim == tbox::Dimension(2))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            const int offset_1_mass_fractions = offset_mass_fractions[1];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min;
                        
                        const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                            (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                        
                        X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[2];
                        X[si][idx_min] = std::max(double(0), X[si][idx_min]);
                        sum[idx_min] += X[si][idx_min];
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_min] -= Y[si][idx_mass_fractions];
                    }
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min;
                    
                    X[d_num_species - 1][idx_min] = Y_last[idx_min]/species_molecular_properties[2];
                    X[d_num_species - 1][idx_min] = std::max(double(0), X[d_num_species - 1][idx_min]);
                    sum[idx_min] += X[d_num_species - 1][idx_min];
                }
            }
        }
        else if (d_dim == tbox::Dimension(3))
        {
            /*
             * Get the local lower indices, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_lo_1 = domain_lo[1];
            const int domain_lo_2 = domain_lo[2];
            const int domain_dim_0 = domain_dims[0];
            const int domain_dim_1 = domain_dims[1];
            const int domain_dim_2 = domain_dims[2];
            
            const int offset_0_min = offset_min[0];
            const int offset_1_min = offset_min[1];
            const int offset_2_min = offset_min[2];
            const int ghostcell_dim_0_min = ghostcell_dims_min[0];
            const int ghostcell_dim_1_min = ghostcell_dims_min[1];
            
            const int offset_0_mass_fractions = offset_mass_fractions[0];
            const int offset_1_mass_fractions = offset_mass_fractions[1];
            const int offset_2_mass_fractions = offset_mass_fractions[2];
            const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
            const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr, si);
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min +
                                (k + offset_2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                                (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                (k + offset_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                    ghostcell_dim_1_mass_fractions;
                            
                            X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[2];
                            X[si][idx_min] = std::max(double(0), X[si][idx_min]);
                            sum[idx_min] += X[si][idx_min];
                            
                            // Compute the mass fraction of the last species.
                            Y_last[idx_min] -= Y[si][idx_mass_fractions];
                        }
                    }
                }
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min +
                            (k + offset_2_min)*ghostcell_dim_0_min*
                                ghostcell_dim_1_min;
                        
                        X[d_num_species - 1][idx_min] = Y_last[idx_min]/species_molecular_properties[2];
                        X[d_num_species - 1][idx_min] = std::max(double(0), X[d_num_species - 1][idx_min]);
                        sum[idx_min] += X[d_num_species - 1][idx_min];
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of components in the data of mass fractions provided is not"
            << " equal to the total number of species or (total number of species - 1)."
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_min = offset_min[0];
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_min = i + offset_0_min;
                
                X[si][idx_min] = X[si][idx_min]/sum[idx_min];
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int offset_0_min = offset_min[0];
        const int offset_1_min = offset_min[1];
        const int ghostcell_dim_0_min = ghostcell_dims_min[0];
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min;
                    
                    X[si][idx_min] = X[si][idx_min]/sum[idx_min];
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int offset_0_min = offset_min[0];
        const int offset_1_min = offset_min[1];
        const int offset_2_min = offset_min[2];
        const int ghostcell_dim_0_min = ghostcell_dims_min[0];
        const int ghostcell_dim_1_min = ghostcell_dims_min[1];
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min +
                            (k + offset_2_min)*ghostcell_dim_0_min*
                                ghostcell_dim_1_min;
                        
                        X[si][idx_min] = X[si][idx_min]/sum[idx_min];
                    }
                }
            }
        }
    }
    
    /*
     * Compute the effective binary diffusivity for each species.
     */
    
    if (domain.empty())
    {
        data_mass_diffusivities->fillAll(double(0));
    }
    else
    {
        data_mass_diffusivities->fillAll(double(0), domain);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_mass_diffusivities = offset_mass_diffusivities[0];
        const int offset_0_min = offset_min[0];
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int sj = 0; sj < d_num_species; sj++)
            {
                if (si != sj)
                {
                    int idx_ij = 0;
                    if (sj > si)
                    {
                        idx_ij = (d_num_species - 1)*d_num_species/2 -
                            (d_num_species - 1 - si)*(d_num_species - si)/2 +
                            (sj - (si + 1));
                    }
                    else
                    {
                        idx_ij = (d_num_species - 1)*d_num_species/2 -
                            (d_num_species - 1 - sj)*(d_num_species - sj)/2 +
                            (si - (sj + 1));
                    }
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_mass_diffusivities = i + offset_0_mass_diffusivities;
                        const int idx_min = i + offset_0_min;
                        
                        D[si][idx_mass_diffusivities] += (X[sj][idx_min] + double(EPSILON))/
                            (D_ij[idx_ij][idx_min] + double(EPSILON));
                    }
                }
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_mass_diffusivities = i + offset_0_mass_diffusivities;
                const int idx_min = i + offset_0_min;
                
                D[si][idx_mass_diffusivities] = (double(1) - X[si][idx_min] + double(EPSILON))/
                    D[si][idx_mass_diffusivities];
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        const int offset_0_mass_diffusivities = offset_mass_diffusivities[0];
        const int offset_1_mass_diffusivities = offset_mass_diffusivities[1];
        const int ghostcell_dim_0_mass_diffusivities = ghostcell_dims_mass_diffusivities[0];
        
        const int offset_0_min = offset_min[0];
        const int offset_1_min = offset_min[1];
        const int ghostcell_dim_0_min = ghostcell_dims_min[0];
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int sj = 0; sj < d_num_species; sj++)
            {
                if (si != sj)
                {
                    int idx_ij = 0;
                    if (sj > si)
                    {
                        idx_ij = (d_num_species - 1)*d_num_species/2 -
                            (d_num_species - 1 - si)*(d_num_species - si)/2 +
                            (sj - (si + 1));
                    }
                    else
                    {
                        idx_ij = (d_num_species - 1)*d_num_species/2 -
                            (d_num_species - 1 - sj)*(d_num_species - sj)/2 +
                            (si - (sj + 1));
                    }
                    
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_mass_diffusivities = (i + offset_0_mass_diffusivities) +
                                (j + offset_1_mass_diffusivities)*ghostcell_dim_0_mass_diffusivities;
                            
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min;
                            
                            D[si][idx_mass_diffusivities] += (X[sj][idx_min] + double(EPSILON))/
                                (D_ij[idx_ij][idx_min] + double(EPSILON));
                        }
                    }
                }
            }
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_mass_diffusivities = (i + offset_0_mass_diffusivities) +
                        (j + offset_1_mass_diffusivities)*ghostcell_dim_0_mass_diffusivities;
                    
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min;
                    
                    D[si][idx_mass_diffusivities] = (double(1) - X[si][idx_min] + double(EPSILON))/
                        D[si][idx_mass_diffusivities];
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        const int offset_0_mass_diffusivities = offset_mass_diffusivities[0];
        const int offset_1_mass_diffusivities = offset_mass_diffusivities[1];
        const int offset_2_mass_diffusivities = offset_mass_diffusivities[2];
        const int ghostcell_dim_0_mass_diffusivities = ghostcell_dims_mass_diffusivities[0];
        const int ghostcell_dim_1_mass_diffusivities = ghostcell_dims_mass_diffusivities[1];
        
        const int offset_0_min = offset_min[0];
        const int offset_1_min = offset_min[1];
        const int offset_2_min = offset_min[2];
        const int ghostcell_dim_0_min = ghostcell_dims_min[0];
        const int ghostcell_dim_1_min = ghostcell_dims_min[1];
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int sj = 0; sj < d_num_species; sj++)
            {
                if (si != sj)
                {
                    int idx_ij = 0;
                    if (sj > si)
                    {
                        idx_ij = (d_num_species - 1)*d_num_species/2 -
                            (d_num_species - 1 - si)*(d_num_species - si)/2 +
                            (sj - (si + 1));
                    }
                    else
                    {
                        idx_ij = (d_num_species - 1)*d_num_species/2 -
                            (d_num_species - 1 - sj)*(d_num_species - sj)/2 +
                            (si - (sj + 1));
                    }
                    
                    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                    {
                        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_mass_diffusivities = (i + offset_0_mass_diffusivities) +
                                    (j + offset_1_mass_diffusivities)*ghostcell_dim_0_mass_diffusivities +
                                    (k + offset_2_mass_diffusivities)*ghostcell_dim_0_mass_diffusivities*
                                        ghostcell_dim_1_mass_diffusivities;
                                
                                const int idx_min = (i + offset_0_min) +
                                    (j + offset_1_min)*ghostcell_dim_0_min +
                                    (k + offset_2_min)*ghostcell_dim_0_min*
                                        ghostcell_dim_1_min;
                                
                                D[si][idx_mass_diffusivities] += (X[sj][idx_min] + double(EPSILON))/
                                    (D_ij[idx_ij][idx_min] + double(EPSILON));
                            }
                        }
                    }
                }
            }
            
            for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_mass_diffusivities = (i + offset_0_mass_diffusivities) +
                            (j + offset_1_mass_diffusivities)*ghostcell_dim_0_mass_diffusivities +
                            (k + offset_2_mass_diffusivities)*ghostcell_dim_0_mass_diffusivities*
                                ghostcell_dim_1_mass_diffusivities;
                        
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min +
                            (k + offset_2_min)*ghostcell_dim_0_min*
                                ghostcell_dim_1_min;
                        
                        D[si][idx_mass_diffusivities] = (double(1) - X[si][idx_min] + double(EPSILON))/
                            D[si][idx_mass_diffusivities];
                    }
                }
            }
        }
    }
}


/*
 * Get the molecular properties of a species.
 */
void
EquationOfMassDiffusivityMixingRulesReid::getSpeciesMolecularProperties(
    std::vector<double*>& species_molecular_properties,
    const int species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >= 3);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    *(species_molecular_properties[0]) = d_species_epsilon_by_k[species_index];
    *(species_molecular_properties[1]) = d_species_sigma[species_index];
    *(species_molecular_properties[2]) = d_species_M[species_index];
}


/*
 * Compute the mass diffusivity of a binary mixture.
 */
double
EquationOfMassDiffusivityMixingRulesReid::getMassDiffusivity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties_1,
    const std::vector<const double*>& molecular_properties_2) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties_1.size()) >= 3);
    TBOX_ASSERT(static_cast<int>(molecular_properties_2.size()) >= 3);
#endif
    
    const double A = double(1.06036);
    const double B = double(-0.1561);
    const double C = double(0.19300);
    const double D = double(-0.47635);
    const double E = double(1.03587);
    const double F = double(-1.52996);
    const double G = double(1.76474);
    const double H = double(-3.89411);
    
    const double& epsilon_by_k_1 = *(molecular_properties_1[0]);
    const double& sigma_1 = *(molecular_properties_1[1]);
    const double& M_1 = *(molecular_properties_1[2]);
    
    const double& epsilon_by_k_2 = *(molecular_properties_2[0]);
    const double& sigma_2 = *(molecular_properties_2[1]);
    const double& M_2 = *(molecular_properties_2[2]);
    
    const double T_epsilon_12 = sqrt(epsilon_by_k_1*epsilon_by_k_2);
    const double sigma_12 = double(1)/double(2)*(sigma_1 + sigma_2);
    const double M_12 = double(2)/(double(1)/M_1 + double(1)/M_2);
    
    const double& p = *pressure;
    const double& T = *temperature;
    
    const double T_star_12 = T/T_epsilon_12;
    const double Omega_D_12 = A*pow(T_star_12, B) + C*exp(D*T_star_12) + E*exp(F*T_star_12) + G*exp(H*T_star_12);
    
    double D_12 = double(0.0266)*pow(T, double(3)/double(2))/(Omega_D_12*p*sqrt(M_12)*sigma_12*sigma_12);
    
    return D_12;
}
