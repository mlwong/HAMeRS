#include "util/mixing_rules/equations_of_mass_diffusivity/constant/EquationOfMassDiffusivityMixingRulesConstant.hpp"

#define EPSILON HAMERS_EPSILON

EquationOfMassDiffusivityMixingRulesConstant::EquationOfMassDiffusivityMixingRulesConstant(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const HAMERS_SHARED_PTR<tbox::Database>& equation_of_mass_diffusivity_mixing_rules_db):
        EquationOfMassDiffusivityMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_mass_diffusivity_mixing_rules_db)
{
    if (d_num_species > 2)
    {
        /*
         * Get the mass diffusivity of each pair of species from the database.
         */
        
        if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("D_ij"))
        {
            size_t D_ij_array_size =
                equation_of_mass_diffusivity_mixing_rules_db->getArraySize("D_ij");
            if (static_cast<int>(D_ij_array_size) == (d_num_species - 1)*d_num_species/2)
            {
                d_D_ij =
                    equation_of_mass_diffusivity_mixing_rules_db->getRealVector("D_ij");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "number of 'D_ij' entries must be equal to '(num_species - 1)*(num_species)/2'."
                    << std::endl);
            }
        }
        else if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("d_D_ij"))
        {
            size_t species_D_ij_array_size =
                equation_of_mass_diffusivity_mixing_rules_db->getArraySize("d_D_ij");
            if (static_cast<int>(species_D_ij_array_size) == (d_num_species - 1)*d_num_species/2)
            {
                d_D_ij =
                    equation_of_mass_diffusivity_mixing_rules_db->getRealVector("d_D_ij");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "number of 'd_D_ij' entries must be equal to '(num_species - 1)*(num_species)/2'."
                    << std::endl);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Key data 'D_ij'/'d_D_ij'"
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
                    equation_of_mass_diffusivity_mixing_rules_db->getRealVector("species_M");
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
                    equation_of_mass_diffusivity_mixing_rules_db->getRealVector("d_species_M");
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
    else
    {
        /*
         * Get the mass diffusivity of each species from the database.
         */
        
        if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("species_D"))
        {
            size_t species_D_array_size =
                equation_of_mass_diffusivity_mixing_rules_db->getArraySize("species_D");
            if (static_cast<int>(species_D_array_size) == d_num_species)
            {
                d_species_D =
                    equation_of_mass_diffusivity_mixing_rules_db->getRealVector("species_D");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "number of 'species_D' entries must be equal to 'num_species'."
                    << std::endl);
            }
        }
        else if (equation_of_mass_diffusivity_mixing_rules_db->keyExists("d_species_D"))
        {
            size_t species_D_array_size =
                equation_of_mass_diffusivity_mixing_rules_db->getArraySize("d_species_D");
            if (static_cast<int>(species_D_array_size) == d_num_species)
            {
                d_species_D =
                    equation_of_mass_diffusivity_mixing_rules_db->getRealVector("d_species_D");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "number of 'd_species_D' entries must be equal to 'd_num_species'."
                    << std::endl);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Key data 'species_D'/'d_species_D'"
                << "not found in data for equation of mass diffusivity mixing rules."
                << std::endl);
        }
    }
}


/*
 * Print all characteristics of the equation of mass diffusivity class.
 */
void
EquationOfMassDiffusivityMixingRulesConstant::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfMassDiffusivityMixingRulesConstant object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfMassDiffusivityMixingRulesConstant: this = "
       << (EquationOfMassDiffusivityMixingRulesConstant *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    if (d_num_species > 2)
    {
        /*
         * Print the mass diffusivity of each pair of species.
         */
        
        os << "d_D_ij = ";
        for (int si = 0; si < (d_num_species - 1)*d_num_species/2 - 1; si++)
        {
            os << d_D_ij[si] << ", ";
            
        }
        os << d_D_ij[(d_num_species - 1)*d_num_species/2 - 1];
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
    else
    {
        /*
         * Print the mass diffusivity of each species.
         */
        
        os << "d_species_D = ";
        for (int si = 0; si < d_num_species - 1; si++)
        {
            os << d_species_D[si] << ", ";
        }
        os << d_species_D[d_num_species - 1];
        os << std::endl;
    }
}


/*
 * Put the characteristics of the equation of mass diffusivity mixing rules class into the restart
 * database.
 */
void
EquationOfMassDiffusivityMixingRulesConstant::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    if (d_num_species > 2)
    {
        restart_db->putRealVector("d_D_ij", d_D_ij);
        restart_db->putRealVector("d_species_M", d_species_M);
    }
    else
    {
        restart_db->putRealVector("d_species_D", d_species_D);
    }
}


/*
 * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfMassDiffusivityMixingRulesConstant::getMassDiffusivities(
    std::vector<Real*>& mass_diffusivities,
    const Real* const pressure,
    const Real* const temperature,
    const std::vector<const Real*>& mass_fractions) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(static_cast<int>(mass_diffusivities.size()) == d_num_species);
    TBOX_ASSERT((static_cast<int>(mass_fractions.size()) == d_num_species) ||
                (static_cast<int>(mass_fractions.size()) == d_num_species - 1));
#endif
    
    if (d_num_species > 2)
    {
        /*
         * Initialize the containers and pointers to the containers for the molecular properties
         * of species i and j.
         */
        
        std::vector<Real> species_molecular_properties_i;
        std::vector<Real*> species_molecular_properties_ptr_i;
        std::vector<const Real*> species_molecular_properties_const_ptr_i;
        
        const int num_molecular_properties = getNumberOfSpeciesMolecularProperties();
        
        species_molecular_properties_i.resize(num_molecular_properties);
        species_molecular_properties_ptr_i.reserve(num_molecular_properties);
        species_molecular_properties_const_ptr_i.reserve(num_molecular_properties);
        
        for (int mi = 0; mi < num_molecular_properties; mi++)
        {
            species_molecular_properties_ptr_i.push_back(&species_molecular_properties_i[mi]);
            species_molecular_properties_const_ptr_i.push_back(&species_molecular_properties_i[mi]);
        }
        
        /*
         * Compute the mole fractions.
         */
        
        std::vector<Real> X;
        X.reserve(d_num_species);
        
        Real sum = Real(0);
        
        if (static_cast<int>(mass_fractions.size()) == d_num_species)
        {
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr_i, si);
                X.push_back((*(mass_fractions[si]))/(species_molecular_properties_i[0]));
                X[si] = std::max(Real(0), X[si]);
                sum += X[si];
            }
        }
        else if (static_cast<int>(mass_fractions.size()) == d_num_species - 1)
        {
            Real Y_last = Real(1);
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr_i, si);
                X.push_back((*(mass_fractions[si]))/(species_molecular_properties_i[0]));
                X[si] = std::max(Real(0), X[si]);
                sum += X[si];
                
                // Compute the mass fraction of the last species.
                Y_last -= *(mass_fractions[si]);
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr_i, d_num_species - 1);
            X.push_back(Y_last/(species_molecular_properties_i[0]));
            X[d_num_species - 1] = std::max(Real(0), X[d_num_species - 1]);
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
            Real& D = *(mass_diffusivities[si]);
            D = Real(0);
            
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
                    
                    D += (X[sj] + Real(EPSILON))/(d_D_ij[idx] + Real(EPSILON));
                }
            }
            
            D = (Real(1) - X[si] + Real(EPSILON))/D;
        }
    }
    else if (d_num_species == 2)
    {
        getMassDiffusivitiesForTwoSpecies(
            mass_diffusivities,
            pressure,
            temperature,
            mass_fractions);
    }
    else // d_num_species == 1
    {
        *(mass_diffusivities[0]) = Real(0);
    }
}


/*
 * Compute the mass diffusivities of the mixture with isothermal and isobaric equilibrium assumptions.
 */
void
EquationOfMassDiffusivityMixingRulesConstant::computeMassDiffusivities(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_diffusivities,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_pressure);
    NULL_USE(data_temperature);
    
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
            data_mass_diffusivities->fillAll(Real(0));
        }
        else
        {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(data_mass_diffusivities->getGhostBox().contains(domain));
#endif
            
            data_mass_diffusivities->fillAll(Real(0), domain);
        }
        
        return;
    }
    
    if (d_num_species > 2)
    {
        // Get the dimensions of the ghost cell boxes.
        const hier::Box ghost_box_mass_diffusivities = data_mass_diffusivities->getGhostBox();
        const hier::IntVector ghostcell_dims_mass_diffusivities = ghost_box_mass_diffusivities.numberCells();
        
        const hier::Box ghost_box_mass_fractions = data_mass_fractions->getGhostBox();
        const hier::IntVector ghostcell_dims_mass_fractions = ghost_box_mass_fractions.numberCells();
        
        // Delcare data containers for mole fractions.
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mole_fractions;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_sum;
        
        // Declare data container for last volume fraction.
        HAMERS_SHARED_PTR<pdat::CellData<Real> > data_mass_fractions_last;
        
        /*
         * Get the local lower index and number of cells in each direction of the domain.
         * Also, get the offsets of all data and dimensions of the ghost cell box for binary
         * mass diffusivities, mole fractions and last mass fraction and allocate memory.
         */
        
        hier::IntVector domain_lo(d_dim);
        hier::IntVector domain_dims(d_dim);
        
        hier::IntVector offset_mass_diffusivities(d_dim);
        hier::IntVector offset_mass_fractions(d_dim);
        hier::IntVector offset_min(d_dim);
        
        hier::IntVector ghostcell_dims_min(d_dim);
        
        if (domain.empty())
        {
            // Get the numbers of ghost cells.
            const hier::IntVector num_ghosts_mass_diffusivities = data_mass_diffusivities->getGhostCellWidth();
            const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
            
            // Get the dimensions of box that covers the interior of patch.
            const hier::Box interior_box = data_mass_diffusivities->getBox();
            const hier::IntVector interior_dims = interior_box.numberCells();
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(data_mass_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
            
            /*
             * Get the minimum number of ghost cells and the dimensions of the ghost cell box for binary
             * mass diffusivities, mole fractions and last mass fraction.
             */
            
            hier::IntVector num_ghosts_min(d_dim);
            
            num_ghosts_min = num_ghosts_mass_diffusivities;
            num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
            
            hier::Box ghost_box = interior_box;
            ghost_box.grow(num_ghosts_min);
            
            domain_lo = -num_ghosts_min;
            domain_dims = ghost_box.numberCells();
            
            offset_min = num_ghosts_min;
            offset_mass_diffusivities = num_ghosts_mass_diffusivities;
            offset_mass_fractions = num_ghosts_mass_fractions;
            
            ghostcell_dims_min = interior_dims + num_ghosts_min*2;
            
            data_mole_fractions = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
                interior_box, d_num_species, num_ghosts_min);
            
            data_sum = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
                interior_box, 1, num_ghosts_min);
            
            if (data_mass_fractions->getDepth() == d_num_species - 1)
            {
                data_mass_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
                    interior_box, 1, num_ghosts_min);
            }
        }
        else
        {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
            TBOX_ASSERT(data_mass_diffusivities->getGhostBox().contains(domain));
            TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
            
            domain_lo = hier::IntVector::getZero(d_dim);
            domain_dims = domain.numberCells();
            
            offset_min = hier::IntVector::getZero(d_dim);
            offset_mass_diffusivities = domain.lower() - ghost_box_mass_diffusivities.lower();
            offset_mass_fractions = domain.lower() - ghost_box_mass_fractions.lower();
            
            ghostcell_dims_min = domain_dims;
            
            data_mole_fractions = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
                domain, d_num_species, hier::IntVector::getZero(d_dim));
            
            data_sum = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
                domain, 1, hier::IntVector::getZero(d_dim));
            
            if (data_mass_fractions->getDepth() == d_num_species - 1)
            {
                data_mass_fractions_last = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
                    domain, 1, hier::IntVector::getZero(d_dim));
            }
        }
        
        data_sum->fillAll(Real(0));
        
        /*
         * Get the pointers to the cell data.
         */
        
        std::vector<Real*> D;
        D.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            D.push_back(data_mass_diffusivities->getPointer(si));
        }
        std::vector<Real*> X;
        X.reserve(d_num_species);
        for (int si = 0; si < d_num_species; si++)
        {
            X.push_back(data_mole_fractions->getPointer(si));
        }
        Real* sum = data_sum->getPointer(0);
        
        /*
         * Compute the mole fractions.
         */
        
        if (data_mass_fractions->getDepth() == d_num_species)
        {
            /*
             * Get the pointers to the cell data of mass fractions.
             */
            
            std::vector<Real*> Y;
            Y.reserve(d_num_species);
            for (int si = 0; si < d_num_species; si++)
            {
                Y.push_back(data_mass_fractions->getPointer(si));
            }
            
            /*
             * Initialize the container and pointer to the container for the molecular properties.
             */
            
            std::vector<Real> species_molecular_properties;
            std::vector<Real*> species_molecular_properties_ptr;
            
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
                    
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = i + offset_0_min;
                        const int idx_mass_fractions = i + offset_0_mass_fractions;
                        
                        X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[0];
                        X[si][idx_min] = std::max(Real(0), X[si][idx_min]);
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
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min;
                            
                            const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                                (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[0];
                            X[si][idx_min] = std::max(Real(0), X[si][idx_min]);
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
                            HAMERS_PRAGMA_SIMD
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
                                
                                X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[0];
                                X[si][idx_min] = std::max(Real(0), X[si][idx_min]);
                                sum[idx_min] += X[si][idx_min];
                            }
                        }
                    }
                }
            }
        }
        else if (data_mass_fractions->getDepth() == d_num_species - 1)
        {
            data_mass_fractions_last->fillAll(Real(1));
            
            /*
             * Get the pointers to the cell data of mass fractions.
             */
            
            std::vector<Real*> Y;
            Y.reserve(d_num_species - 1);
            for (int si = 0; si < d_num_species - 1; si++)
            {
                Y.push_back(data_mass_fractions->getPointer(si));
            }
            
            Real* Y_last = data_mass_fractions_last->getPointer(0);
            
            /*
             * Initialize the container and pointer to the container for the molecular properties.
             */
            
            std::vector<Real> species_molecular_properties;
            std::vector<Real*> species_molecular_properties_ptr;
            
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
                    
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_min = i + offset_0_min;
                        const int idx_mass_fractions = i + offset_0_mass_fractions;
                        
                        X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[0];
                        X[si][idx_min] = std::max(Real(0), X[si][idx_min]);
                        sum[idx_min] += X[si][idx_min];
                        
                        // Compute the mass fraction of the last species.
                        Y_last[idx_min] -= Y[si][idx_mass_fractions];
                    }
                }
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_min = i + offset_0_min;
                    
                    X[d_num_species - 1][idx_min] = Y_last[idx_min]/species_molecular_properties[0];
                    X[d_num_species - 1][idx_min] = std::max(Real(0), X[d_num_species - 1][idx_min]);
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
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min;
                            
                            const int idx_mass_fractions = (i + offset_0_mass_fractions) +
                                (j + offset_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[0];
                            X[si][idx_min] = std::max(Real(0), X[si][idx_min]);
                            sum[idx_min] += X[si][idx_min];
                            
                            // Compute the mass fraction of the last species.
                            Y_last[idx_min] -= Y[si][idx_mass_fractions];
                        }
                    }
                }
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr, d_num_species - 1);
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min;
                        
                        X[d_num_species - 1][idx_min] = Y_last[idx_min]/species_molecular_properties[0];
                        X[d_num_species - 1][idx_min] = std::max(Real(0), X[d_num_species - 1][idx_min]);
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
                            HAMERS_PRAGMA_SIMD
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
                                
                                X[si][idx_min] = Y[si][idx_mass_fractions]/species_molecular_properties[0];
                                X[si][idx_min] = std::max(Real(0), X[si][idx_min]);
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
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_min = (i + offset_0_min) +
                                (j + offset_1_min)*ghostcell_dim_0_min +
                                (k + offset_2_min)*ghostcell_dim_0_min*
                                    ghostcell_dim_1_min;
                            
                            X[d_num_species - 1][idx_min] = Y_last[idx_min]/species_molecular_properties[0];
                            X[d_num_species - 1][idx_min] = std::max(Real(0), X[d_num_species - 1][idx_min]);
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
                HAMERS_PRAGMA_SIMD
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
                    HAMERS_PRAGMA_SIMD
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
                        HAMERS_PRAGMA_SIMD
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
            data_mass_diffusivities->fillAll(Real(0));
        }
        else
        {
            data_mass_diffusivities->fillAll(Real(0), domain);
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
                        
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_mass_diffusivities = i + offset_0_mass_diffusivities;
                            const int idx_min = i + offset_0_min;
                            
                            D[si][idx_mass_diffusivities] += (X[sj][idx_min] + Real(EPSILON))/
                                (d_D_ij[idx_ij]+ Real(EPSILON));
                        }
                    }
                }
                
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_mass_diffusivities = i + offset_0_mass_diffusivities;
                    const int idx_min = i + offset_0_min;
                    
                    D[si][idx_mass_diffusivities] = (Real(1) - X[si][idx_min] + Real(EPSILON))/
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
                            HAMERS_PRAGMA_SIMD
                            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_mass_diffusivities = (i + offset_0_mass_diffusivities) +
                                    (j + offset_1_mass_diffusivities)*ghostcell_dim_0_mass_diffusivities;
                                
                                const int idx_min = (i + offset_0_min) +
                                    (j + offset_1_min)*ghostcell_dim_0_min;
                                
                                D[si][idx_mass_diffusivities] += (X[sj][idx_min] + Real(EPSILON))/
                                    (d_D_ij[idx_ij]+ Real(EPSILON));
                            }
                        }
                    }
                }
                
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_mass_diffusivities = (i + offset_0_mass_diffusivities) +
                            (j + offset_1_mass_diffusivities)*ghostcell_dim_0_mass_diffusivities;
                        
                        const int idx_min = (i + offset_0_min) +
                            (j + offset_1_min)*ghostcell_dim_0_min;
                        
                        D[si][idx_mass_diffusivities] = (Real(1) - X[si][idx_min] + Real(EPSILON))/
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
                                HAMERS_PRAGMA_SIMD
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
                                    
                                    D[si][idx_mass_diffusivities] += (X[sj][idx_min] + Real(EPSILON))/
                                        (d_D_ij[idx_ij]+ Real(EPSILON));
                                }
                            }
                        }
                    }
                }
                
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            D[si][idx_mass_diffusivities] = (Real(1) - X[si][idx_min] + Real(EPSILON))/
                                D[si][idx_mass_diffusivities];
                        }
                    }
                }
            }
        }
        
        
        
    }
    else if (d_num_species == 2)
    {
        computeMassDiffusivitiesForTwoSpecies(
            data_mass_diffusivities,
            data_pressure,
            data_temperature,
            data_mass_fractions,
            domain);
    }
}


/*
 * Get the molecular properties of a species.
 */
void
EquationOfMassDiffusivityMixingRulesConstant::getSpeciesMolecularProperties(
    std::vector<Real*>& species_molecular_properties,
    const int species_index) const
{
    if (d_num_species > 2)
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >= 1);
        TBOX_ASSERT(species_index >= 0);
        TBOX_ASSERT(species_index < d_num_species);
#endif
        *(species_molecular_properties[0]) = d_species_M[species_index];
    }
    else
    {
        getSpeciesMolecularPropertiesForTwoSpecies(
            species_molecular_properties,
            species_index);
    }
}


void
EquationOfMassDiffusivityMixingRulesConstant::getMassDiffusivitiesForTwoSpecies(
    std::vector<Real*>& mass_diffusivities,
    const Real* const pressure,
    const Real* const temperature,
    const std::vector<const Real*>& mass_fractions) const
{
    NULL_USE(pressure);
    NULL_USE(temperature);
    NULL_USE(mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(static_cast<int>(mass_diffusivities.size()) == d_num_species);
#endif
    
    for (int si = 0; si < d_num_species; si++)
    {
        *(mass_diffusivities[si]) = d_species_D[si];
    }
}


void
EquationOfMassDiffusivityMixingRulesConstant::computeMassDiffusivitiesForTwoSpecies(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_diffusivities,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_mass_fractions,
    const hier::Box& domain) const
{
    NULL_USE(data_pressure);
    NULL_USE(data_temperature);
    NULL_USE(data_mass_fractions);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    
    TBOX_ASSERT(data_mass_diffusivities->getDepth() == d_num_species);
#endif
    
    if (domain.empty())
    {
        for (int si = 0; si < d_num_species; si++)
        {
            data_mass_diffusivities->fill(d_species_D[si], si);
        }
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mass_diffusivities->getGhostBox().contains(domain));
#endif
        for (int si = 0; si < d_num_species; si++)
        {
            data_mass_diffusivities->fill(d_species_D[si], domain, si);
        }
    }
}


void
EquationOfMassDiffusivityMixingRulesConstant::getSpeciesMolecularPropertiesForTwoSpecies(
    std::vector<Real*>& species_molecular_properties,
    const int species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(species_molecular_properties.size()) >= 1);
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    *(species_molecular_properties[0]) = d_species_D[species_index];
}
