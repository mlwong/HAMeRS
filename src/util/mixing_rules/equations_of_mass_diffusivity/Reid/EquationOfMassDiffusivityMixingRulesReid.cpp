#include "util/mixing_rules/equations_of_mass_diffusivity/Reid/EquationOfMassDiffusivityMixingRulesReid.hpp"

#include <cmath>

#define EPSILON 1e-40

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
 * Compute the mass diffusivities of the mixture with isothermal and isobaric assumptions.
 */
void
EquationOfMassDiffusivityMixingRulesReid::getMassDiffusivities(
    std::vector<double*>& mass_diffusivities,
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& mass_fraction) const
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(mass_diffusivities.size()) == d_num_species);
    TBOX_ASSERT((d_mixing_closure_model == MIXING_CLOSURE_MODEL::ISOTHERMAL_AND_ISOBARIC) ||
                (d_mixing_closure_model == MIXING_CLOSURE_MODEL::NO_MODEL && d_num_species == 1));
    TBOX_ASSERT((static_cast<int>(mass_fraction.size()) == d_num_species) ||
                (static_cast<int>(mass_fraction.size()) == d_num_species - 1));
#endif
    
    if (d_num_species > 1)
    {
        /*
         * Initialize the containers and pointers to the containers for the molecular properties
         * of species 1 and 2.
         */
        
        std::vector<double> species_molecular_properties_1;
        std::vector<double> species_molecular_properties_2;
        std::vector<double*> species_molecular_properties_ptr_1;
        std::vector<double*> species_molecular_properties_ptr_2;
        std::vector<const double*> species_molecular_properties_const_ptr_1;
        std::vector<const double*> species_molecular_properties_const_ptr_2;
        
        const int num_molecular_properties = getNumberOfSpeciesMolecularProperties();
        
        species_molecular_properties_1.resize(num_molecular_properties);
        species_molecular_properties_2.resize(num_molecular_properties);
        species_molecular_properties_ptr_1.reserve(num_molecular_properties);
        species_molecular_properties_ptr_2.reserve(num_molecular_properties);
        species_molecular_properties_const_ptr_1.reserve(num_molecular_properties);
        species_molecular_properties_const_ptr_2.reserve(num_molecular_properties);
        
        for (int mi = 0; mi < num_molecular_properties; mi++)
        {
            species_molecular_properties_ptr_1.push_back(&species_molecular_properties_1[mi]);
            species_molecular_properties_ptr_2.push_back(&species_molecular_properties_2[mi]);
            species_molecular_properties_const_ptr_1.push_back(&species_molecular_properties_1[mi]);
            species_molecular_properties_const_ptr_2.push_back(&species_molecular_properties_2[mi]);
        }
        
        /*
         * Get the mass diffusion coefficient for each pair of species.
         */
        
        std::vector<double> D_ij(0.5*(d_num_species - 1)*d_num_species);
        
        for (int i = 0; i < d_num_species; i++)
        {
            for (int j = i + 1; j < d_num_species; j++)
            {
                const int idx = 0.5*(d_num_species - 1)*d_num_species -
                    0.5*(d_num_species - 1 - i)*(d_num_species - i) +
                    (j - (i + 1));
                
                getSpeciesMolecularProperties(species_molecular_properties_ptr_1, i);
                getSpeciesMolecularProperties(species_molecular_properties_ptr_2, j);
                
                D_ij[idx] = getMassDiffusivity(
                        pressure,
                        temperature,
                        species_molecular_properties_const_ptr_1,
                        species_molecular_properties_const_ptr_2);
            }
        }
        
        /*
         * Compute the mole fractions.
         */
        
        std::vector<double> X;
        X.reserve(d_num_species);
        
        double sum = 0.0;
        
        if (static_cast<int>(mass_fraction.size()) == d_num_species - 1)
        {
            double Y_last = 1.0;
            
            for (int si = 0; si < d_num_species - 1; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr_1, si);
                X.push_back((*(mass_fraction[si]))/(species_molecular_properties_1[2]));
                sum += X[si];
                
                // Compute the mass fraction of the last species.
                Y_last -= *(mass_fraction[si]);
            }
            
            getSpeciesMolecularProperties(species_molecular_properties_ptr_1, d_num_species - 1);
            X.push_back(Y_last/(species_molecular_properties_1[2]));
            sum += X[d_num_species - 1];
        }
        else
        {
            for (int si = 0; si < d_num_species; si++)
            {
                getSpeciesMolecularProperties(species_molecular_properties_ptr_1, si);
                X.push_back((*(mass_fraction[si]))/(species_molecular_properties_1[2]));
                X[si] = std::max(0.0, X[si]);
                sum += X[si];
            }
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            X[si] = X[si]/sum;
        }
        
        /*
         * Compute the effective binary diffusion coefficients for each species.
         */
        
        for (int si = 0; si < d_num_species; si++)
        {
            double& D = *(mass_diffusivities[si]);
            D = 0.0;
            
            for (int sj = 0; sj < d_num_species; sj++)
            {
                if (si != sj)
                {
                    int idx = 0.0;
                    if (sj > si)
                    {
                        idx = 0.5*(d_num_species - 1)*d_num_species -
                            0.5*(d_num_species - 1 - si)*(d_num_species - si) +
                            (sj - (si + 1));
                    }
                    else
                    {
                        idx = 0.5*(d_num_species - 1)*d_num_species -
                            0.5*(d_num_species - 1 - sj)*(d_num_species - sj) +
                            (si - (sj + 1));
                    }
                    
                    D += (X[sj] + EPSILON)/(D_ij[idx] + EPSILON);
                }
            }
            
            D = (1.0 - X[si] + EPSILON)/D;
        }
    }
    else
    {
        *(mass_diffusivities[0]) = 0.0;
    }
}


/*
 * Get the molecular properties of a species.
 */
void
EquationOfMassDiffusivityMixingRulesReid::getSpeciesMolecularProperties(
    std::vector<double*>& species_molecular_properties,
    const int& species_index) const
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
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
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties_1.size()) >= 3);
    TBOX_ASSERT(static_cast<int>(molecular_properties_2.size()) >= 3);
#endif
    
    double D_12 = 0.0;
    
    const double& epsilon_by_k_1 = *(molecular_properties_1[0]);
    const double& sigma_1 = *(molecular_properties_1[1]);
    const double& M_1 = *(molecular_properties_1[2]);
    
    const double& epsilon_by_k_2 = *(molecular_properties_2[0]);
    const double& sigma_2 = *(molecular_properties_2[1]);
    const double& M_2 = *(molecular_properties_2[2]);
    
    const double& T = *temperature;
    const double& p = *pressure;
    
    const double A = 1.06036;
    const double B = -0.1561;
    const double C = 0.19300;
    const double D = -0.47635;
    const double E = 1.03587;
    const double F = -1.52996;
    const double G = 1.76474;
    const double H = -3.89411;
    
    const double T_epsilon_12 = sqrt(epsilon_by_k_1*epsilon_by_k_2);
    const double T_star_12 = T/T_epsilon_12;
    const double Omega_D_12 = A*pow(T_star_12, B) + C*exp(D*T_star_12) + E*exp(F*T_star_12) + G*exp(H*T_star_12);
    
    const double M_12 = 2.0/(1.0/M_1 + 1.0/M_2);
    const double sigma_12 = 0.5*(sigma_1 + sigma_2);
    
    D_12 = 0.0266*pow(T, 1.5)/(Omega_D_12*p*sqrt(M_12)*sigma_12*sigma_12);
    
    return D_12;
}
