#include "util/mixing_rules/equations_of_bulk_viscosity/Cramer/EquationOfBulkViscosityCramer.hpp"

#include <cmath>

/*
 * Print all characteristics of the equation of bulk viscosity class.
 */
void
EquationOfBulkViscosityCramer::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfBulkViscosityCramer object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfBulkViscosityCramer: this = "
       << (EquationOfBulkViscosityCramer *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the bulk viscosity.
 */
double
EquationOfBulkViscosityCramer::getBulkViscosity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 7);
#endif
    
    const double& T = *temperature;
    
    const double& gamma = *(molecular_properties[0]);
    
    const double& A_r = *(molecular_properties[1]);
    const double& B_r = *(molecular_properties[2]);
    
    const double& c_v_v = *(molecular_properties[3]);
    const double& A_v   = *(molecular_properties[4]);
    const double& B_v   = *(molecular_properties[5]);
    const double& C_v   = *(molecular_properties[6]);
    
    double mu_v = A_r + B_r*T +
        (gamma - 1.0)*(gamma - 1.0)*c_v_v*A_v*exp(B_v/(pow(T, 1.0/3.0)) + C_v/(pow(T, 2.0/3.0)));
    
    return mu_v;
}


/*
 * Compute the bulk viscosity.
 */
void
EquationOfBulkViscosityCramer::computeBulkViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& molecular_properties,
    const hier::Box& domain) const
{
}


/*
 * Compute the bulk viscosity.
 */
void
EquationOfBulkViscosityCramer::computeBulkViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_molecular_properties,
    const hier::Box& domain) const
{
}
