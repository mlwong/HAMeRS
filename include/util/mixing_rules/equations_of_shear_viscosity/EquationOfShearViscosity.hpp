#ifndef EQUATION_OF_SHEAR_VISCOSITY_HPP
#define EQUATION_OF_SHEAR_VISCOSITY_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class EquationOfShearViscosity
{
    public:
        EquationOfShearViscosity(
            const std::string& object_name,
            const tbox::Dimension& dim):
                d_object_name(object_name),
                d_dim(dim)
        {}
        
        /*
         * Print all characteristics of the equation of shear viscosity class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Compute the shear viscosity.
         */
        virtual double
        getShearViscosity(
            const double* const pressure,
            const double* const temperature,
            const std::vector<const double*>& molecular_properties) const = 0;
        
        /*
         * Compute the shear viscosity.
         */
        void
        computeShearViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& molecular_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeShearViscosity(
                data_shear_viscosity,
                data_pressure,
                data_temperature,
                molecular_properties,
                empty_box);
        }
        
        /*
         * Compute the shear viscosity.
         */
        virtual void
        computeShearViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const std::vector<const double*>& molecular_properties,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the shear viscosity.
         */
        void
        computeShearViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_molecular_properties) const
        {
            const hier::Box empty_box(d_dim);
            computeShearViscosity(
                data_shear_viscosity,
                data_pressure,
                data_temperature,
                data_molecular_properties,
                empty_box);
        }
        
        /*
         * Compute the shear viscosity.
         */
        virtual void
        computeShearViscosity(
            boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
            const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
            const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
            const boost::shared_ptr<pdat::CellData<double> >& data_molecular_properties,
            const hier::Box& domain) const = 0;
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
};

#endif /* EQUATION_OF_SHEAR_VISCOSITY_HPP */
