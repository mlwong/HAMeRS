#ifndef GRADIENT_SENSOR_DUCROS_HPP
#define GRADIENT_SENSOR_DUCROS_HPP

#include "util/gradient_sensors/GradientSensor.hpp"

class GradientSensorDucros: public GradientSensor
{
    public:
        /*
         * If useStrainRateInsteadOfDilatation is true, the idea of the sensor is simiar to
         * Travin, A., Shur, M., Strelets, M., & Spalart, P. R. (2002) of the paper
         * "Physical and numerical upgrades in the detached-eddy simulation of complex turbulent flows."
         * In Advances in LES of Complex Flows: Proceedings of the Euromech Colloquium 412, Springer Netherlands.
         */
        GradientSensorDucros(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const bool use_dilatation,
            const bool use_strain_rate);
        
        ~GradientSensorDucros() {}
        
        /*
         * Compute the gradient with the given cell data.
         */
        void
        computeGradient(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& gradient,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& cell_data_velocity,
            hier::Patch& patch,
            const int depth = 0);
        
    private:
        const bool d_use_dilatation;
        const bool d_use_strain_rate;
};

#endif /* GRADIENT_SENSOR_DUCROS_HPP */