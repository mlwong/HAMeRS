#ifndef GRADIENT_SENSOR_SECOND_DERIVATIVE_HPP
#define GRADIENT_SENSOR_SECOND_DERIVATIVE_HPP

#include "util/gradient_sensors/GradientSensor.hpp"

class GradientSensorSecondDerivative: public GradientSensor
{
    public:
        GradientSensorSecondDerivative(
            const std::string& object_name,
            const tbox::Dimension& dim);
        
        ~GradientSensorSecondDerivative() {}
        
        /*
         * Compute the gradient with the given cell data.
         */
        void
        computeGradient(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data,
            boost::shared_ptr<pdat::CellData<double> > gradient,
            int depth = 0);
        
        /*
         * Compute the gradient and the local mean of the given cell data.
         */
        void
        computeGradientWithVariableLocalMean(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data,
            boost::shared_ptr<pdat::CellData<double> > gradient,
            boost::shared_ptr<pdat::CellData<double> > variable_local_mean,
            int depth = 0);
        
};

#endif /* GRADIENT_SENSOR_SECOND_DERIVATIVE_HPP */
