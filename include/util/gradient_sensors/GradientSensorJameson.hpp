#ifndef GRADIENT_SENSOR_JAMESON_HPP
#define GRADIENT_SENSOR_JAMESON_HPP

#include "util/gradient_sensors/GradientSensor.hpp"

class GradientSensorJameson: public GradientSensor
{
    public:
        GradientSensorJameson(
            const std::string& object_name,
            const tbox::Dimension& dim);
        
        ~GradientSensorJameson() {}
        
        /*
         * Compute the gradient with the given cell data.
         */
        void
        computeGradient(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& gradient,
            const HAMERS_SHARED_PTR<pdat::CellData<double> >& cell_data,
            hier::Patch& patch,
            const int depth = 0);
        
};

#endif /* GRADIENT_SENSOR_JAMESON_HPP */
