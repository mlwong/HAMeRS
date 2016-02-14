#ifndef GRADIENT_SENSOR_JAMESON_HPP
#define GRADIENT_SENSOR_JAMESON_HPP

#include "utils/gradient_sensors/GradientSensor.hpp"

class GradientSensorJameson: public GradientSensor
{
    public:
        GradientSensorJameson(
            const std::string& object_name,
            const tbox::Dimension& dim);
        
        /*
         * Compute the gradient with the given cell data.
         */
        boost::shared_ptr<pdat::CellData<double> >
        ComputeGradient(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data);
        
};

#endif /* GRADIENT_SENSOR_JAMESON_HPP */