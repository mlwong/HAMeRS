#ifndef GRADIENT_SENSOR_HPP
#define GRADIENT_SENSOR_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class GradientSensor
{
    public:
        GradientSensor(
            const std::string& object_name,
            const tbox::Dimension& dim):
                d_object_name(object_name),
                d_dim(dim),
                d_num_gradient_ghosts(hier::IntVector::getZero(d_dim)),
                d_num_ghosts(hier::IntVector::getZero(d_dim))
        {}
        
        virtual ~GradientSensor() {}
        
        /*
         * Get the number of ghost cells needed by the gradient sensor.
         */
        hier::IntVector
        getGradientSensorNumberOfGhostCells(void) const
        {
            return d_num_gradient_ghosts;
        }
        
        /*
         * Set the number of ghost cells needed.
         */
        void
        setNumberOfGhostCells(const hier::IntVector& num_ghosts)
        {
            d_num_ghosts = num_ghosts;
        }
        
        /*
         * Compute the gradient.
         */
        virtual void
        computeGradient(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<double> > cell_data,
            boost::shared_ptr<pdat::CellData<double> > gradient,
            int depth = 0) = 0;
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * Number of ghost cells needed by the gradient sensor.
         */
        hier::IntVector d_num_gradient_ghosts;
        
        /*
         * Number of ghost cells for time-independent variables.
         */
        hier::IntVector d_num_ghosts;
        
};

#endif /* GRADIENT_SENSOR_HPP */
