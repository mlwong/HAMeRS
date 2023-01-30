#ifndef IMMERSED_BOUNDARIES_HPP
#define IMMERSED_BOUNDARIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class ImmersedBoundaries
{
    public:
        ImmersedBoundaries(
            const std::string& object_name,
            const tbox::Dimension& dim):
                d_object_name(object_name),
                d_dim(dim)
        {}
        
        virtual ~ImmersedBoundaries() {}
        
        void setImmersedBoundaryMask();
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
};

#endif /* IMMERSED_BOUNDARIES_HPP */
