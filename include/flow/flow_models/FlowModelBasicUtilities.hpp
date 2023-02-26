#ifndef FLOW_MODEL_BASIC_UTILITIES_HPP
#define FLOW_MODEL_BASIC_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include <string>

namespace AVERAGING
{
    enum TYPE { SIMPLE,
                ROE };
}

class FlowModel;

class FlowModelBasicUtilities
{
    public:
        FlowModelBasicUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species),
                d_num_eqn(num_eqn),
                d_proj_var_conservative_averaging_type(AVERAGING::SIMPLE),
                d_proj_var_primitive_averaging_type(AVERAGING::SIMPLE)
        {}
        
        virtual ~FlowModelBasicUtilities() {}
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const HAMERS_WEAK_PTR<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Convert conservative variables to primitive variables.
         */
        virtual void
        convertConservativeVariablesToPrimitiveVariables(
            const std::vector<const Real*>& conservative_variables,
            const std::vector<Real*>& primitive_variables) = 0;
        
        /*
         * Convert conservative variables to primitive variables.
         */
        virtual void
        convertConservativeVariablesToPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables) = 0;
        
        /*
         * Convert primitive variables to conservative variables.
         */
        virtual void
        convertPrimitiveVariablesToConservativeVariables(
            const std::vector<const Real*>& primitive_variables,
            const std::vector<Real*>& conservative_variables) = 0;
        
        /*
         * Convert primitive variables to conservative variables.
         */
        virtual void
        convertPrimitiveVariablesToConservativeVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables) = 0;
        
        /*
         * Check whether the given cell conservative variables are within the bounds.
         */
        virtual void
        checkCellDataOfConservativeVariablesBounded(
            HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables) = 0;
        
        /*
         * Check whether the given side conservative variables are within the bounds.
         */
        virtual void
        checkSideDataOfConservativeVariablesBounded(
            HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables) = 0;
        
        /*
         * Check whether the given cell primitive variables are within the bounds.
         */
        virtual void
        checkCellDataOfPrimitiveVariablesBounded(
            HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& primitive_variables) = 0;
        
        /*
         * Check whether the given side primitive variables are within the bounds.
         */
        virtual void
        checkSideDataOfPrimitiveVariablesBounded(
            HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables) = 0;
        
        /*
         * Register the required derived variables for transformation between conservative
         * variables and characteristic variables.
         */
        virtual void
        registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING::TYPE& averaging_type) = 0;
        
        /*
         * Register the required derived variables for transformation between primitive variables
         * and characteristic variables.
         */
        virtual void
        registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING::TYPE& averaging_type) = 0;
        
        /*
         * Get the number of projection variables for transformation between conservative
         * variables and characteristic variables.
         */
        virtual int
        getNumberOfProjectionVariablesForConservativeVariables() const = 0;
        
        /*
         * Get the number of projection variables for transformation between primitive variables
         * and characteristic variables.
         */
        virtual int
        getNumberOfProjectionVariablesForPrimitiveVariables() const = 0;
        
        /*
         * Compute the side data of the projection variables for transformation between conservative variables and
         * characteristic variables.
         */
        virtual void
        computeSideDataOfProjectionVariablesForConservativeVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables) = 0;
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        virtual void
        computeSideDataOfProjectionVariablesForPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables) = 0;
        
        /*
         * Compute the side data of characteristic variables from conservative variables.
         */
        virtual void
        computeSideDataOfCharacteristicVariablesFromConservativeVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables,
            const int& idx_offset) = 0;
        
        /*
         * Compute the side data of characteristic variables from primitive variables.
         */
        virtual void
        computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables,
            const int& idx_offset) = 0;
        
        /*
         * Compute the side data of conservative variables from characteristic variables.
         */
        virtual void
        computeSideDataOfConservativeVariablesFromCharacteristicVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables) = 0;
        
        /*
         * Compute the side data of primitive variables from characteristic variables.
         */
        virtual void
        computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables) = 0;
        
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
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * HAMERS_WEAK_PTR to FlowModel.
         */
        HAMERS_WEAK_PTR<FlowModel> d_flow_model;
        
        /*
         * Settings for projection variables.
         */
        AVERAGING::TYPE d_proj_var_conservative_averaging_type;
        AVERAGING::TYPE d_proj_var_primitive_averaging_type;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_HPP */
