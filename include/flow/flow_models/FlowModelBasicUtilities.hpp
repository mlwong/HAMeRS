#ifndef FLOW_MODEL_BASIC_UTILITIES_HPP
#define FLOW_MODEL_BASIC_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "flow/flow_models/FlowModel.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include "boost/weak_ptr.hpp"
#include <string>

namespace AVERAGING_TMP
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
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species),
                d_num_eqn(num_eqn),
                d_proj_var_conservative_averaging_type(AVERAGING_TMP::SIMPLE),
                d_proj_var_primitive_averaging_type(AVERAGING_TMP::SIMPLE)
        {}
        
        virtual ~FlowModelBasicUtilities() {}
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const boost::weak_ptr<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Convert conservative variables to primitive variables.
         */
        virtual void
        convertConservativeVariablesToPrimitiveVariables(
            const std::vector<const double*>& conservative_variables,
            const std::vector<double*>& primitive_variables) = 0;
        
        /*
         * Convert conservative variables to primitive variables.
         */
        virtual void
        convertConservativeVariablesToPrimitiveVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables) = 0;
        
        /*
         * Convert primitive variables to conservative variables.
         */
        virtual void
        convertPrimitiveVariablesToConservativeVariables(
            const std::vector<const double*>& primitive_variables,
            const std::vector<double*>& conservative_variables) = 0;
        
        /*
         * Convert primitive variables to conservative variables.
         */
        virtual void
        convertPrimitiveVariablesToConservativeVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables) = 0;
        
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
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables) = 0;
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        virtual void
        computeSideDataOfProjectionVariablesForPrimitiveVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables) = 0;
        
        /*
         * Compute the side data of characteristic variables from conservative variables.
         */
        virtual void
        computeSideDataOfCharacteristicVariablesFromConservativeVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
            const int& idx_offset) = 0;
        
        /*
         * Compute the side data of characteristic variables from primitive variables.
         */
        virtual void
        computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
            const int& idx_offset) = 0;
        
        /*
         * Compute the side data of conservative variables from characteristic variables.
         */
        virtual void
        computeSideDataOfConservativeVariablesFromCharacteristicVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables) = 0;
        
        /*
         * Compute the side data of primitive variables from characteristic variables.
         */
        virtual void
        computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables) = 0;
        
        /*
         * Check whether the given cell conservative variables are within the bounds.
         */
        virtual void
        checkCellDataOfConservativeVariablesBounded(
            boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables) = 0;
        
        /*
         * Check whether the given side conservative variables are within the bounds.
         */
        virtual void
        checkSideDataOfConservativeVariablesBounded(
            boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables) = 0;
        
        /*
         * Check whether the given cell primitive variables are within the bounds.
         */
        virtual void
        checkCellDataOfPrimitiveVariablesBounded(
            boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables) = 0;
        
        /*
         * Check whether the given side primitive variables are within the bounds.
         */
        virtual void
        checkSideDataOfPrimitiveVariablesBounded(
            boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables) = 0;
        
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
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * boost::weak_ptr to FlowModel.
         */
        boost::weak_ptr<FlowModel> d_flow_model;
        
        /*
         * Settings for projection variables.
         */
        AVERAGING_TMP::TYPE d_proj_var_conservative_averaging_type;
        AVERAGING_TMP::TYPE d_proj_var_primitive_averaging_type;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_HPP */
