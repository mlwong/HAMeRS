#ifndef FLOW_MODEL_HELPER_HPP
#define FLOW_MODEL_HELPER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "util/MPI_helpers/MPIHelper.hpp"

#include <string>

using namespace SAMRAI;

class FlowModelMPIHelper: public MPIHelper
{
    public:
        FlowModelMPIHelper(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const bool use_diffusive_flux_utilities = false):
                MPIHelper(
                    object_name,
                    dim,
                    grid_geometry,
                    patch_hierarchy),
                d_flow_model(flow_model),
                d_use_diffusive_flux_utilities(use_diffusive_flux_utilities)
        {}
        
    protected:
        /*
         * Register the patch and data in the flow model and compute the corresponding
         * average.
         */
        void setupFlowModelAndRegisterPatchWithDataContext(
            const hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        void registerDerivedVariables(
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data);
        
        void allocateMemoryForDerivedCellData();
        
        void computeDerivedCellData();
        
        /*
         * Get the cell data of one cell variable in the registered patch.
         */
        HAMERS_SHARED_PTR<pdat::CellData<double> > getCellData(const std::string& variable_key);
        
        void unregisterPatch();
        
        /*
         * HAMERS_SHARED_PTR to FlowModel.
         */
        HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * Whether to use diffusive flux utilities.
         */
        const bool d_use_diffusive_flux_utilities;
        
        HAMERS_SHARED_PTR<FlowModelDiffusiveFluxUtilities> d_diffusive_flux_utilities;
        
};

#endif /* FLOW_MODEL_HELPER_HPP */
