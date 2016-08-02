#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

#include "flow/flow_models/FlowModels.hpp"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructor(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_conv_ghosts(hier::IntVector::getZero(d_dim)),
                d_num_eqn(num_eqn),
                d_num_species(num_species),
                d_flow_model(flow_model),
                d_convective_flux_reconstructor_db(convective_flux_reconstructor_db)
        {}
        
        virtual ~ConvectiveFluxReconstructor() {}
        
        /*
         * Get the number of ghost cells needed by the convective flux
         * reconstructor.
         */
        hier::IntVector
        getConvectiveFluxNumberOfGhostCells(void) const
        {
            return d_num_conv_ghosts;
        }
        
        /*
         * Print all characteristics of the convective flux reconstruction class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the convective flux reconstruction class
         * into the restart database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the convective fluxes and sources due to hyperbolization
         * of the equations.
         */
        virtual void
        computeConvectiveFluxesAndSources(
            hier::Patch& patch,
            const double time,
            const double dt,
            const int RK_step_number,
            const boost::shared_ptr<pdat::FaceVariable<double> >& variable_convective_flux,
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
            const boost::shared_ptr<hier::VariableContext>& data_context) = 0;
    
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
         * Number of ghost cells needed by the shock capturing scheme.
         */
        hier::IntVector d_num_conv_ghosts;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Flow model.
         */
        const boost::shared_ptr<FlowModel> d_flow_model;
        
        /*
         * boost::shared_ptr to database of the convective_flux_reconstructor.
         */
        const boost::shared_ptr<tbox::Database> d_convective_flux_reconstructor_db;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_HPP */
