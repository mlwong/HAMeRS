#ifndef FLOW_MODEL_SPONGE_HPP
#define FLOW_MODEL_SPONGE_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class FlowModelSponge
{
    public:
        FlowModelSponge(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<tbox::Database>& source_terms_db,
            const int& num_eqn):
                d_object_name(object_name),
                d_project_name(project_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_source_terms_db(source_terms_db),
                d_num_eqn(num_eqn)
        {
            if (source_terms_db->keyExists("sponge_box_lo"))
            {
                d_sponge_box_lo = source_terms_db->getDoubleVector("sponge_box_lo");
            }
            else if (source_terms_db->keyExists("d_sponge_box_lo"))
            {
                d_sponge_box_lo = source_terms_db->getDoubleVector("d_sponge_box_lo");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'sponge_box_lo' or 'd_sponge_box_lo' found in data for source terms."
                    << std::endl);
            }
            
            if (source_terms_db->keyExists("sponge_box_hi"))
            {
                d_sponge_box_hi = source_terms_db->getDoubleVector("sponge_box_hi");
            }
            else if (source_terms_db->keyExists("d_sponge_box_hi"))
            {
                d_sponge_box_hi = source_terms_db->getDoubleVector("d_sponge_box_hi");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'sponge_box_hi' or 'd_sponge_box_hi' found in data for source terms."
                    << std::endl);
            }
            
            d_sponge_exterior = source_terms_db->getBoolWithDefault("sponge_exterior", true);
            d_sponge_exterior = source_terms_db->getBoolWithDefault("d_sponge_exterior", d_sponge_exterior);
            
            if (source_terms_db->keyExists("sponge_rate"))
            {
                d_sponge_rate = source_terms_db->getDouble("sponge_rate");
            }
            else if (source_terms_db->keyExists("d_sponge_rate"))
            {
                d_sponge_rate = source_terms_db->getDouble("d_sponge_rate");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'sponge_rate' or 'd_sponge_rate' found in data for source terms."
                    << std::endl);
            }
            
            if (d_sponge_box_lo.size() != d_dim.getValue())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSponge::FlowModelSponge\n"
                    << "Size of 'sponge_box_lo' or 'd_sponge_box_lo' is not consistent with problem dimension."
                    << std::endl);
            }
            
            if (d_sponge_box_hi.size() != d_dim.getValue())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSponge::FlowModelSponge\n"
                    << "Size of 'sponge_box_hi' or 'd_sponge_box_hi' is not consistent with problem dimension."
                    << std::endl);
            }
        }
        
        void
        computeSpongeSourceTermsOnPatch(
            HAMERS_SHARED_PTR<pdat::CellData<double> >& source,
            const hier::Patch& patch,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
            const double time,
            const double dt,
            const int RK_step_number);
        
        std::vector<double> getSpongeBoxLo() const
        {
            return d_sponge_box_lo;
        }
        
        std::vector<double> getSpongeBoxHi() const
        {
            return d_sponge_box_hi;
        }
        
        bool getSpongeExterior() const
        {
            return d_sponge_exterior;
        }
        
        double getSpongeRate() const
        {
            return d_sponge_rate;
        }
        
        void
        putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db)
        {
            restart_source_terms_db->putDoubleVector("d_sponge_box_lo", d_sponge_box_lo);
            restart_source_terms_db->putDoubleVector("d_sponge_box_hi", d_sponge_box_hi);
        
            restart_source_terms_db->putBool("d_sponge_exterior", d_sponge_exterior);
            restart_source_terms_db->putDouble("d_sponge_rate",   d_sponge_rate);
        }
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Name of the project.
         */
        std::string d_project_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Database of source terms.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_source_terms_db;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Sponge box definition.
         */
        std::vector<double> d_sponge_box_lo; // Lower spatial coordinates.
        std::vector<double> d_sponge_box_hi; // Upper spatial coordinates.
        bool d_sponge_exterior;
        double d_sponge_rate;
        
};

#endif /* FLOW_MODEL_SPONGE_HPP */