#ifndef FLOW_MODEL_SPECIAL_SOURCE_TERMS_HPP
#define FLOW_MODEL_SPECIAL_SOURCE_TERMS_HPP

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

class FlowModelSpecialSourceTerms
{
    public:
        FlowModelSpecialSourceTerms(
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
            if (source_terms_db->keyExists("special_source_box_lo"))
            {
                d_special_source_box_lo = source_terms_db->getRealVector("special_source_box_lo");
            }
            else if (source_terms_db->keyExists("d_special_source_box_lo"))
            {
                d_special_source_box_lo = source_terms_db->getRealVector("d_special_source_box_lo");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'special_source_box_lo' or 'd_special_source_box_lo' found in data for source terms."
                    << std::endl);
            }
            
            if (source_terms_db->keyExists("special_source_box_hi"))
            {
                d_special_source_box_hi = source_terms_db->getRealVector("special_source_box_hi");
            }
            else if (source_terms_db->keyExists("d_special_source_box_hi"))
            {
                d_special_source_box_hi = source_terms_db->getRealVector("d_special_source_box_hi");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'special_source_box_hi' or 'd_special_source_box_hi' found in data for source terms."
                    << std::endl);
            }

            // AFK 8/16/23 bottom sponge region 
            if (source_terms_db->keyExists("special_source_box_bo"))
            {
                d_special_source_box_bo = source_terms_db->getRealVector("special_source_box_bo");
            }
            else if (source_terms_db->keyExists("d_special_source_box_bo"))
            {
                d_special_source_box_bo = source_terms_db->getRealVector("d_special_source_box_bo");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'special_source_box_bo' or 'd_special_source_box_bo' found in data for source terms."
                    << std::endl);
            }
            
            // AFK 8/16/23 top sponge region 
            if (source_terms_db->keyExists("special_source_box_to"))
            {
                d_special_source_box_to = source_terms_db->getRealVector("special_source_box_to");
            }
            else if (source_terms_db->keyExists("d_special_source_box_to"))
            {
                d_special_source_box_to = source_terms_db->getRealVector("d_special_source_box_to");
            }
            else
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "No key 'special_source_box_to' or 'd_special_source_box_to' found in data for source terms."
                    << std::endl);
            }

            d_special_source_exterior = source_terms_db->getBoolWithDefault("special_source_exterior", true);
            d_special_source_exterior = source_terms_db->getBoolWithDefault("d_special_source_exterior", d_special_source_exterior);
            
            if (d_special_source_box_lo.size() != d_dim.getValue())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSpecialSourceTerms::FlowModelSpecialSourceTerms\n"
                    << "Size of 'special_source_box_lo' or 'd_special_source_box_lo' is not consistent with problem dimension."
                    << std::endl);
            }
            
            if (d_special_source_box_hi.size() != d_dim.getValue())
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSpecialSourceTerms::FlowModelSpecialSourceTerms\n"
                    << "Size of 'special_source_box_hi' or 'd_special_source_box_hi' is not consistent with problem dimension."
                    << std::endl);
            }

            if (d_special_source_box_bo.size() != d_dim.getValue()) // AFK 8/16/23
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSpecialSourceTerms::FlowModelSpecialSourceTerms\n"
                    << "Size of 'special_source_box_bo' or 'd_special_source_box_bo' is not consistent with problem dimension."
                    << std::endl);
            }

            if (d_special_source_box_to.size() != d_dim.getValue()) // AFK 8/16/23
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelSpecialSourceTerms::FlowModelSpecialSourceTerms\n"
                    << "Size of 'special_source_box_to' or 'd_special_source_box_to' is not consistent with problem dimension."
                    << std::endl);
            }
        }
        
        /*
         * Add the effects of the special source terms.
         */
        void
        computeSpecialSourceTermsOnPatch(
            HAMERS_SHARED_PTR<pdat::CellData<Real> >& source,
            const hier::Patch& patch,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables,
            const std::unordered_map<std::string, Real>& monitoring_statistics_map,
            const double time,
            const double dt,
            const int RK_step_number);
        
        std::vector<Real> getSpecialSourceBoxLo() const
        {
            return d_special_source_box_lo;
        }
        
        std::vector<Real> getSpecialSourceBoxHi() const
        {
            return d_special_source_box_hi;
        }

        std::vector<Real> getSpecialSourceBoxBo() const // AFK 8/16/23
        {
            return d_special_source_box_bo;
        }

        std::vector<Real> getSpecialSourceBoxTo() const // AFK 8/16/23
        {
            return d_special_source_box_to;
        }
        
        bool getSpecialSourceExterior() const
        {
            return d_special_source_exterior;
        }
        
        virtual void putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db);
        
    private:
        void
        putToRestartBase(const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db)
        {
            restart_source_terms_db->putRealVector("d_special_source_box_lo", d_special_source_box_lo);
            restart_source_terms_db->putRealVector("d_special_source_box_hi", d_special_source_box_hi);
            restart_source_terms_db->putRealVector("d_special_source_box_bo", d_special_source_box_bo); // AFK 8/16/23
            restart_source_terms_db->putRealVector("d_special_source_box_to", d_special_source_box_to); // AFK 8/16/23
        
            restart_source_terms_db->putBool("d_special_source_exterior", d_special_source_exterior);
        }
        
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
         * Special source box definition.
         */
        std::vector<Real> d_special_source_box_lo; // Lower spatial coordinates.
        std::vector<Real> d_special_source_box_hi; // Uppe/r spatial coordinates.
        std::vector<Real> d_special_source_box_bo; // AFK 8/16/23 Bottom spatial coordinates.
        std::vector<Real> d_special_source_box_to; // AFK 8/16/23 Top spatial coordinates.
        bool d_special_source_exterior;
        
};

#endif /* FLOW_MODEL_SPECIAL_SOURCE_TERMS_HPP */