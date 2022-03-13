#ifndef FILTER_HPP
#define FILTER_HPP

#include <string>

class EnsembleStatistics
{
    public:
        EnsembleStatistics(const std::string& object_name):
            d_object_name(object_name),
            d_num_ensembles(0)
        {
        }
        
        int getNumberOfEnsembles() const
        {
            return d_num_ensembles;
        }
        
        void incrementNumberOfEnsembles()
        {
            d_num_ensembles++;
        }
        
        virtual void clearAllData()
        {
        }
        
    protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Number of ensembles.
         */
        int d_num_ensembles;
};

#endif /* ENSEMBLE_STATISTICS_HPP */
