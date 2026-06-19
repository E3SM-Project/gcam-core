#ifndef __DEGREE_DAYS__
#define __DEGREE_DAYS__

/*
 * LEGAL NOTICE
 * This computer software was prepared by Battelle Memorial Institute,
 * hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
 * with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
 * CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this
 * sentence must appear on any copies of this computer software.
 *
 * EXPORT CONTROL
 * User agrees that the Software will not be shipped, transferred or
 * exported into any country or used in any manner prohibited by the
 * United States Export Administration Act or any other applicable
 * export laws, restrictions or regulations (collectively the "Export Laws").
 * Export of the Software may require some form of license or other
 * authority from the U.S. Government, and failure to obtain such
 * export control license may result in criminal liability under
 * U.S. laws. In addition, if the Software is identified as export controlled
 * items under the Export Laws, User represents and warrants that User
 * is not a citizen, or otherwise located within, an embargoed nation
 * (including without limitation Iran, Syria, Sudan, Cuba, and North Korea)
 *     and that User is not otherwise prohibited
 * under the Export Laws from receiving the Software.
 *
 * Copyright 2011 Battelle Memorial Institute.  All Rights Reserved.
 * Distributed as open-source under the terms of the Educational Community
 * License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php
 *
 * For further details, see: http://www.globalchange.umd.edu/models/gcam/
 *
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "util/base/include/auto_file.h"
#include "../include/aspatial_data.h"

class DegreeDays : public ASpatialData 
{
    public:
        DegreeDays(const int aNumLon, const int aNumLat, const bool aReadSubregions);

        ~DegreeDays();

        void aggregateDegreeDays(const int aGCAMYear, const double *const aELMArea, const double *const aELMHDD, const double *const aELMCDD,
            const double *const aPopDensity, const double *const aELMLandFrac, std::vector<int> &aYears, std::vector<std::string> &aRegions,
            std::vector<double> &aELMHDDVector, std::vector<double> &aELMCDDVector, int &aNumValues, std::string aBaseHDDFileName, std::string aBaseCDDFileName);

        void createDegreeDaysVectors(const int aGCAMYear, std::vector<int> &aYears, std::vector<std::string> &aRegions,
                            std::vector<double> &aDegreeDaysVector, const std::map<std::string, double> &aDegreeDaysMap);

        void writeDegreeDays(const std::string &aFileName, const std::vector<int> &aYears, const std::vector<std::string> &aRegions,
                const std::vector<double> &aHDDVector, const std::vector<double> &aCDDVector, const int aLength);

        int readDegreeDays(const std::string &aFileName, std::vector<int> &aYears, std::vector<std::string> &aRegions,
                std::vector<double> &aHDDVector, std::vector<double> &aCDDVector);

        void readBaseDegreeDays(const std::string aBaseHDDFileName, const std::string aBaseCDDFileName);


    private:    
        // Number of latitude and longitude values. Storing these so they do not have to be passed to every method
        int mNumLat;
        int mNumLon;

        // Baseline degree days for scaling (gridded data at lat/lon level)
        std::vector<double> mBaseHDDVector;
        std::vector<double> mBaseCDDVector;

        // Scalars recalculated every coupled year from the current ELM aggregation and the baseline.
        // Maps region ID to scalar value (scalar = ELM_current_year / baseline), bounded to [0.125, 2.0].
        // These are applied to GCAM's degree-days values.
        std::map<std::string, double> mHDDScalars;
        std::map<std::string, double> mCDDScalars;
};

#endif // __DEGREE_DAYS__
