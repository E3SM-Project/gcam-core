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

        void readRegionalMappingData(const std::string &aFileName);

        void aggregateDegreeDays(const int aGCAMYear, const double *const aELMArea, const double *const aELMDegreeDays,  
            const double *const aPopDensity, const double *const aELMLandFrac, std::vector<int> &aYears, std::vector<std::string> &aRegions, 
            std::vector<double> &aDegreeDaysVector, int &aNumValues);

        void createDegreeDayVectors(const int aGCAMYear, std::vector<int> &aYears, std::vector<std::string> &aRegions,
                            std::vector<double> &aDegreeDaysVector, const std::map<std::string, double> &aDegreeDaysMap);

        void writeDegreeDays(const std::string &aFileName, const std::vector<int> &aYears, const std::vector<std::string> &aRegions,
                const std::vector<double> &aDegreeDaysVector, const int aLength);

        int readDegreeDays(const std::string &aFileName, std::vector<int> &aYears, std::vector<std::string> &aRegions,
                std::vector<double> &aDegreeDaysVector);
    private:    
        // Boolean indicating whether to read subregions in the mapping file. If false, subregions will be aggregated to their parent regions
        bool mReadSubregions;

        // Number of latitude and longitude values. Storing these so they do not have to be passed to every method
        int mNumLat;
        int mNumLon;
        
        // Map grid cells to regions. Key is a string with longitude and latitude ("lon_lat").
        // Key maps to a vector of strings containing the region and subregion. For example, "RegionA.Subregion1". 
        // Note that this map will contain more than lat * lon entries since some grid cells map to multiple regions
        std::map<std::string, std::vector<std::string>> mRegionMapping;
        
        // Map region weights (these indicate the fraction of a grid cell assigned to each region)
        // Key is a pair indicating the grid cell and the region/subregion.
        // Key maps to a double representing the fraction of the grid cell in that region/subregion
        std::map<std::pair<std::string,std::string>, double> mRegionWeights;
};

#endif // __DEGREE_DAYS__
