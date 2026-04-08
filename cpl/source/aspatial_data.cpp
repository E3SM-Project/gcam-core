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

/*!
 * \file aspatial_data.cpp
 * \brief This file processes gridded data either from file or from memory. It is an abstract
 * class from which the carbon scaler and emissions downscale code is derived.
 *
 * \author Kate Calvin and Alan Di Vittorio and Philip Myint
 */

#include <iostream>
#include <fstream>
#include <iomanip>

#include "util/base/include/auto_file.h"
#include "../include/aspatial_data.h"

using namespace std;

// Constructor
ASpatialData::ASpatialData(int aSize):
mLatVector(aSize, 0),
mLonVector(aSize, 0),
mIDVector(aSize, 0),
mValueVector(aSize, 0) {
}

// Constructor with option to read subregions in the mapping file. If false, subregions will be aggregated to their parent regions.
// Uses constructor delegation to call the main constructor to initialize the vectors
ASpatialData::ASpatialData(const int aSize, const bool aReadSubregions):
ASpatialData(aSize)
{
    mReadSubregions = aReadSubregions;
}

// Destructor
ASpatialData::~ASpatialData() {
}

// Read in spatial data from a space separated value file
// this is currently used to read in the base co2 spatial data
double ASpatialData::readSpatialData(std::string aFileName, bool aHasLatLon, bool aHasID, bool aCalcTotal) {
    // Create a double to store totals (if aCalcTotal == true)
    double total = 0.0;
    
    ifstream data(aFileName);
    if (!data.is_open())
    {
        ILogger& coupleLog = ILogger::getLogger( "coupling_log" );
        coupleLog.setLevel( ILogger::ERROR );
        coupleLog << "File not found: " << aFileName << endl;
        exit(EXIT_FAILURE);
    }
    string str;
    getline(data, str); // skip the first line
    int row = 0;
    while (getline(data, str))
    {
        istringstream iss(str);
        string token;
        double value;
        
        if ( aHasID ) {
            int id;
            
            // Parse ID
            getline(iss, token, ' ');
            id = std::stoi(token);
            mIDVector.at(row) = id;
         }
        
        if ( aHasLatLon ) {
            double lon;
            int lat;
            
            // Parse longitude
            getline(iss, token, ' ');
            lon = std::stod(token);
            mLonVector[row] = lon;
            
            // Parse latitude
            getline(iss, token, ' ');
            lat = std::stod(token);
            mLatVector[row] = lat;
        }
        
        // Parse Value
        getline(iss, token, ' ');
        value = std::stod(token);
        mValueVector[row] = value;
         
        // if aCalcTotal == true, then add this to the total
        total += value;
                
        row++;
    }
    
    return total;
}

// Read in spatial data from a space separated value file directly into an array
double ASpatialData::readSpatialData(std::string aFileName, bool aHasLatLon, bool aHasID, bool aCalcTotal, double *aValueArray) {
    // Create a double to store totals (if aCalcTotal == true)
    double total = 0.0;
    
    ifstream data(aFileName);
    if (!data.is_open())
    {
        ILogger& coupleLog = ILogger::getLogger( "coupling_log" );
        coupleLog.setLevel( ILogger::ERROR );
        coupleLog << "File not found: " << aFileName << endl;
        exit(EXIT_FAILURE);
    }
    string str;
    getline(data, str); // skip the first line
    int row = 0;
    while (getline(data, str))
    {
        istringstream iss(str);
        string token;
        double value;
        
        if ( aHasID ) {
            int id;
            
            // Parse ID
            getline(iss, token, ' ');
            id = std::stoi(token);
            mIDVector.at(row) = id;
        }
        
        if ( aHasLatLon ) {
            double lon;
            int lat;
            
            // Parse longitude
            getline(iss, token, ' ');
            lon = std::stod(token);
            mLonVector[row] = lon;
            
            // Parse latitude
            getline(iss, token, ' ');
            lat = std::stod(token);
            mLatVector[row] = lat;
        }
        
        // Parse Value
        getline(iss, token, ' ');
        value = std::stod(token);
        aValueArray[row] = value;
        
        // if aCalcTotal == true, then add this to the total
        total += value;
        
        row++;
    }

    return total;
}


// Read in spatial data from a csv file
// this is used for reading base elm data for scalar calculation
//   and also for the same in testing in the main program
double ASpatialData::readSpatialDataCSV(std::string aFileName, bool aHasLatLon, bool aHasID, bool aCalcTotal) {
    // Create a double to store totals (if aCalcTotal == true)
    double total = 0.0;
    
    ifstream data(aFileName);
    if (!data.is_open())
    {
        ILogger& coupleLog = ILogger::getLogger( "coupling_log" );
        coupleLog.setLevel( ILogger::ERROR );
        coupleLog << "File not found: " << aFileName << endl;
        exit(EXIT_FAILURE);
    }
    string str;
    getline(data, str); // skip the first line
    int row = 0;
    while (getline(data, str))
    {
        istringstream iss(str);
        string token;
        double value;
        
        if ( aHasID ) {
            int id;

            // Parse ID
            getline(iss, token, ',');
            id = std::stoi(token);
            mIDVector.at(row) = id;
         }
        
        if ( aHasLatLon ) {
            double lon;
            int lat;

            // Parse longitude
            getline(iss, token, ',');
            lon = std::stod(token);
            mLonVector[row] = lon;

            // Parse latitude
            getline(iss, token, ',');
            lat = std::stod(token);
            mLatVector[row] = lat;
        }

        // Parse value
        getline(iss, token, ',');
        value = std::stod(token);
        mValueVector[row] = value;

        // if aCalcTotal == true, then add this to the total
        total += value;

        row++;
    }
    
    return total;
}


// Read in spatial data from a csv file directly into an array
double ASpatialData::readSpatialDataCSV(std::string aFileName, bool aHasLatLon, bool aHasID, bool aCalcTotal, double *aValueArray) {
    // Create a double to store totals (if aCalcTotal == true)
    double total = 0.0;
    
    ifstream data(aFileName);
    if (!data.is_open())
    {
        ILogger& coupleLog = ILogger::getLogger( "coupling_log" );
        coupleLog.setLevel( ILogger::ERROR );
        coupleLog << "File not found: " << aFileName << endl;
        exit(EXIT_FAILURE);
    }
    string str;
    getline(data, str); // skip the first line
    int row = 0;
    while (getline(data, str))
    {
        istringstream iss(str);
        string token;
        double value;
        
        if ( aHasID ) {
            int id;

            // Parse ID
            getline(iss, token, ',');
            id = std::stoi(token);
            mIDVector.at(row) = id;
        }
        
        if ( aHasLatLon ) {
            double lon;
            int lat;

            // Parse longitude
            getline(iss, token, ',');
            lon = std::stod(token);
            mLonVector[row] = lon;

            // Parse latitude
            getline(iss, token, ',');
            lat = std::stod(token);
            mLatVector[row] = lat;
        }

        // Parse value
        getline(iss, token, ',');
        value = std::stod(token);
        aValueArray[row] = value;
        
        // if aCalcTotal == true, then add this to the total
        total += value;

        row++;
    }

    return total;
}


void ASpatialData::writeSpatialData(std::string aFileName, bool aWriteID) {
    ofstream oFile;
    oFile.open(aFileName);

    if( aWriteID ) {
       oFile << "yyyymm<ll>,lon_deg,lat_deg,";
    }
    oFile << "co2_kg/m2/s" << endl;

    for (int i = 0; i < mValueVector.size(); i++) {
        if( aWriteID ) {
            oFile << mIDVector[i] << "," <<  mLonVector[i] << "," << mLatVector[i] << ",";
        }
        oFile << fixed << setprecision(20); 
        oFile << mValueVector[i] << endl;
        oFile << defaultfloat;
    }
    oFile.close();
    
    return;
}

/*!
 * \brief Read regional mapping data from CSV file
 * 
 * Reads a CSV file that maps E3SM grid cells to GCAM regions. Each grid cell can map to 
 * one or more GCAM regions (for cells that straddle regional boundaries), with fractional 
 * weights indicating what portion of the cell belongs to each region.
 * 
 * \param aFileName Path to the regional mapping CSV file
 * 
 * \details
 * **CSV File Format:**
 * The input file must be a comma-separated file with a header row (skipped) and the following columns:
 * - Column 1: Region ID (skipped)
 * - Column 2: GLU ID (skipped)
 * - Column 3: Longitude index (integer)
 * - Column 4: Latitude index (integer)
 * - Column 5: Region name (string, may include quotes)
 * - Column 6: Subregion name (string, may include quotes)
 * - Column 7: Weight (double, 0.0 to 1.0)
 * 
 * **Example CSV:**
 * \code
 * region_id,glu_id,lon,lat,region,subregion,weight
 * 1,1,144,96,"USA","Northwest",1.0
 * 2,2,144,97,"USA","Northwest",0.6
 * 2,2,144,97,"Canada","BritishColumbia",0.4
 * \endcode
 * 
 * **Data Structures Populated:**
 * 
 * 1. `mRegionMapping`: Maps grid cell ID to list of regions
 *    - Key: "lon_lat" (e.g., "144_96")
 *    - Value: Vector of region IDs (e.g., ["USA.Northwest", "Canada.BritishColumbia"])
 *    - A grid cell may map to multiple regions if it straddles a boundary
 * 
 * 2. `mRegionWeights`: Maps (grid cell, region) pair to fractional weight
 *    - Key: pair<string, string> (gridID, regionID)
 *    - Value: double (0.0 to 1.0) representing fraction of cell in that region
 *    - Weights for a single grid cell across all regions should sum to 1.0
 * 
 */
// TODO: Modify all subclasses of ASpatialData (e.g., CarbonScalers, EmissDownscale, DegreeDays) to use this 
// method to read in the regional mapping data
void ASpatialData::readRegionalMappingData(const std::string &aFileName) 
{
    // Open the file and check that it opened successfully
    std::ifstream data(aFileName);
    if (!data.is_open())
    {
        exit(EXIT_FAILURE);
    }
    std::string line;
    // Skip the first line (i.e., the header)
    std::getline(data, line); 

    // Read each line of the file one by one
    while (std::getline(data, line))
    {
        std::istringstream iss(line);
        std::string token;
        int lon;
        int lat;
        std::string region;
        std::string subregion;
        double weight;
            
        // Skip region & GLU ID
        std::getline(iss, token, ',');
        std::getline(iss, token, ',');
        
        // Parse longitude
        std::getline(iss, token, ',');
        lon = std::stoi(token);
        
        // Parse latitude
        std::getline(iss, token, ',');
        lat = std::stoi(token);
        
        // Example gridID: "144_96", where 144 is the longitude index and 96 is the latitude index
        std::string gridID = std::to_string(lon) + "_" + std::to_string(lat);
        
        // Parse Region ID
        std::getline(iss, token, ',');
        region = token;
        // Remove quotes from region name, if they exist. TODO: Also remove spaces?
        region.erase( remove( region.begin(), region.end(), '\"' ), region.end() );
        
        std::string regID;
        if (mReadSubregions)
        {
            // Parse Subregion ID
            std::getline(iss, token, ',');
            subregion = token;
            // Remove quotes from subregion name, if they exist. TODO: Also remove spaces?
            subregion.erase( remove( subregion.begin(), subregion.end(), '\"' ), subregion.end() );

            // Example region ID: "USA.Northwest", where "USA" is the region and "Northwest" is the subregion (e.g., basin, state, or other sub-regional unit)
            regID = region + "." + subregion;

            // Add region ID to the mapping vector. Note that there maybe more than one regID per gridID (hence, a vector). If gridID is not found, this will 
            // automatically create a new vector with this regID as the vector's first element. If gridID is found, this will add to the existing regID vector
            mRegionMapping[gridID].push_back(regID);

            // Parse weight -- this is the fraction of the grid cell in a particular GCAM region
            std::getline(iss, token, ',');
            weight = std::stod(token);

            // If reading subregions, then the key is the specific region-subregion combination (e.g., "USA.Northwest"), 
            // and it will be unique, so we can just add the entry to the map. This will also create the map entry if it does not already exist
            mRegionWeights[std::make_pair(gridID, regID)] = weight;
        }
        else
        {
            // If not reading subregions, then skip reading subregions so that the region ID is just the main region (e.g., "USA")
            std::getline(iss, token, ',');
            regID = region;

            // Parse weight -- this is the fraction of the grid cell in a particular GCAM region
            std::getline(iss, token, ',');
            weight = std::stod(token);

            // Populate the region mapping and weights maps
            auto it = mRegionMapping.find(gridID);
            if (it == mRegionMapping.end())
            {
                // If gridID is not found, then create a vector with this regID and add it to the map
                mRegionMapping[gridID] = std::vector<std::string>{regID};

                // Set the region weight for this grid cell and region combination
                mRegionWeights[std::make_pair(gridID, regID)] = weight;
            }
            else
            {
                // If gridID is found, then add to the existing regID vector only if it is not already in the vector (to avoid duplicates when 
                // subregions are being aggregated to the same parent region)
                auto &regIDs = it->second;
                if (std::find(regIDs.begin(), regIDs.end(), regID) == regIDs.end())
                {
                    regIDs.push_back(regID);
                }

                // If regID is not already in the map, the [] operator will create a new map entry with a default value of 0, so we can just add the weight
                // in both the case where the regID is new or already exists in the map
                mRegionWeights[std::make_pair(gridID, regID)] += weight;
            }

        } // if (mReadSubregions)

    } // while (std::getline(data, line))

    data.close();
}

void ASpatialData::readMapping(std::string aFileName) {
    return;
}


void ASpatialData::setValueVector(std::vector<double> aValueVector) {
    mValueVector = aValueVector;
    return;
}

void ASpatialData::setIDVector(std::vector<int> aIDVector) {
    mIDVector = aIDVector;
    return;
}

void ASpatialData::setLonVector(std::vector<double> aLonVector) {
    mLonVector = aLonVector;
    return;
}

void ASpatialData::setLatVector(std::vector<double> aLatVector) {
    mLatVector = aLatVector;
    return;
}

std::vector<double> ASpatialData::getValueVector() {
    return mValueVector;
}

std::vector<int> ASpatialData::getIDVector() {
    return mIDVector;
}

std::vector<double> ASpatialData::getLatVector() {
    return mLatVector;
}

std::vector<double> ASpatialData::getLonVector() {
    return mLonVector;
}

