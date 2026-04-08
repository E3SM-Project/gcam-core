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
 * \file runoff_data.cpp
 * \brief This file aggregates gridded data on runoff from E3SM into basin-specific sums for GCAM
 *
 * \author Philip Myint
 */

#include "../include/runoff_data.h"

// Constructor
RunoffData::RunoffData(const int aNumLon, const int aNumLat, const bool aReadSubregions):
ASpatialData(aNumLon * aNumLat, aReadSubregions),
mNumLon(aNumLon),
mNumLat(aNumLat)
{
}

// Destructor (could just use the default destructor since there are no dynamically allocated raw pointers, 
// but including an explicit destructor in case we want to add any cleanup code in the future)
RunoffData::~RunoffData() 
{
}

/*!
 * \brief Aggregate gridded runoff data to basin scale using area and land fraction weighting
 * 
 * This method performs aggregation of runoff data from E3SM's gridded output to GCAM's basin scale. 
 * 
 * \param aGCAMYear [INPUT] GCAM year for this data (e.g., 2025, 2050)
 * \param aELMArea [INPUT] Grid cell areas in km² (mNumLat × mNumLon elements)
 * \param aELMRunoffData [INPUT] Runoff data for each grid cell (mNumLat × mNumLon elements)
 * \param aELMLandFrac [INPUT] Land fraction for each grid cell (mNumLat × mNumLon elements, values between 0 and 1)
 * \param aYears [OUTPUT] Vector of years (one entry per basin, all equal to aGCAMYear)
 * \param aRegions [OUTPUT] Vector of region names (e.g., "USA", "China", "India")
 * \param aBasins [OUTPUT] Vector of basin names (e.g., "Mississippi", "Yangtze", "Ganges")
 * \param aRunoffDataVector [OUTPUT] Vector of runoff data by basin
 * \param aNumValues [OUTPUT] Number of basin values written to output vectors
 * 
 * \details
 * 
 * \pre mRegionMapping and mRegionWeights must be populated by ASpatialData::readRegionalMappingData()
 * \pre Input arrays must have mNumLat × mNumLon elements
 * \pre aELMArea values must be positive (km²)
 * \pre aELMLandFrac values must be between 0 and 1
 * 
 * \post Output vectors contain one entry per basin
 * \post aNumValues equals the number of basins processed
 * 
 * \note Ocean cells (not in mRegionMapping) are automatically skipped
 * 
 */
void RunoffData::aggregateRunoffData(const int aGCAMYear, const double *const aELMArea, const double *const aELMRunoffData, 
    const double *const aELMLandFrac, std::vector<int> &aYears, std::vector<std::string> &aRegions, std::vector<std::string> &aBasins, 
    std::vector<double> &aRunoffDataVector, int &aNumValues) 
{
    // TODO: Generate any diagnostic files or exclude outliers? Probably not necessary or relevant for runoff data, but can add if desired

    // Create map to store information on runoff by basin, which will be used to create the output vectors for GCAM
    std::map<std::string, double> aRunoffDataMap;
    
    // Note: E3SM data will have longitude moving fastest, then latitude, so loop so have longitude is the inner loop and latitude is the outer loop
    int gridIndex; 

    // Conversion factor from mm/s (E3SM units) to km^3/year (GCAM units); need to multiply this number by area in km^2, as done below, to get km^3/year
    constexpr double secsPerYear = 365.25 * 24.0 * 3600.0; // seconds in a year
    constexpr double mmTOkm = 1.0e-6; // conversion from mm to km
    constexpr double runoffConversionFactor = secsPerYear * mmTOkm; // overall conversion factor

    for (int k = 1; k <= mNumLat; k++) 
    {
        for (int j = 1; j <= mNumLon; j++) 
        {
            gridIndex = (k - 1) * mNumLon + (j - 1);

            // Get the gridID for this grid cell based on the longitude and latitude indices
            std::string gridID = std::to_string(j) + "_" + std::to_string(k);
            
            auto tempGrid = mRegionMapping.find(gridID);
            if (tempGrid == mRegionMapping.end()) 
            {
                // Grid isn't found in the mapping. Currently, this probably means it is an ocean grid.
                // TODO: set up loop only over land grids, either using the mRegionMapping or one of the files from ELM (same for carbon scalers)
            } 
            else 
            {
                // Get the vector of regionIDs for this grid cell
                const std::vector<std::string> &regInGrd = tempGrid->second;

                // Loop over all regions this grid is mapped to and calculate the contribution to the total runoff of each region.
                for (const auto &regID : regInGrd) 
                {
                    // Calculate the contribution of this grid cell to the total runoff in each region that the grid cell is a part of, which involves
                    // the fraction of the grid cell in this region (from the mapping file), the area of the grid cell, and the fraction of the grid cell 
                    // that is actually land. TODO: verify whether land fraction needs to be included. Perform sum (as assumed here) or weighted average?
                    aRunoffDataMap[regID] += mRegionWeights[std::make_pair(gridID, regID)] * aELMArea[gridIndex] * 
                        aELMLandFrac[gridIndex] * aELMRunoffData[gridIndex] * runoffConversionFactor; 

                } // for (const auto &regID : regInGrd) 

            } // if (tempGrid == mRegionMapping.end()) 

        } // for (int j = 1; j <= mNumLon; j++)

    } // (int k = 1; k <= mNumLat; k++)
    
    // Create output vectors for GCAM. Loop through the maps and create the vectors as needed
    createRunoffDataVectors(aGCAMYear, aYears, aRegions, aBasins, aRunoffDataVector, aRunoffDataMap);

    // Set the number of actual records
    aNumValues = static_cast<int>(aRunoffDataMap.size());
}


/*!
 * \brief Convert map structures to parallel vectors for GCAM interface
 * 
 * This helper method packages runoff data from map structures (used for efficient
 * aggregation) into parallel vectors (required by GCAM's interface). The output vectors
 * are organized such that the same index across all vectors represents data for one region + basin combination.
 * 
 * \param aGCAMYear [INPUT] GCAM year to assign to all entries
 * \param aYears [OUTPUT] Vector of years (all entries equal to aGCAMYear)
 * \param aRegions [OUTPUT] Vector of region names extracted from map keys
 * \param aBasins [OUTPUT] Vector of basin names extracted from map keys
 * \param aRunoffDataVector [OUTPUT] Vector of runoff data in same order as aRegions
 * \param aRunoffDataMap [INPUT] Map of region name → runoff data value
 * 
 * \details
 * 
 * **Output Vector Structure (Parallel Vectors):**
 * \code
 * Index:       0         1          2
 * aYears:     [2050,    2050,      2050]
 * aRegions:   ["USA",   "China",   "India"]
 * aBasins:    ["Mississippi", "Yangtze", "Ganges"]
 * aRunoffDataVector: [2845.3,  3421.9,    456.8]
 * \endcode
 * 
 * **Parallel Vector Concept:**
 * 
 * All vectors have the same length, and index `i` across all vectors represents the
 * complete data for region+basin `i`:
 * - `aYears[i]` = year for region i
 * - `aRegions[i]` = name of region in region+basin i  
 * - `aBasins[i]` = name of basin in in region+basin i
 * - `aRunoffDataVector[i]` = runoff data for region i
 * 
 */
void RunoffData::createRunoffDataVectors(const int aGCAMYear, std::vector<int> &aYears, std::vector<std::string> &aRegions, 
    std::vector<std::string> &aBasins, std::vector<double> &aRunoffDataVector, const std::map<std::string, double> &aRunoffDataMap) 
{
    // Clear the vectors first to ensure they are empty before adding data
    aYears.clear();
    aRegions.clear();
    aBasins.clear();
    aRunoffDataVector.clear();

    // Token and string stream for parsing region names (i.e., to split region and basin combinations into separate vectors)
    std::string token;
    std::istringstream tokenStream;

    // Loop through the map and create the vectors
    for (const auto &pair : aRunoffDataMap) 
    {
        const std::string &regionIDPlusBasin = pair.first;

        // Parse the region name to extract the region and basin
        tokenStream.clear();
        tokenStream.str(regionIDPlusBasin);
        // Get the region name (before the dot)
        std::getline(tokenStream, token, '.'); 
        std::string regionID = token;
        // Get the basin name (after the dot)
        std::getline(tokenStream, token, '.'); 
        std::string basinName = token;

        aRunoffDataVector.push_back(pair.second);
        
        aYears.push_back(aGCAMYear);

        aRegions.push_back(regionID);

        aBasins.push_back(basinName);
    }
}

/*!
 * \brief Write runoff data to CSV file
 * 
 * Writes runoff data to a comma-separated values (CSV) file with maximum
 * precision to ensure exact round-trip conversion between binary and text representations.
 * This is important for reproducibility and restart capabilities.
 * 
 * \param aFileName [INPUT] Path to output CSV file (will be created/overwritten)
 * \param aYears [INPUT] Vector of years (parallel to other vectors)
 * \param aRegions [INPUT] Vector of region names (parallel to other vectors)
 * \param aBasins [INPUT] Vector of basin names (parallel to other vectors)
 * \param aRunoffDataVector [INPUT] Vector of runoff data by region+basin
 * \param aLength [INPUT] Number of entries to write (should be ≤ vector sizes)
 * 
 * \details
 * **CSV Output Format:**
 * \code
 * Year,Region,Basin,RunoffData
 * 2050,USA,Mississippi,2845.342145678901234
 * 2050,China,Yangtze,3421.876543210987
 * 2050,India,Ganges,456.7890123456789
 * \endcode
 * 
 * **Precision Handling:**
 * 
 * The file is written with `std::numeric_limits<double>::max_digits10` precision
 * (typically 17 decimal digits for double). This ensures that:
 * - Values can be read back in exactly as they were written
 * - No precision is lost in the round-trip (binary → text → binary)
 * - Results are reproducible across different systems
 * 
 * **Column Descriptions:**
 * - **Year**: GCAM year (integer, e.g., 2020, 2025, 2050)
 * - **Region**: GCAM region name (string, e.g., "USA", "China", "India")
 * - **Basin**: Basin name (string, e.g., "Mississippi", "Yangtze", "Ganges")
 * - **RunoffData**: Runoff data (double, mm/year or similar units)
 *
 */
void RunoffData::writeRunoffData(const std::string &aFileName, const std::vector<int> &aYears, const std::vector<std::string> &aRegions, 
    const std::vector<std::string> &aBasins, const std::vector<double> &aRunoffDataVector, const int aLength) 
{
    std::ofstream oFile;
    oFile.open(aFileName);
    if (!oFile.is_open())
    {
        exit(EXIT_FAILURE);
    }
    // Write runoff data to max precision for exact conversion between binary and text
    oFile << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);

    // Include a header line
    oFile << "Year" << ",Region" << ",Basin" << ",RunoffData" << std::endl;

    // Loop through the vectors and write each entry to the file. Note that the same index across all vectors corresponds to the same region
    for (int i = 0; i < aLength; i++) 
    {
        oFile << aYears[i] << "," << aRegions[i] << "," << aBasins[i] << "," << aRunoffDataVector[i] << std::endl;
    }
    oFile.close();
}

/*!
 * \brief Read runoff data from CSV file
 * 
 * Reads basin runoff data from a comma-separated values (CSV) file and populates
 * parallel vectors with the data. This is the complement to writeRunoffData() and can
 * read files created by that method or similarly formatted CSV files.
 * 
 * \param aFileName [INPUT] Path to input CSV file to read
 * \param aYears [OUTPUT] Vector of years (cleared then populated)
 * \param aRegions [OUTPUT] Vector of region names (cleared then populated)
 * \param aBasins [OUTPUT] Vector of basin names (cleared then populated)
 * \param aRunoffDataVector [OUTPUT] Vector of runoff data (cleared then populated)
 * 
 * \return Number of data rows successfully read (equal to vector sizes)
 * 
 * \details
 * **Expected CSV Input Format:**
 * \code
 * Year,Region,Basin,RunoffData
 * 2050,USA,Mississippi,2845.342145678901234
 * 2050,China,Yangtze,3421.876543210987
 * 2050,India,Ganges,456.7890123456789
 * \endcode
 * 
 * **Column Descriptions:**
 * - **Column 1**: Year (integer, e.g., 2020, 2025, 2050)
 * - **Column 2**: Region name (string, e.g., "USA", "China", "India")
 * - **Column 3**: Basin name (string, e.g., "Mississippi", "Yangtze", "Ganges")
 * - **Column 4**: Runoff data (double, mm/year or similar units)
 * 
 * **File Processing:**
 * 
 * The method:
 * 1. Clears all output vectors to ensure clean state
 * 2. Skips the header row (first line)
 * 3. Reads each data row and parses the comma-separated fields
 * 4. Appends parsed values to corresponding vectors
 * 5. Returns the count of rows read
 * 
 * **Use Cases:**
 * - Read baseline runoff data for initialization
 * - Load runoff data from previous runs for restart
 * - Import runoff data from external sources
 * - Read archived runoff data for comparison
 * 
 * **Error Handling:**
 * - If file cannot be opened: exits program with EXIT_FAILURE
 * - No validation of data values (assumes well-formed CSV)
 * - No checking for missing or extra columns
 * 
 */
int RunoffData::readRunoffData(const std::string &aFileName, std::vector<int> &aYears, std::vector<std::string> &aRegions, 
    std::vector<std::string> &aBasins, std::vector<double> &aRunoffDataVector) 
{
    std::ifstream data(aFileName);
    if (!data.is_open())
    {
        exit(EXIT_FAILURE);
    }
    std::string line;
    // Skip the first line (i.e., the header)
    std::getline(data, line); 

    // Clear the vectors first to ensure they are empty before adding data
    aYears.clear();
    aRegions.clear();
    aBasins.clear();
    aRunoffDataVector.clear();

    while (std::getline(data, line))
    {
        std::istringstream iss(line);
        std::string token;
        int year;
        std::string region;
        std::string basin;
        double runoffData;
        
        // Parse current year
        std::getline(iss, token, ',');
        year = std::stoi(token);
        
        // Parse region
        std::getline(iss, region, ',');
        
        // Parse basin
        std::getline(iss, basin, ',');
        
        // Parse average runoff data    
        std::getline(iss, token, ',');
        runoffData = std::stod(token);
        
        aYears.push_back(year);
        aRegions.push_back(region);
        aBasins.push_back(basin);
        aRunoffDataVector.push_back(runoffData);       
    }

    data.close();
    return static_cast<int>(aYears.size());
}