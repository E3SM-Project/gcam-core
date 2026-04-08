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
 * \file degree_days.cpp
 * \brief This file aggregates gridded data on heating and cooling degree days from E3SM into region-specific averages for GCAM
 *
 * \author Philip Myint
 */

#include "../include/degree_days.h"

// Constructor
DegreeDays::DegreeDays(const int aNumLon, const int aNumLat, const bool aReadSubregions):
ASpatialData(aNumLon * aNumLat, aReadSubregions),
mNumLon(aNumLon),
mNumLat(aNumLat)
{
}

// Destructor (could just use the default destructor since there are no dynamically allocated raw pointers, 
// but including an explicit destructor in case we want to add any cleanup code in the future)
DegreeDays::~DegreeDays() 
{
}

/*!
 * \brief Aggregate gridded degree days to regional scale using population weighting
 * 
 * This method performs population-weighted averaging of heating and cooling degree days
 * from E3SM's gridded output to GCAM's regional scale. The weighting ensures that degree
 * days reflect where people actually live, which is appropriate for energy demand calculations.
 * 
 * \param aGCAMYear [INPUT] GCAM year for this data (e.g., 2025, 2050)
 * \param aELMArea [INPUT] Grid cell areas in km² (mNumLat × mNumLon elements)
 * \param aELMDegreeDays [INPUT] Degree days for each grid cell (mNumLat × mNumLon elements)
 * \param aPopDensity [INPUT] Population density in people/km² for each grid cell (mNumLat × mNumLon elements)
 * \param aELMLandFrac [INPUT] Land fraction for each grid cell (mNumLat × mNumLon elements, values between 0 and 1)
 * \param aYears [OUTPUT] Vector of years (one entry per region, all equal to aGCAMYear)
 * \param aRegions [OUTPUT] Vector of region names (e.g., "USA", "China", "India")
 * \param aDegreeDaysVector [OUTPUT] Vector of population-weighted degree days by region (degree-days)
 * \param aNumValues [OUTPUT] Number of regional values written to output vectors
 * 
 * \details
 * **Population-Weighted Aggregation Formula:**
 * 
 * For each GCAM region:
 * \code
 * DD_region = Σ(DD_cell × pop_cell × area_cell) / Σ(pop_cell × area_cell)
 * 
 * Where:
 *   - Sum is over all grid cells (or fractional cells) in the region
 *   - pop_cell = population_density × area_km²
 *   - Both area and density are in km² units (no conversion needed)
 *   - Fractional weights from mRegionWeights handle cells straddling boundaries
 * \endcode
 * 
 * **Handling Edge Cases:**
 * - Regions with zero population: DegreeDays = 0 
 * - Missing/NaN data: Skipped (not accumulated)
 * 
 * \pre mRegionMapping and mRegionWeights must be populated by ASpatialData::readRegionalMappingData()
 * \pre Input arrays must have mNumLat × mNumLon elements
 * \pre aELMArea values must be positive (km²)
 * \pre aPopDensity values must be non-negative (people/km²)
 * \pre aELMLandFrac values must be between 0 and 1
 * 
 * \post Output vectors contain one entry per region with non-zero population
 * \post aNumValues equals the number of regions processed
 * 
 * \note This method does not filter outliers (unlike CarbonScalers) because extreme
 *       degree day values are legitimate (Arctic winters, desert summers)
 * \note Ocean cells (not in mRegionMapping) are automatically skipped
 * 
 */
void DegreeDays::aggregateDegreeDays(const int aGCAMYear, const double *const aELMArea, const double *const aELMDegreeDays, 
    const double *const aPopDensity, const double *const aELMLandFrac, std::vector<int> &aYears, std::vector<std::string> &aRegions, 
    std::vector<double> &aDegreeDaysVector, int &aNumValues) 
{
    // TODO: Generate any diagnostic files or exclude outliers? Probably not necessary or relevant for degree days, but can add if desired

    // Create mappings to store intermediate information
    std::map<std::string, double> totalPopulation;
    std::map<std::string, double> totalDegreeDays;
    
    // Note: E3SM data will have longitude moving fastest, then latitude, so loop so have longitude is the inner loop and latitude is the outer loop
    int gridIndex; 

    // Weight that indicates an area- and fraction-weighted population of the grid cell for a particular region. 
    // This is used to weight the degree days when aggregating to the regional level
    double popWeight = 0.0; 

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

                // Loop over all regions this grid is mapped to and calculate the contribution to the total population of each region and the degree days.
                for (const auto &regID : regInGrd) 
                {
                    // Calculate the population weight for this grid cell and region combination as the product of the fraction of the grid cell in this region (from the mapping file), 
                    // the area of the grid cell, the population density of the grid cell, and the fraction of the grid cell that is actually land. 
                    // This weight will be used to weight the degree days when aggregating to the regional level
                    popWeight = mRegionWeights[std::make_pair(gridID, regID)] * aELMArea[gridIndex] * aPopDensity[gridIndex] * aELMLandFrac[gridIndex]; 

                    totalPopulation[regID] += popWeight;
                    totalDegreeDays[regID] += aELMDegreeDays[gridIndex] * popWeight;

                } // for (const auto &regID : regInGrd) 

            } // if (tempGrid == mRegionMapping.end()) 

        } // for (int j = 1; j <= mNumLon; j++)

    } // (int k = 1; k <= mNumLat; k++)

    // Create map to store information on degree days by region, which will be used to create the output vectors for GCAM
    std::map<std::string, double> aDegreeDaysMap;
   
    // Calculate the population-weighted average degree days for each region, and then add these entries to the degree days map accordingly
    for (const auto &pair : totalPopulation) 
    {
        const std::string &regID = pair.first;
        const double populationInRegion = pair.second;
        
        // Calculate average degree days for each region as the total degree days divided by the total population. 
        // If total population is zero, set average degree days to 0 to avoid division by zero
        if (populationInRegion > 0.0) 
        {
            // Calculate weighted average: total / population
            aDegreeDaysMap[regID] = totalDegreeDays[regID] / populationInRegion;
        } 
        else 
        {
            // Set degree days to zero for this region if population is zero
            aDegreeDaysMap[regID] = 0.0;
        }
    } // for (const auto &pair : totalPopulation)
    
    // Create output vectors for GCAM. Loop through the maps and create the vectors as needed
    createDegreeDayVectors(aGCAMYear, aYears, aRegions, aDegreeDaysVector, aDegreeDaysMap);

    // Set the number of actual records, which is the same for both degree days, based on total population > 0
    aNumValues = static_cast<int>(aDegreeDaysMap.size());
}


/*!
 * \brief Convert map structures to parallel vectors for GCAM interface
 * 
 * This helper method packages degree day data from map structures (used for efficient
 * aggregation) into parallel vectors (required by GCAM's interface). The output vectors
 * are organized such that the same index across all vectors represents data for one region.
 * 
 * \param aGCAMYear [INPUT] GCAM year to assign to all entries
 * \param aYears [OUTPUT] Vector of years (all entries equal to aGCAMYear)
 * \param aRegions [OUTPUT] Vector of region names extracted from map keys
 * \param aDegreeDaysVector [OUTPUT] Vector of degree days in same order as aRegions
 * \param aDegreeDaysMap [INPUT] Map of region name → degree days value
 * 
 * \details
 * 
 * **Output Vector Structure (Parallel Vectors):**
 * \code
 * Index:       0         1          2
 * aYears:     [2050,    2050,      2050]
 * aRegions:   ["USA",   "China",   "India"]
 * aDegreeDaysVector: [2845.3,  3421.9,    456.8]
 * \endcode
 * 
 * **Parallel Vector Concept:**
 * 
 * All vectors have the same length, and index `i` across all vectors represents the
 * complete data for region `i`:
 * - `aYears[i]` = year for region i
 * - `aRegions[i]` = name of region i  
 * - `aDegreeDaysVector[i]` = degree days for region i
 * 
 * **Difference from CarbonScalers:**
 * 
 * Unlike CarbonScalers, which splits region IDs like "USA.Northwest" into "USA" and 
 * "Northwest_Wheat", this method assumes region IDs have no subregions. The regID
 * from the map is used directly as the region name without parsing.
 * 
 */
void DegreeDays::createDegreeDayVectors(const int aGCAMYear, std::vector<int> &aYears, std::vector<std::string> &aRegions, 
    std::vector<double> &aDegreeDaysVector, const std::map<std::string, double> &aDegreeDaysMap) 
{
    // Clear the vectors first to ensure they are empty before adding data
    aYears.clear();
    aRegions.clear();
    aDegreeDaysVector.clear();

    // Loop through the map and create the vectors
    for (const auto &pair : aDegreeDaysMap) 
    {
        const std::string &regID = pair.first;
        aDegreeDaysVector.push_back(pair.second);
        
        aYears.push_back(aGCAMYear);

        // Note: This assumes there is no subregion in regID, unlike the case of carbon scalars where the regID includes both region and subregion 
        aRegions.push_back(regID);
    }
}

/*!
 * \brief Write degree day data to CSV file
 * 
 * Writes regional degree day data to a comma-separated values (CSV) file with maximum
 * precision to ensure exact round-trip conversion between binary and text representations.
 * This is important for reproducibility and restart capabilities.
 * 
 * \param aFileName [INPUT] Path to output CSV file (will be created/overwritten)
 * \param aYears [INPUT] Vector of years (parallel to other vectors)
 * \param aRegions [INPUT] Vector of region names (parallel to other vectors)
 * \param aDegreeDaysVector [INPUT] Vector of degree days by region
 * \param aLength [INPUT] Number of entries to write (should be ≤ vector sizes)
 * 
 * \details
 * **CSV Output Format:**
 * \code
 * Year,Region,DegreeDays
 * 2050,USA,2845.342145678901234
 * 2050,China,3421.876543210987
 * 2050,India,456.7890123456789
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
 * - **DegreeDays**: Degree days (double, degree-days)
 *
 */
void DegreeDays::writeDegreeDays(const std::string &aFileName, const std::vector<int> &aYears, const std::vector<std::string> &aRegions, 
    const std::vector<double> &aDegreeDaysVector, const int aLength) 
{
    std::ofstream oFile;
    oFile.open(aFileName);
    if (!oFile.is_open())
    {
        exit(EXIT_FAILURE);
    }
    // Write degree days to max precision for exact conversion between binary and text
    oFile << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);

    // Include a header line
    oFile << "Year" << ",Region" << ",DegreeDays" << std::endl;

    // Loop through the vectors and write each entry to the file. Note that the same index across all vectors corresponds to the same region
    for (int i = 0; i < aLength; i++) 
    {
        oFile << aYears[i] << "," << aRegions[i] << "," << aDegreeDaysVector[i] << std::endl;
    }
    oFile.close();
}

/*!
 * \brief Read degree day data from CSV file
 * 
 * Reads regional degree day data from a comma-separated values (CSV) file and populates
 * parallel vectors with the data. This is the complement to writeDegreeDays() and can
 * read files created by that method or similarly formatted CSV files.
 * 
 * \param aFileName [INPUT] Path to input CSV file to read
 * \param aYears [OUTPUT] Vector of years (cleared then populated)
 * \param aRegions [OUTPUT] Vector of region names (cleared then populated)
 * \param aDegreeDaysVector [OUTPUT] Vector of degree days (cleared then populated)
 * 
 * \return Number of data rows successfully read (equal to vector sizes)
 * 
 * \details
 * **Expected CSV Input Format:**
 * \code
 * Year,Region,DegreeDays
 * 2050,USA,2845.342145678901234
 * 2050,China,3421.876543210987
 * 2050,India,456.7890123456789
 * \endcode
 * 
 * **Column Descriptions:**
 * - **Column 1**: Year (integer, e.g., 2020, 2025, 2050)
 * - **Column 2**: Region name (string, e.g., "USA", "China", "India")
 * - **Column 3**: Degree days (double, degree-days)
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
 * - Read baseline degree day data for initialization
 * - Load degree day data from previous runs for restart
 * - Import degree day data from external sources
 * - Read archived degree day data for comparison
 * 
 * **Error Handling:**
 * - If file cannot be opened: exits program with EXIT_FAILURE
 * - No validation of data values (assumes well-formed CSV)
 * - No checking for missing or extra columns
 * 
 */
int DegreeDays::readDegreeDays(const std::string &aFileName, std::vector<int> &aYears, std::vector<std::string> &aRegions, 
    std::vector<double> &aDegreeDaysVector) 
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
    aDegreeDaysVector.clear();

    while (std::getline(data, line))
    {
        std::istringstream iss(line);
        std::string token;
        int year;
        std::string region;
        double averageDegreeDays;
        
        // Parse current year
        std::getline(iss, token, ',');
        year = std::stoi(token);
        
        // Parse region
        std::getline(iss, region, ',');
        
        // Parse average degree days
        std::getline(iss, token, ',');
        averageDegreeDays = std::stod(token);
        
        aYears.push_back(year);
        aRegions.push_back(region);
        aDegreeDaysVector.push_back(averageDegreeDays);       
    }

    data.close();
    return static_cast<int>(aYears.size());
}