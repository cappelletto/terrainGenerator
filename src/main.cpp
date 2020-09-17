/**
 * @file terrainGenerator.cpp
 * @author Jose Cappelletto (cappelletto@gmail.com)
 * @brief Synthetic terrain generator based on Geotiff, OpenCV, CGAL & GDAL. Creates basic regular bathymetry/surface maps and export them as geoTIFF
 * @version 0.1
 * @date 2020-09-14
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#include "terrainGenerator.h"
#include "options.h"
#include "geotiff.hpp" // Geotiff class definitions

#include <opencv2/core.hpp>
#include <yaml-cpp/yaml.h>

using namespace std;
using namespace cv;
using namespace teg;

/*!
	@fn		int main(int argc, char* argv[])
	@brief	Main function
*/
int main(int argc, char *argv[])
{
    int retval = initParser(argc, argv);   // initial argument validation, populates arg parsing structure args
    if (retval != 0)  // some error ocurred, we have been signaled to stop
        return retval;
    // Parameters hierarchy
    // ARGS > CONFIG > DEFAULT (this)
    canvasParametersStruct canvas; //= getDefaultParams(); // structure to hold configuration (populated with defaults).
    fnParametersStruct     function;
    // Config structures will be updated if a config file or command line arguments are provided
    YAML::Node config;
    if (argConfig)     // check if config YAML file is provided
        config = teg::readConfiguration(args::get(argConfig), &canvas, &function); // populates params structure with content of the YAML file

    // Input file priority: must be defined either by the config.yaml or --input argument
    string inputFileName    = ""; // command arg or config defined
    string inputFilePath    = ""; // can be retrieved from the fully qualified inputFileName 
    string outputFilePrefix = ""; // none, output filenames will be the same as the standard
    string outputFilePath   = ""; // same relative folder

    // // override defaults or config file with command provided values (DEFAULT < CONFIG < ARGUMENT)
    // if (argAlphaRadius)     params.alphaShapeRadius = args::get(argAlphaRadius);
    // if (argGroundThreshold) params.groundThreshold  = args::get(argGroundThreshold);
    // if (argHeightThreshold) params.heightThreshold  = args::get(argHeightThreshold);
    // if (argSlopeThreshold)  params.slopeThreshold   = args::get(argSlopeThreshold);
    // if (argRobotHeight)     params.robotHeight      = args::get(argRobotHeight);
    // if (argRobotLength)     params.robotLength      = args::get(argRobotLength);
    // if (argRobotWidth)      params.robotWidth       = args::get(argRobotWidth);
    // if (argProtrusionSize)  params.protrusionSize   = args::get(argProtrusionSize);
    // if (argRotation){
    //                         params.rotation         = args::get(argRotation);
    //                         params.fixRotation      = true;
    // }   
    //**************************************************************************
    /* Summary list parameters */
    // cout << yellow << "****** Summary **********************************" << reset << endl;
    // cout << "Input file:   \t" << inputFileName << endl;
    // cout << "Input path:   \t" << inputFilePath << endl;
    // cout << "Output prefix:\t" << outputFilePrefix << endl;
    // cout << "Output path:  \t" << outputFilePath << endl;
    // cout << "fParam:       \t" << fParam << endl;
    // cout << "iParam:       \t" << iParam << endl;
    // params.robotDiagonal = sqrt(params.robotWidth*params.robotWidth + params.robotLength*params.robotLength); 
    // lad::printParams(&params);

    // cout << "Verbose level:\t\t" << pipeline.verbosity << endl;    
    // cout << "Multithreaded version, max concurrent threads: [" << yellow << nThreads << reset << "]" << endl;
    // cout << yellow << "*************************************************" << reset << endl << endl;

    // waitKey(0);
    return 0;
}
