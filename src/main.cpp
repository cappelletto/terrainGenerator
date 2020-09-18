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
    // First, we populate the structures with default values
    // Second, we parse the config file (YAML) if any
    // Third, we parse the CLI arguments, if available

    // TODO : add get default values for both structures
    canvasParametersStruct canvas; //= getDefaultParams(); // structure to hold configuration (populated with defaults).
    fnParametersStruct     function;
    teg::useDefaults (&canvas, &function);
    // Config structures will be updated if a config file or command line arguments are provided
    YAML::Node config;
    if (argConfig)     // check if config YAML file is provided
        config = teg::readConfiguration(args::get(argConfig), &canvas, &function); // populates params structure with content of the YAML file

    // Input file priority: must be defined either by the config.yaml or --input argument
    string templateFileName = ""; // command arg or config defined
    string outputFileName   = ""; // same relative folder
    int verbosityLevel      = teg::NO_VERBOSE;
    int numThreads          = teg::N_MAX_THREAD;

    // 

    if (argTemplate)       templateFileName    = args::get(argTemplate);
    if (argOutput)         outputFileName      = args::get(argOutput);
    if (argVerbose)        verbosityLevel      = args::get(argVerbose);
    if (argNThreads)       numThreads          = args::get(argNThreads);
    if (argParamType){ // would love to have  aswitch statement but the arg is string (for user simplicity)
        string option = args::get(argParamType);
        if      (option == "constant")   function.type = teg::constant;
        else if (option == "step")       function.type = teg::step;
        else if (option == "ramp")       function.type = teg::ramp;
        else if (option == "pulse")      function.type = teg::pulse;
        else if (option == "square")     function.type = teg::square;
        else if (option == "pwm")        function.type = teg::pwm;
        else if (option == "saw")        function.type = teg::saw;
        else if (option == "triangular") function.type = teg::triangular;
        else if (option == "sine")       function.type = teg::sine;
        else // uh oh, it appears to be an invalid option... let's inform the user and fall back to the existing option
            cout << "[main] unrecognized user provide waveform type: [" << yellow << argParamType << reset <<"] " << endl;
    }
    if (argParamAmplitude) function.amplitude  = args::get(argParamAmplitude);
    if (argParamFrequency) function.frequency  = args::get(argParamFrequency);
    if (argParamPeriod)    function.period     = args::get(argParamPeriod);
    if (argParamPhase)     function.phase      = args::get(argParamPhase);
    if (argParamOffset)    function.offset     = args::get(argParamOffset);

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
