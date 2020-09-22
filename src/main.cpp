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

#include "options.h"
#include "terrainGenerator.h"
#include "ogrsf_frmts.h"
// #include "geotiff.hpp" // Geotiff class definitions

#include <geotiff.hpp>
#include <gdal_priv.h>
#include <cpl_conv.h> // for CPLMalloc()

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
    // int retval;
    int retval = initParser(argc, argv);   // initial argument validation, populates arg parsing structure args
    if (retval != 0){  // some error ocurred, we have been signaled to stop
        cout << red << "[main] Some error ocurred when parsing provided arguments" << endl;
        return retval;
    }
    if (argc<=1){
        cout << "[main] No argument provided, please use --help to show available options" << endl;
        return 0;
    }
    // Parameters hierarchy
    // ARGS > CONFIG > DEFAULT (this)
    // First, we populate the structures with default values
    // Second, we parse the config file (YAML) if any
    // Third, we parse the CLI arguments, if available

    // FIRST +++++++++++++++++++++++++++++++++++++++++++++++++++
    canvasParametersStruct canvas;      // structure to hold configuration (to be populated with defaults).
    fnParametersStruct     function;
    teg::useDefaults (&canvas, &function);
    // Config structures will be updated if a config file or command line arguments are provided
    // SECOND +++++++++++++++++++++++++++++++++++++++++++++++++++
    YAML::Node config;
    // cout << "Do we have argConfig?" << endl;
    if (argConfig)     // check if config YAML file is provided
    {
        // cout << "YES" << endl;
        config = teg::readConfiguration(args::get(argConfig), &canvas, &function); // populates params structure with content of the YAML file
    }
    // Input file priority: must be defined either by the config.yaml or --input argument
    string templateFileName = ""; // command arg or config defined
    string outputFileName   = "default.tif"; // same relative folder
    int verbosityLevel      = teg::NO_VERBOSE;
    int numThreads          = teg::N_MAX_THREAD;
    if (config["general"]) {
        if (config["general"]["verbosity"])
            verbosityLevel = config["general"]["verbosity"].as<int>(); //verbosity level
    }
    if (config["output"]) {
        if (config["output"]["filename"])
            outputFileName = config["output"]["filename"].as<string>(); //verbosity level
    }
    // THIRD +++++++++++++++++++++++++++++++++++++++++++++++++++
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
    cout << yellow << "****** Summary **********************************" << reset << endl;
    cout << "Template file:   \t" << templateFileName << endl;
    cout << "Output file:     \t" << outputFileName << endl;
    cout << "Verbosity level: \t" << verbosityLevel << endl;
    cout << "Num. threads:    \t" << numThreads << endl;
    teg::printParams(canvas, function);

    // At this point, we have combined default and user defined parameters
    // Let's proceed to create the canvas and start populating it. First, we need to determine the size of the storage unit
    // As we are planning to use geoTIFF as containers, we should point directly to its matrix structure
    canvas.cols = ceil(fabs((canvas.xmax - canvas.xmin) / canvas.resolution));
    canvas.rows = ceil(((canvas.ymax - canvas.ymin) / canvas.resolution));
    cout << "[main] Canvas size: [ " << canvas.cols << " x " << canvas.rows << "]" << endl;
    // this is the cvMat container
    cv::Mat rasterData(canvas.rows, canvas.cols, CV_64FC1, canvas.nodata);

    //let's populate the image with the corresponding map
    teg::generateTerrain (canvas, function, rasterData);

    GDALAllRegister();
    GDALDataset     *geotiffDataset;
    GDALDriver      *driverGeotiff;
    GDALRasterBand  *geotiffBand; // also declare pointers for Geotiff

    int nrows  = rasterData.rows; //layerDimensions[1]
    int ncols  = rasterData.cols; //layerDimensions[0]

    cout << "[main] Dataset dimensions (COL x ROW): [" << ncols << "] x [" << nrows << "]\tNoData = [" << canvas.nodata << "]" << endl; 

    double transformMatrix[6];      // 6-element geotranform array.
    transformMatrix[GEOTIFF_PARAM_SX] = canvas.resolution;
    transformMatrix[GEOTIFF_PARAM_SY] = canvas.resolution;
    transformMatrix[GEOTIFF_PARAM_CX] = canvas.xmin;
    transformMatrix[GEOTIFF_PARAM_CY] = canvas.ymin;

    char **optionsForTIFF = NULL;
    optionsForTIFF = CSLSetNameValue(optionsForTIFF, "COMPRESS", "LZW");
    driverGeotiff = GetGDALDriverManager()->GetDriverByName("GTiff");
    geotiffDataset = driverGeotiff->Create(outputFileName.c_str(), ncols, nrows, 1, GDT_Float64, optionsForTIFF);
    geotiffDataset->SetGeoTransform(transformMatrix);
        // // cout << "[r.writeLayer] Projection string:" << endl;
        // // cout << layerProjection.c_str() << endl;
        // geotiffDataset->SetProjection(layerProjection.c_str());
        // // \todo figure out if we need to convert/cast the cvMat to float/double for all layers
    int errcode;
    double *rowBuff = (double*) CPLMalloc(sizeof(double)*canvas.cols);
    geotiffDataset->GetRasterBand(1)->SetNoDataValue (canvas.nodata);       
    for(int row=0; row<canvas.rows; row++) {
        for(int col=0; col<canvas.cols; col++) {
            rowBuff[col] = (double) rasterData.at<double>(cv::Point(col,row)); // tempData should be CV_64F
        }
        errcode = geotiffDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, row, canvas.cols, 1, rowBuff, canvas.cols, 1, GDT_Float64, 0, 0);
    }

    GDALClose(geotiffDataset) ;
    return 0;
}
