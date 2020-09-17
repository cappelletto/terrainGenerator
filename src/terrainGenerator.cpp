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
// #include "options.h"

using namespace std;
using namespace teg;

/**
 * @brief 
 * 
 * @param filename 
 * @param param 
 * @param func 
 * @return YAML::Node 
 */
YAML::Node teg::readConfiguration(string filename, canvasParametersStruct *canvas, fnParametersStruct *func){

    cout << "Processing user defined configuration: [" <<  filename  << "]" << endl;
    YAML::Node config = YAML::LoadFile(filename);

    int verb = 0;
    if (config["general"]) {
        if (config["general"]["verbosity"])
            verb = config["general"]["verbosity"].as<int>(); //verbosity level
    }

    if (config["template"]){
        if (verb > teg::NO_VERBOSE)
            cout << "[readConfiguration] template section present" << endl;
        if (config["template"]["filename"])
            cout << "Template file: " << config["template"]["filename"].as<string>() << endl;
    }

    if (config["canvas"]){
        if (verb > 0)
            cout << "[readConfiguration] canvas section present" << endl;
        if (config["canvas"]["resolution"])
            canvas->resolution = config["canvas"]["resolution"].as<double>();
        if (config["canvas"]["rotation"])
            canvas->rotation = config["canvas"]["rotation"].as<double>();
        if (config["canvas"]["mode"])
            canvas->mode = config["canvas"]["mode"].as<int>();
        if (config["canvas"]["nodata"])
            canvas->nodata = config["canvas"]["nodata"].as<double>();
    }

    if (config["parameters"]){
        if (verb > 0)
            cout << "[readConfiguration] parameters section present" << endl;
        if (config["parameters"]["amplitude"])
            func->amplitude =   config["parameters"]["amplitude"].as<double>();
        if (config["parameters"]["frequency"])
            func->frequency =   config["parameters"]["frequency"].as<double>();
        if (config["parameters"]["phase"])
            func->phase =       config["parameters"]["phase"].as<double>();
        if (config["parameters"]["offset"])
            func->offset =      config["parameters"]["offset"].as<double>();
        if (config["parameters"]["period"])
            func->period =      config["parameters"]["period"].as<double>();
    }

    return config;

} 