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

using namespace std;
using namespace teg;

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 * @param f 
 * @return double 
 */
double teg::parametricTransform (double x, double y, teg::scaleParamStruct f){
    return (f.ax*x + f.ay*y + f.c);
}


/**
 * @brief 
 * 
 * @param canvas 
 * @param func 
 */
void teg::printParams(teg::canvasParametersStruct  canvas, teg::fnParametersStruct  func){
    //CANVAS
    cout << cyan << "Canvas: " << reset << endl;
    cout << "\tmode:      " << ((canvas.mode == teg::depthMap) ? "depthMap" : "heightMap") << endl;
    cout << "\tresolution:" << canvas.resolution   << "\t[m/pixel]" << endl;
    cout << "\trotation:  " << canvas.rotation     << "\t[degrees]" << endl;
    cout << "\tnodata:    " << canvas.nodata       << endl;
    cout << cyan << "\tExtent:" << reset <<"\tX=[ " << canvas.xmin << ", " << canvas.xmax << " ]\tY=[ " <<  canvas.ymin << ", " << canvas.ymax << " ]" << endl;

    cout << cyan << "Function: " << reset << endl;
    cout << "\ttype:       " << func.type << endl;
    cout << "\tamplitude:  " << func.amplitude  << "\t[m]" << endl;
    cout << "\tperiod:     " << func.period     << "\t[m]" << endl;
    cout << "\tphase:      " << func.phase      << "\t[m]" << endl;
    cout << "\toffset:     " << func.offset     << "\t[m]" << endl;
    cout << "\tfrequency:  " << func.frequency  << "\t[1/m]" << endl;

    cout << endl;
}

/**
 * @brief 
 * 
 * @param canvas 
 * @param func 
 */
void teg::useDefaults(teg::canvasParametersStruct *canvas, teg::fnParametersStruct *func){
    canvas->mode = teg::heightMap; //positive up
    canvas->nodata = -9999;
    canvas->resolution = 0.01;  //10 mm per pixel
    canvas->rotation = 0; // no rotation
    canvas->xmin = 0;
    canvas->ymin = 0;
    canvas->xmax = 1;
    canvas->ymax = 1;   // 1m x 1m rectangle

    func->type = teg::constant;
    func->amplitude = 1.0;
    func->frequency = 1.0;
    func->offset    = 0.0;
    func->period    = 1.0;
    func->phase     = 0.0;
}

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
        if (config["canvas"]["extent"]){
            canvas->xmin = config["canvas"]["extent"]["xmin"].as<double>();
            canvas->xmax = config["canvas"]["extent"]["xmax"].as<double>();
            canvas->ymin = config["canvas"]["extent"]["ymin"].as<double>();
            canvas->ymax = config["canvas"]["extent"]["ymax"].as<double>();
        }
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