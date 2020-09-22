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
#include <opencv2/core.hpp>
#include "terrainGenerator.h"

using namespace std;
using namespace teg;

void teg::generateTerrain(teg::canvasParametersStruct canvas, teg::fnParametersStruct func, cv::Mat img){
    cout << "[gt] Generating [" << canvas.cols << " x " << canvas.rows << "] map" << endl;
        // constant = 1,   // fixed value y = D
        // step,           // unitary step response, y=1 if t>0
        // ramp,           // slope=1 ramp defined for t>0
        // pulse,          // square pulse defined for a given period 0<= t <= T, where T:1/B
        // square,         // periodic squared signal 
        // pwm,            // similar to square, but with user defined duty-cycle
        // triangular,     // periodic triangular signal
        // saw,            // saw-shaped periodic signal
        // sine            // pure sinusoidal
    int col, row;
    double x, y, z, t;
    double rx = cos(canvas.rotation);
    double ry = sin(canvas.rotation);

    cout << "rx " << rx << endl;
    cout << "ry " << ry << endl;


    switch (func.type){
        case teg::constant:
            cout << yellow << "[gt] Creating constant map [ z=" << func.offset << " ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                for (row=0; row<canvas.rows; row++){
                    z = func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::step:

            cout << yellow << "[gt] Creating step map [ t>0, z=1" << func.offset << " ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + row*canvas.resolution;
                    t = transform (x*rx, y*ry , func.period, func.phase); 
                    z = 0;
                    if (t>=0) z=1;
                    z = func.amplitude*z + func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        default:
            cout << red << "\t Unknown terrain type specified [" << func.type << "]" << reset << endl; 
            break;
    }
}

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
        if (config["parameters"]["waveform"]){
            string option = config["parameters"]["waveform"].as<string>();
            if      (option == "constant")   func->type = teg::constant;
            else if (option == "step")       func->type = teg::step;
            else if (option == "ramp")       func->type = teg::ramp;
            else if (option == "pulse")      func->type = teg::pulse;
            else if (option == "square")     func->type = teg::square;
            else if (option == "pwm")        func->type = teg::pwm;
            else if (option == "saw")        func->type = teg::saw;
            else if (option == "triangular") func->type = teg::triangular;
            else if (option == "sine")       func->type = teg::sine;
            else // uh oh, it appears to be an invalid option... let's inform the user and fall back to the existing option
                cout << "[readConfiguration] unrecognized YAML provide waveform type: [" << yellow << option << reset <<"] " << endl;
        }

    }

    return config;

} 