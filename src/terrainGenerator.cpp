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
#include <math.h>
#include <cmath>
#include <opencv2/core.hpp>
#include "terrainGenerator.h"

using namespace std;
using namespace teg;

void teg::generateTerrain(teg::canvasParametersStruct canvas, teg::fnParametersStruct func, cv::Mat img){
    cout << "[gt] Generating [" << canvas.cols << " x " << canvas.rows << "] map" << endl;
        // pwm,            // similar to square, but with user defined duty-cycle
    int col, row, c;
    double x, y, z, t, t1, t2;
    double rx = cos(canvas.rotation * M_PI/180);
    double ry = sin(canvas.rotation * M_PI/180);
    double r;
    switch (func.type){
        case teg::constant:
            cout << yellow << "[gt] Creating CONSTANT map [ z=" << func.offset << " ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                for (row=0; row<canvas.rows; row++){
                    z = func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::step:
            cout << yellow << "[gt] Creating STEP map [ t>0, z=" << func.offset << " ]" << reset << endl;
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

        case teg::ramp:
            cout << yellow << "[gt] Creating RAMP map [ t>0, z=t ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + row*canvas.resolution;
                    t = transform (x*rx, y*ry , func.period, func.phase); 
                    z = 0;
                    if (t>=0) z=t;
                    z = func.amplitude*z + func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::pulse:
            cout << yellow << "[gt] Creating PULSE map [ t>0, z=t ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + row*canvas.resolution;
                    t = transform (x*rx, y*ry , func.period, func.phase); 
                    z = 0;
                    if (t>=0 && t<=func.period) z=1;
                    z = func.amplitude*z + func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::sine:
            cout << yellow << "[gt] Creating SINE map [ z=sin(t) ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + row*canvas.resolution;
                    t = transform (x*rx, y*ry , func.period, func.phase); 
                    z = sin(2*M_PI*t);
                    z = func.amplitude*z + func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::square:
            cout << yellow << "[gt] Creating SQUARE map [ z=sq(t) ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + row*canvas.resolution;
                    t = transform (x*rx, y*ry , func.period, func.phase); 
                    z = 0;
                    if (sin(2*M_PI*t)>=0) z=1;
                    z = func.amplitude*z + func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::saw:
            cout << yellow << "[gt] Creating SAW map [ z=saw(t) ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + row*canvas.resolution;
                    t = transform (x*rx, y*ry , func.period, func.phase); 
                    z = t - floor(t);
                    z = func.amplitude*z + func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::triangular:
            cout << yellow << "[gt] Creating TRIANGULAR map [ z=tri(t) ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + row*canvas.resolution;
                    t = transform (x*rx, y*ry , func.period, func.phase); 
                    t1 = t - floor(t);
                    // t = transform (x*rx, y*ry , func.period, func.phase); 
                    t2 = t - floor(t + 0.5);    // we displace half a period (theta += T/2)
                    if (t1 > t2)
                        z = 2*(t1 - 0.5);
                    else
                        z = 2*(0.5 - t1);
                    z = func.amplitude*z + func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::gaussian:
            cout << yellow << "[gt] Creating GAUSSIAN map [ z=N(t,period) ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + row*canvas.resolution;
                    t = transform (x*x, y*y , func.period, func.phase); 
                    z = exp(-(t*t)/(2*func.period*func.period));
                    z = func.amplitude*z + func.offset;
                    img.at<double>(cv::Point(col,row)) = z;
                }
            }
            break;

        case teg::circle:
            cout << yellow << "[gt] Creating CIRCLE map [ z = k -(x^2 + y^2)) ]" << reset << endl;
            for (col=0; col<canvas.cols; col++){
                x = canvas.xmin + (double)col*canvas.resolution;
                for (row=0; row<canvas.rows; row++){
                    y = canvas.ymin + (double)row*canvas.resolution;
                    t = x*x + y*y;
                    r = func.period*func.period;
                    if (t>=r) z = 0;
                    else z = sqrt(r - t);
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
            else if (option == "gaussian")   func->type = teg::gaussian;
            else if (option == "circle")     func->type = teg::circle;
            else // uh oh, it appears to be an invalid option... let's inform the user and fall back to the existing option
                cout << "[readConfiguration] unrecognized YAML provide waveform type: [" << yellow << option << reset <<"] " << endl;
        }

    }

    return config;

} 