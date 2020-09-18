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

#ifndef _TERRAIN_GENERATOR_H_
#define _TERRAIN_GENERATOR_H_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <yaml-cpp/yaml.h>

namespace teg{

    //constant definitions
    const int N_MAX_THREAD = 24;
    
    const std::string green("\033[1;32m");
    const std::string yellow("\033[1;33m");
    const std::string cyan("\033[1;36m");
    const std::string red("\033[1;31m");
    const std::string reset("\033[0m");
    const std::string highlight("\033[30;43m");

    enum VerbosityLevel{
        NO_VERBOSE = 0,
        LOW_VERBOSITY,
        HIGH_VERBOSITY
    };

    enum CanvasType{
        depthMap =  0,  // the map is positive down (standard bathymetry)
        heightMap = 1   // the map is positive up (standard altimetry)
    };

    enum WaveFormType{
        constant = 1,   // fixed value y = D
        step,           // unitary step response, y=1 if t>0
        ramp,           // slope=1 ramp defined for t>0
        pulse,          // square pulse defined for a given period 0<= t <= T, where T:1/B
        square,         // periodic squared signal 
        pwm,            // similar to square, but with user defined duty-cycle
        triangular,     // periodic triangular signal
        saw,            // saw-shaped periodic signal
        sine            // pure sinusoidal
    };

    // configuration parameters
    typedef struct{         // Waveform parameters as defined by template:
                            // Y = A.func(B.x + C) + D : amplitude x func(frequency.t + phase) + offset 
        WaveFormType type;           // waveform type
        double amplitude;   // (A)   peak amplitude of the waveform
        double period;      // (T:1/B) Spatial period of the wavefor if periodic. Signal duration if non periodic 
        double phase;       // (C)   Spatial phase (offset)
        double offset;      // (D)   Vertical offset (ordinate axis)
        double frequency;   // (B)   Reciprocal of period
    }fnParametersStruct;

    typedef struct{
        int mode;           // +1 Height, -1 Depth
        double resolution;  // as meters per pixel
        double rotation;    // default : 0
        double nodata;      // for geoTIFF export
        double xmin;        // boundariers in meters of the map canvas
        double xmax;
        double ymin;
        double ymax;        // no min < max validation is performed (this allows reverse-definition of the canvas)
    }canvasParametersStruct;

    void printParams(teg::canvasParametersStruct , teg::fnParametersStruct ); // populate paras structure with default values
    void useDefaults(teg::canvasParametersStruct*, teg::fnParametersStruct*); // populate paras structure with default values
    YAML::Node readConfiguration(std::string file, teg::canvasParametersStruct*, teg::fnParametersStruct*); // populates params structure with content of the YAML file


};


#endif //_TERRAIN_GENERATOR_H_