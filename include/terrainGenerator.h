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

// #include "options.h"

#include <yaml-cpp/yaml.h>

namespace teg{

    //constant definitions
    const int N_MAX_THREAD = 24;
    
    enum VerbosityLevel{
        NO_VERBOSE = 0,
        LOW_VERBOSITY,
        HIGH_VERBOSITY
    };

    enum CanvasType{
        depthMap = -1,  // the map is positive down (standard bathymetry)
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
        double resolution;  // meter per pixel
        double rotation;    // default : 0
        double nodata;      // for geoTIFF export
    }canvasParametersStruct;


    YAML::Node readConfiguration(std::string file, teg::canvasParametersStruct*, teg::fnParametersStruct*); // populates params structure with content of the YAML file


};


#endif //_TERRAIN_GENERATOR_H_