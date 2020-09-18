/**
 * @file options.h
 * @brief Argument parser options based on args.hxx
 * @version 1.1
 * @date 18/06/2020
 * @author Jose Cappelletto
 */

#ifndef _TEG_OPTIONS_H_
#define _TEG_OPTIONS_H_

#include <opencv2/core.hpp> // for the OpenCV version info
#include <args.hxx>
#include <iostream>
#include "terrainGenerator.h"

using namespace std;
using namespace teg;

args::ArgumentParser argParser("","");
args::HelpFlag 	     argHelp(argParser, "help", "Display this help menu", {'h', "help"});
// args::CompletionFlag completion(argParser, {"complete"});	//TODO: figure out why is missing in current version of args.hxx

args::ValueFlag <std::string> 	argTemplate(argParser,     "template", "Provides a file that will be used as template for the canvas definition, typ. a geoTIFF", {"template"});
args::ValueFlag	<std::string> 	argOutput(argParser,    "output",   "Output file",{'o',"output"});
args::ValueFlag	<int> 	        argVerbose(argParser,   "level",  "Define verbosity level",          {"verbose"});
args::ValueFlag	<int> 	        argNThreads(argParser,  "number",   "Define max number of threads",  {"nthreads"});
args::ValueFlag <std::string>   argConfig(argParser,    "file.yaml","Provides path to file with user defied configuration", {"config"});

// Free parameters for debugging
args::ValueFlag	<int> 	argIntParam(argParser,  "param",    "User defined parameter INTEGER for testing purposes",  {"int"});
args::ValueFlag	<float> argFloatParam(argParser,"param",    "User defined parameter FLOAT for testing purposes",    {"float"});

// Waveform parameters
args::ValueFlag	<string> argParamType      (argParser,"type",    "Type of waveform: constant, step, pulse, ramp, square, pwm, triangular, saw, sine", {"type"});
args::ValueFlag	<double> argParamPeriod    (argParser,"period",  "Waveform spatial period. This overrides the frequency when both are provided", {"period"});
args::ValueFlag	<double> argParamFrequency (argParser,"frequency",  "Waveform spatial frequency, reciprocal definition of the period", {"frequency"});
args::ValueFlag	<double> argParamAmplitude (argParser,"peak",    "Waveform peak amplitude", {"amplitude"});
args::ValueFlag	<double> argParamOffset    (argParser,"offset",  "Offset in the vertizal (Z) axis", {"offset"});
args::ValueFlag	<double> argParamPhase     (argParser,"phase",   "Waveform offset in the horizontal axis (XY)", {"phase"});


// int initParser(int argc, char *argv[]);
int initParser(int argc, char *argv[]){
        //*********************************************************************************
    /* PARSER section */
    std::string descriptionString =
        "terrainGenerator - testing module part of [landing-area-detection] pipeline \
         Compatible interface with geoTIFF bathymetry datasets via GDAL + OpenCV";

    argParser.Description(descriptionString);
    argParser.Epilog("Author: J. Cappelletto (GitHub: @cappelletto)\n");
    argParser.Prog(argv[0]);
    argParser.helpParams.width = 120;

    try
    {
        argParser.ParseCLI(argc, argv);
    }
    catch (const args::Completion &e)
    {
        cout << e.what();
        cout << "tried to autocomplete" << endl;
        return 0;
    }

    catch (args::Help)
    { // if argument asking for help, show this message
        cout << argParser;
        return -1;
    }
    catch (args::ParseError e)
    { //if some error ocurr while parsing, show summary
        std::cerr << e.what() << std::endl;
        std::cerr << "Use -h, --help command to see usage" << std::endl;
        return -1;
    }
    catch (args::ValidationError e)
    { // if some error at argument validation, show
        std::cerr << "Bad input commands" << std::endl;
        std::cerr << "Use -h, --help command to see usage" << std::endl;
        return -1;
    }
    cout << cyan << "terrainGenerator" << reset << endl; // CREATE OUTPUT TEMPLATE STRING
    cout << "\tOpenCV version:\t" << yellow << CV_VERSION << reset << endl;
    cout << "\tGit commit:\t" << yellow << GIT_COMMIT << reset << endl
         << endl;
    // cout << "\tBuilt:\t" << __DATE__ << " - " << __TIME__ << endl;   // TODO: solve, make is complaining about this
    return 0;
}

#endif //_PROJECT_OPTIONS_H_