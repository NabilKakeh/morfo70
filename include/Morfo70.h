/*********************************************************************************
 * 
 *                  Morfo70.h (Morfo70)
 * 
 *  Global functions  and variables in order to perform simulation 
 *  Implementation in:
 *      - Morfo70.cpp 
 * 
 * *******************************************************************************/

#ifndef MORFO70_H_
#define MORFO70_H_

//Linux/Windows Path separator
#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif



#include "Sediment.h"
#include "Math.h"
#include "MorfoIniFile.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <string>


namespace Morfo70{

//////////////////////////
/// Offshore wave values  
//////////////////////////  
extern double H_off;
extern double T_off;
extern double Theta_off;
extern double Zs_off;

extern string inputDir;
extern string outputDir;
extern list<string> iniFiles;
extern string batFile;

/**
 *  Run morfo70 simulation 
 * @param outputDir folder wich stores output files
 * */
extern void Run(Grid* grid,Bathymetry* bathy,Wave* wave, Hydro* hydro, Sediment* sediment,string outputDir);

/**
 *   Save time step in netcdf file, file mask -> m70_<time>.nc
 * */
void saveTimeStep(string outputFile,double time,Grid* grid,Bathymetry* bathy,  Wave* wave, Hydro* hydro,Sediment* sediment);

/**
 *  Search in the the oputput folder for previous simulations, load data and build simulation objects 
 *  in otder to restart simulation 
 * */
double resume(string outputDir,Bathymetry* bathy,  Wave* wave, Hydro* hydro,Sediment* sediment);



/////////////////////////////////////////////////////////////////////////////////////
////  I/O Auxiliary functions
/////////////////////////////////////////////////////////////////////////////////////
/**
 * Chech if a file exists
 * */
inline bool checkFile(string path){
    struct stat info;    
     if(stat( path.c_str(), &info ) != 0)
        return 0;
    else if(info.st_mode & S_IFDIR)
        return 0;
    else
        return 1;
}

/**
 * Chech if a folder exists
 * */
inline bool checkDir(string path){
    struct stat info;
    
     if(stat( path.c_str(), &info ) != 0)
        return 0;
    else if(info.st_mode & S_IFDIR)
        return 1;
    else
        return 0;
}

/**
 * Read promt command arguments (input dir, output dir, ini file adn bathy file), and check them
 * */
extern bool readAndCheckParameters( int argc, const char** argv );
}

/////////////////////////////////////////////////////////////////////////////////////
////  END I/O Auxiliary functions
/////////////////////////////////////////////////////////////////////////////////////



#endif