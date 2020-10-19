/*********************************************************************************
 * 
 *                  OpenMP.h (Morfo70)
 * 
 *  functions to manage wave parallization  (ONLY WAVE)
 *  inf NO USE OPENMP, the functions return serial values 
 *  Implementation in:
 *      - OpenMP.cpp 
 * 
 * *******************************************************************************/
#ifndef OPENMP_H_
#define OPENMP_H_

#include <omp.h>

//OpenMP flaf
#ifndef USE_OMP 
    #define USE_OMP 0
#endif

namespace Omp{
    /**
     *  Return number of system processors if OMP flag is turned on, otherwise return 1
     * */    
    extern int get_num_processors();

    /**
     *  Return current proccesor id if OMP flag is turned on, otherwise return 0
     * */    
    extern int get_pid();
}


#endif