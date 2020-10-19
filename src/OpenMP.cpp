/*********************************************************************************
 * 
 *                  OpenMP.cpp (Morfo70)
 *
 *  functions to manage wave parallization  (ONLY WAVE)
 *  inf NO USE OPENMP, the functions return serial values 
 *  Implementation in:
 *
 *  Methods and parameters description in "OpenMP.h"  
 * 
 * *******************************************************************************/

#include "OpenMP.h"

namespace Omp{
    int get_num_processors(){
        #if USE_OMP
            return omp_get_max_threads();
        #else
            return 1;
        #endif
    }

    int get_pid(){
        #if USE_OMP
            return omp_get_thread_num();
        #else
            return 0;
        #endif
    }
}