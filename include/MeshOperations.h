/*********************************************************************************
 * 
 *                  MeshOperations.h (Morfo70)
 * 
 *  Inline function to performs matrix operations 
 * 
 * *******************************************************************************/

#ifndef MESH_OPERATIONS_H_
#define MESH_OPERATIONS_H_
#include <iostream>
namespace MeshOperations{
/**
 *  Interpolates values in center nodes to x-nodes
 * @param var_x[out] values in x-nodes
 * @param var_c values in c-nodes 
*/
inline void C2X(double  __restrict_arr  * var_x, double const __restrict_arr*  var_c,int nx, int ny)
{
   
    for (int i=0;i<nx-1;i++)
    {
             for (int j=0;j<ny;j++)
                var_x[j+i*ny]=0.5*(var_c[j+i*ny]+var_c[j+(i+1)*ny]);      
    }
    
}

/**
 *  Interpolates values in x nodes to center nodes
 *  @param var_c[out] values in x-nodes
 *  @param var_x values in c-nodes 
 *  @param x2c_d weight for current postion F_c()=d*F_x() + up*F_x(+)
 *  @param x2c_d weight for x+1 position F_c()=d*F_x() + up*F_x(+)
*/
inline void X2C(double  __restrict_arr  * var_c, double const __restrict_arr*  var_x, double const __restrict_arr* x2c_d,const double __restrict_arr* x2c_up,int nx, int ny)
{
    
    for (int i=0;i<nx;i++)
    {
        const double xd=x2c_d[i];
        const double xup=x2c_up[i];
        for (int j=0;j<ny;j++)
            var_c[j+i*ny]=xd*var_x[j+i*ny] +xup*var_x[j+(i+1)*ny];
    
    }
}


}

#endif