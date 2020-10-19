/*********************************************************************************
 * 
 *                  Grid.h (Morfo70)
 * 
 *  Grid class stores mesh and nodes information. 
 *  Calculate and store coeffcients for performing  x-derivatives and interpolation
 *  Implementation in:
 *      - Grid.cpp 
 * 
 * *******************************************************************************/

#ifndef GRID_H_
#define GRID_H_

#include "Netcdf.h"
using namespace std;

class Grid{

public:
  double* x;
  double* dx; //vector with cells x-size
  double* dx_x;  //vector with x-cells x-size

  //Coefficients for x-derivatives in x-cells
  // dF_dx()=low*F(-) + d*F() + up*F(+) 
  double* difx_x_up; 
  double* difx_x_d;
  double* difx_x_low;
  //Coefficients for second x-derivatives in x-cells
  // d2F_dxx()=low_2*F(-) + d_2*F() + up_2*F(+) 
  double* dif2x_x_up;
  double* dif2x_x_d;
  double* dif2x_x_low;
  //Coeficcients for x-derivatives in centers cells 
  double* difx_c_up;
  double* difx_c_d;
  double* difx_c_low;
  //Coeficcients for second x-derivatives in x-cells
  double* dif2x_c_up;
  double* dif2x_c_d;
  double* dif2x_c_low;  
  //Coeficcients for interpolation X to C
  //F_c()=d*F_x() + up*F_x(+)
  double* x2c_up;
  double* x2c_d;
  
  double* y;
  double dy;
  int nx;
  int ny;
  int size;
  int* indexUpY;  //Only por waves indexes for lateral boundary conditions TO IMPROVE
  int* indexLowY;
public:

/**
 *  Constructor, read x,y coordinates from netcdf file, and build Morfo70 mesh
 * */
Grid(string bathyFile);

/**
 *  Destructor, free memory. Only called when Morfo70 end
 * */
~Grid();

private:

/**
 *  calculate first x-derivatives weights
 *  dF_dx()=low*F(-) + d*F() + up*F(+) 
 *  @param dx0 x-size  of the x-previous cell
 *  @param dx1 x-size  of the x-current cell
 *  @param up[out] up coefficient
 *  @param d[out] center coefficient
 *  @param low[out] low coefficient
 * */
void calculateDifCoefficients(double dx0,double dx1,double & up, double & d, double & low);

/**
 *  calculate second x-derivativew weights
 *  d2F_dxx()=low*F(-) + d*F() + up*F(+) 
 *  @param dx0 x-size  of the x-previous cell
 *  @param dx1 x-size  of the x-current cell
 *  @param up[out] up coefficient
 *  @param d[out] center coefficient
 *  @param low[out] low coefficient
 * */
void calculateDif2Coefficients(double dx0,double dx1,double & up, double & d, double & low);
};


#endif