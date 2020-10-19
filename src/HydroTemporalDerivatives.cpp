/*********************************************************************************
 * 
 *                  HydroTemporalDerivatives.cpp (Morfo70)
 * 
 *          DEPRECATED
 * 
 * *******************************************************************************/
#include "Hydro.h"

// void Hydro::getdQx_dt(int column){
//     const int ny=this->grid->ny;
//     const int i=column;

//     double __restrict_arr const  *qx=this->Qx+(i+i)*ny;
//     double __restrict_arr const  *u=this->U+(i+i)*ny;
//     double __restrict_arr const  *v=this->V+(i+i)*ny;
//     double __restrict_arr const  *forcing_x=this->Forcing_x+i*ny;
//     double __restrict_arr const  *d=this->bathy->D+i*ny;
//     double __restrict_arr const  *urms_x=this->Urms_x+i*ny;
// }

// void Hydro::calculateTurbulence_x(){
//     int nx=this->grid->nx;    
//     int ny=this->grid->ny;    
//     double * turb_x=this->Turbulence_x;
//     double * u_low=this->U;
//     double * u=this->U+ny;
//     double * u_up=this->U+2*ny;
//     double * v=this->V+ny;
//     double * v_up=this->V+2*ny;
//     double * dx=this->grid->dx;
//     double * difx_x_up=this->grid->difx_x_up;
//     double * difx_x_d=this->grid->difx_x_d;
//     double * difx_x_low=this->grid->difx_x_low;
//     double * dif2x_x_up=this->grid->dif2x_x_up;
//     double * dif2x_x_d=this->grid->dif2x_x_d;
//     double * dif2x_x_low=this->grid->dif2x_x_low;

//     for (int i=0;i<nx-1;i++)
//     {
//         double inv_dx=1/dx[i];
//         for (int j=1;j<ny;j++)
//         {
//             //double dv_dx=0.5*inv_dx(v[j]+v[j-1]-);
//         }
//     }
// }