/*********************************************************************************
 * 
 *               HydroBoundaryConditions.cpp (Morfo70) 
 * 
 *      Methods to apply Onshore and Offshore boundary conditions in linear System
 *      and post solver
 * 
 *   Methods and parameters description in "Hydro.h"  
 * 
 * *******************************************************************************/

#include "Hydro.h"
void Hydro::setQx_BC(){    
    const int ny=this->grid->ny;    
    const double f1=this->bcOffshoreQx_f1;
    const double f0=this->bcOffshoreQx_f0;    
    double __restrict_arr  * qx1=this->Qx+ny;    
    double __restrict_arr * qx0=this->Qx;    
    for (int j=0;j<ny;j++)    
        qx0[j]=f1*qx1[j]+f0;
    
    
}

void Hydro::setSolverQxx_BC(){    
    int ny=this->grid->ny;
    double f1=this->bcOffshoreQx_f1;
    double f0=this->bcOffshoreQx_f0;
    double __restrict_arr  * low=this->low;
    double __restrict_arr * d=this->d;    
    double __restrict_arr * f=this->f;    
    for (int j=0;j<ny;j++)    
    {
        d[j]+=f1*low[j];
        f[j]-=f0*low[j];
    }
       
    
}

void Hydro::setSolverQyx_BC(){
    const int size=this->grid->size;
    const int ny=this->grid->ny;

    //Offshore
    const double bc=this->bcOffshoreQy;    
    double  __restrict_arr * d=this->d;
    double __restrict_arr   * low=this->low;
    
    for (int j=0;j<ny;j++)    
        d[j]+=bc*low[j];       

    //Onshore
    const double inv_3=1.0/3.0;
    d+=size-2*ny;
    double __restrict_arr * up = this->up+size-2*ny;
    for (int j=0;j<ny;j++)    
        d[j]+=up[j]*inv_3;
    
}

void Hydro::setQy_BC(){ 
    const int size=this->grid->size;
    const int ny=this->grid->ny;        
    
    const double bc=this->bcOffshoreQy;    
    double  __restrict_arr  * qy1=this->Qy+ny;
    double __restrict_arr * qy0=this->Qy;
    for (int j=0;j<ny;j++)
        qy0[j]=bc*qy1[j];
    
    //Onshore
    const double inv_3=1.0/3.0;
    qy0=this->Qy+size;
    qy1=qy0-ny;
    for (int j=0;j<ny;j++)
        qy0[j]=qy1[j]*inv_3;
    
}








