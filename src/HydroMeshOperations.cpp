/*********************************************************************************
 * 
 *                  HydroMeshOperations.cpp (Morfo70) (TO REMOVE IN A FUTURE VERSION)
 * 
 *      Methods for performing matrix operations in order to:
 *      - Interpolate values from  c,x,y nodes to c,x,y nodes. 
 *      - Calculates x, and y derivatives
 *      - Computes inverse and square inverse function
 *
 *   Methods and parameters description in "Hydro.h"  
 *  
 * *******************************************************************************/
#include "Hydro.h"

void Hydro::C2X(double*   var_x, double *  var_c,int nx, int ny)
{
    
    for (int i=0;i<nx-1;i++)
    {
    
    
        for (int j=0;j<ny;j++)    
            var_x[j+i*ny]=0.5*(var_c[j+i*ny]+var_c[j+(i+1)*ny]);    

        
    }    

}

void Hydro::X2C(double* var_c, const double * var_x,const double* x2c_d,const double* x2c_up,int nx, int ny){
    const double * v_x=var_x;
    const double * v_x_up=v_x+ny;
    double* v= var_c;
    for (int i=0;i<nx;i++)
    {
        double xd=x2c_d[i];
        double xup=x2c_up[i];
        for (int j=0;j<ny;j++)
            v[j]=xd*v_x[j] +xup*v_x_up[j];

        v_x=v_x_up;
        v_x_up+=ny;
        v+=ny;
    }
}

void Hydro::C2Y(double* var_y, const double* var_c, int nx, int ny)
{
    const double * v=var_c;    
    double* v_y= var_y;
    for (int i=0;i<nx;i++)
    {
        for (int j=0;j<ny-1;j++)
            v_y[j]=0.5*(v[j+1]+v[j]);
        
        v_y[ny-1]=0.5*(v[0]+v[ny-1]);

        v+=ny;        
        v_y+=ny;
    }
}

void Hydro::C2Y(double* var_y, const double* var_c, int ny)
{
    for (int j=0;j<ny-1;j++)
        var_y[j]=0.5*(var_c[j+1]+var_c[j]);
    
    var_y[ny-1]=0.5*(var_c[0]+var_c[ny-1]);
    
}

void Hydro::Y2C(double* var_c, const double* var_y, int nx, int ny)
{
    const double * v=var_y;    
    double* v_c= var_c;
    for (int i=0;i<nx;i++)
    {
        for (int j=1;j<ny;j++)
            v_c[j]=0.5*(v[j]+v[j-1]);
        
        v_c[0]=0.5*(v[0]+v[ny-1]);

        v+=ny;        
        v_c+=ny;
    }
}

void Hydro::Y2C(double* var_c, const double* var_y, int ny)
{
    
    for (int j=1;j<ny;j++)
        var_c[j]=0.5*(var_y[j]+var_y[j-1]);
    
    var_c[0]=0.5*(var_y[0]+var_y[ny-1]);

    
}

void Hydro::difC2X(double* df_dx, const double * var_c,const double* dx,int nx, int ny)
{
    const double * v=var_c;
    const double * v_up=v+ny;
    double* d_dx= df_dx;
    for (int i=0;i<nx-1;i++)
    {
        double inv_x=1/dx[i];
        for (int j=0;j<ny;j++)
            d_dx[j]=inv_x*(v_up[j]-v[j]);           

        v=v_up;
        v_up+=ny;
        d_dx+=ny;
    }
}

void Hydro::difC2Y(double* _df_dy, const double* var_c, double dy, int nx, int ny){
    const double * v=var_c;    
    double* df_dy= _df_dy;
    double inv_dy=1/dy;
    for (int i=0;i<nx;i++)
    {
        for (int j=0;j<ny-1;j++)
            df_dy[j]=(v[j+1]-v[j])*inv_dy;
        
        df_dy[ny-1]=(v[0]-v[ny-1])*inv_dy;

        v+=ny;        
        df_dy+=ny;
    }
}

void Hydro::difY2C(double* _df_dy, const double* var_y, double dy, int nx, int ny){
    const double * v=var_y;    
    double* df_dy= _df_dy;
    double inv_dy=1/dy;
    for (int i=0;i<nx;i++)
    {
        for (int j=1;j<ny;j++)
            df_dy[j]=(v[j]-v[j-1])*inv_dy;
        
        df_dy[0]=(v[0]-v[ny-1])*inv_dy;

        v+=ny;        
        df_dy+=ny;
    }
}

void Hydro::difY2C(double* df_dy, const double* var_y, double dy,int ny){
        
    double inv_dy=1/dy;
    
        for (int j=1;j<ny;j++)
            df_dy[j]=(var_y[j]-var_y[j-1])*inv_dy;
        
        df_dy[0]=(var_y[0]-var_y[ny-1])*inv_dy;

    
    
}

void Hydro::difC2C_x_addBoundary(double * _df_dx, const double* var_c, const double* difx_d,const double* difx_up,const double* difx_low,int nx, int ny){
    
    double dup=difx_up[0];
    double dlow=difx_low[0];
    double dd=difx_d[0];
    const double* v=var_c;
    const double* v_up=var_c+ny;
    const double* v_low;
    double * df_dx=_df_dx;

    //Primera fila, continiodad v0=2v1-v2
    for (int j=0;j<ny;j++)    
        df_dx[j]=dlow*(2*v[j]-v_up[j])+dd*v[j]+v_up[j]*dup;    
    
    for (int i=1;i<nx-1;i++)
    {        
        v_low=v;
        v=v_up;
        v_up+=ny;
        dup=difx_up[i];
        dd=difx_d[i];
        dlow=difx_low[i];
        df_dx+=ny;
        for (int j=0;j<ny;j++)  
            df_dx[j]=dlow*v_low[j]+dd*v[j]+v_up[j]*dup;    
        
    }

    
    v_low=v;
    v=v_up;    
    dd=difx_d[nx-1];
    dlow=difx_low[nx-1];
    df_dx+=ny;
    //Ultima fila, v_up=0;
    for (int j=0;j<ny;j++)  
            df_dx[j]=dlow*v_low[j]+dd*v[j];    


}

void Hydro::difC2C_y(double* _df_dy, const double * var_c, double dy,int nx, int ny)
{
    double inv_2dy=0.5/dy;
    const double* v=var_c;
    double * df_dy=_df_dy;
    for (int i=0;i<nx;i++)
    {
        for (int j=1;j<ny-1;j++)        
            df_dy[j]=inv_2dy*(v[j+1]-v[j-1]);
        
        df_dy[0]=inv_2dy*(v[1]-v[ny-1]);
        df_dy[ny-1]=inv_2dy*(v[0]-v[ny-2]);
        
        v+=ny;
        df_dy+=ny;



    }
}

/**
 * 1./F
 * */ 
void Hydro::inv(double * inv_F, const double * F, int nx, int ny){
    
    double * inv_f=inv_F;
    const double * f=F;
    for (int i=0;i<nx;i++)
    {
        for (int j=0;j<ny;j++)
        {
            inv_f[j]=1/f[j];
        }
        f+=ny;
        inv_f+=ny;
    }
    

}

/**
 * 1./F^2
 * */
void Hydro::inv2(double * inv2_F, const double * F, int nx, int ny){
    
    double * inv2_f=inv2_F;
    const double * f=F;
    for (int i=0;i<nx;i++)
    {
        for (int j=0;j<ny;j++)
        {
            inv2_f[j]=1/(f[j]*f[j]);
        }
        f+=ny;
        inv2_f+=ny;
    }
    

}