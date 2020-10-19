/*********************************************************************************
 * 
 *                  Hydro (Morfo70)
 * 
 *  Hydro class stores currents, and free surface data for  the entire grid [nx*ny]
 *  Computes the Non Linear Shallow Water Equations forced by wave and rollers radiation
 *  stresses
 *  Implementation in:
 *      - Hydro.cpp
 *      - HydroQx.cpp
 *      - HydroQy.cpp
 *      - HydroParametrizations.cpp
 *      - HydroMeshOperations.cpp (TO REMOVE)
 *      - HydroBoundaryConditions.cpp 
 * 
 * *******************************************************************************/
#ifndef HYDRO_H_
#define HYDRO_H_
#include "Wave.h"
#include <cmath>
#include "Math.h"
#include "MeshOperations.h"

class Hydro{
    public:
        double * U;     //Crosshore current in x-node
        double * U_c;   //Crosshore current in center-node
        double * V;     //Lonshore current in y-node
        double * V_c;   //Lonshore current in center-node
        double * Zs;    //Free surface elevation in c-node    
        double * Qx;    //water volume flux (Qx=U*D) in x-node
        double * Qy;    //water volume flux (Qy=U*D) in y-node
        double * dQx_dx;   //water volume flux x-derivative
        double * dQy_dy;   //water volume flux y-derivative
        

    private:    
        
        Wave* wave;
        Bathymetry* bathy;
        Grid* grid;

        //Dummie node-variables for computation
        double * cD;   //Drag
        double * nuD;   //nu*D
        double * Forcing_x;  //Wave and rollers forcing  (dSij_di+...)
        double * Forcing_y;
        double * Urms_x;    // Urms
        double * Urms_y;                
        double * inv_Urms2_x;                
        double * inv_Urms2_y;                
        double * cD_x;  
        double * nuD_x;
        double * dnuD_dx;
        double * dnuDx_dy;
        double * cD_y;
        double * nuD_y;
        double * dnuD_dy;
        double * dnuDy_dx;

        //Attributes specfic for each hydro iteration
        double * guess_Zs;   //Predictor Zs for middle computation step
        double * guess_D;    //Predictor Zs for middle computation step    
        double * Di;      //i=x,y -> water Depth in x, y nodes
        double * inv_Di;  //i=x,y -> 1/Di        
        double * Ui_j;    //i,j=x,y  ->Currenter componente in other postion   i.e. Ux_y
        double * DragTerm; //drag term
        
        bool * IsWater;  //wet region mask

        /**
         *  Solver dummy variables:
         *  A*x=f-> A=[low,d,up] tridiagonal matrix
         * */
        double* v_tmp;
        double * low;
        double * d;   
        double * up;
        double * f;


        //Variables for turbulence 
        double * Turbulence_x;
        double * Turbulence_y;


        
        ///Boundary conditions        
        double bcOffshoreQx_f1;    //Qx0=Qx1*bcOffshoreQx_f1+bcOffshoreQx_f0
        double bcOffshoreQx_f0;
        double bcOffshoreQy;    //Qy0=Qy1*bcOffshoreQy
        

    public:
        
        
        /**
        *  Constructor: Allocate memory and Initalizes all variables
        *  @param grid pointer to Morfo70 grid
        *  @param bathy pointer to Bathymetry (water Depth)
        *  @param wave poiner to Wave Field (forcing)
        * */
        Hydro(Grid* grid,Bathymetry* bathy,Wave* wave);
        
        
        /**
        *  Destructor, free memory, it is only called when Morfo70 ends
        * */
        ~Hydro();


        /**
        *  load data from previous computation file and set necesary variables in order to
        *  resume simulation
        *  @param resumeFile morfo70 netCDF file with time-step results 
        * */
        void setHydroFromResumeFile(string resumeFile);

        /**
         * Get the wave forcing and solves currents and fres surface for  n time steps
         * (n=nIter). The time step (dt) is set by Parameters::Hydro::dt
         * @param nIter , number of iterations (time steps) to run under the same wave conditions
         * */
        void Run(int nIter);
    
    
    private:

        /**
         * Solves in two steps ( 1. x-integration and 2. y-integration) the x-momentum NLSW equation
         * */
        void solveQx();

        /**
         * Solves in two steps ( 1. y-integration and 2. x-integration) the y-momentum NLSW equation
         * */
        void solveQy();
        
        /**
         *  Set the wave and rollers forcing term of the  x-momentum NLSW equation. (dSij_di)
         * */
        void getForcing_x();        
        
        /**
         *  Set the wave and rollers forcing term of the  y-momentum NLSW equation. (dSij_di)
         * */
        void getForcing_y();        
        
        
        
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    Interpolate Water Depth in x-y cells 
        ////
        /////////////////////////////////////////////////////////////////////////////////////
        /**
         * Get the water detph ant its inverse in the offshore boundary position x=-dx/2, by linear extrapolation
         * */
        void getD_x0();
        
        /**
         * Interpolates the water detph from center to x-cells (and computes 1/D)
         * @param xIndex index of the grid x position of the c-nodes
         * */
        void getD_x(int xIndex);
        
        /**
         * Get the water detph ant its inverse in the offshore boundary position x=-dx, by linear extrapolation
         * */
        void getD_y0();        

        /**
         * Interpolates the water detph from center to y-cells (and computes 1/D)
         * @param xIndex index of the grid x position of the c-nodes
         * */
        void getD_y(int xIndex);
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    Interpolate Water Depth Methodos
        ////
        /////////////////////////////////////////////////////////////////////////////////////
                     

        
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    Solver Methods Qx predictor        
        ////
        /////////////////////////////////////////////////////////////////////////////////////  
        /**
         * set the wet and dry x-nodes of as x position
         * @param xIndex index of the grid x position of the c-nodes
         * */
        void setWaterNodes_x(int xIndex);

        /**
         * Set f terms of tridiagonal system of the prediction step A*Qxx=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */
        void setSolverQxx_f(int xIndex);
        
        /**
         * Set diagonal terms of tridiagonal system of the prediction step A*Qxx=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */
        void setSolverQxx_d(int xIndex);        
        
        /**
         * Set upper terms of tridiagonal system of the prediction step A*Qxx=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */
        void setSolverQxx_up(int xIndex);        
        
        /**
         * Set xIndex=0 lower terms of tridiagonal system of the prediction step A*Qxx=f, A=[up d low]        
         * */
        void setSolverQxx0_low();   
        
        /**
         * Set lower terms of tridiagonal system of the prediction step A*Qxx=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */
        void setSolverQxx_low(int xIndex);   

        /**
         * Set Solver  Offshore Boundary conditions in Qx prediction step
         * */
        void setSolverQxx_BC();         
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    END Solver Methods Qx predictor        
        ////
        /////////////////////////////////////////////////////////////////////////////////////

        
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    Solver Methods Qx projector
        ////
        /////////////////////////////////////////////////////////////////////////////////////  
        /**
         * Set f terms of tridiagonal system of the projector step A*Qxy=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */      
        void setSolverQxy_f(int xIndex);
        
        /**
         * Set diagonal terms of tridiagonal system of the projector step A*Qxy=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */      
        void setSolverQxy_d(int xIndex);
        
        /**
         * Set upper diagonal terms of tridiagonal system of the projector step A*Qxy=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */      
        void setSolverQxy_up(int xIndex);
        
        /**
         * Set lower diagonal terms of tridiagonal system of the projector step A*Qxy=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */      
        void setSolverQxy_low(int xIndex);
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    END Solver Methods Qx projector
        ////
        /////////////////////////////////////////////////////////////////////////////////////        
       

        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    Solver Methods Qy predictor
        ////
        /////////////////////////////////////////////////////////////////////////////////////       
        /**
         * set the wet and dry y-nodes of a x position
         * @param xIndex index of the grid x position of the c-nodes
         * */
        void setWaterNodes_y(int Xindex);

        /**
         * Set f terms of tridiagonal system of the Qy prediction step A*Qyy=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */       
        void setSolverQyy_f(int xIndex);

        /**
         * Set diagonal terms of tridiagonal system of the Qy prediction step A*Qyy=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */       
        void setSolverQyy_d(int xIndex);

        /**
         * Set upper diagonal terms of tridiagonal system of the Qy prediction step A*Qyy=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */       
        void setSolverQyy_up(int xIndex);

        /**
         * Set lower diagonal terms of tridiagonal system of the Qy prediction step A*Qyy=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */       
        void setSolverQyy_low(int xIndex); 
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    END Solver Methods Qy predictor
        ////
        /////////////////////////////////////////////////////////////////////////////////////
        

        
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    Solver Methods Qy projector
        ////
        /////////////////////////////////////////////////////////////////////////////////////
        /**
         * Set f terms of tridiagonal system of the Qy projector step A*Qyx=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */       
        void setSolverQyx_f(int xIndex);       
        
        /**
         * Set diagonal terms of tridiagonal system of the Qy projector step A*Qyx=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */       
        void setSolverQyx_d(int xIndex);

        /**
         * Set upper diagonal terms of tridiagonal system of the Qy projector step A*Qyx=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */       
        void setSolverQyx_up(int xIndex);

        /**
         * Set lower diagonal terms of tridiagonal system of the Qy projector step A*Qyx=f, A=[up d low]
         * @param xIndex index of the grid x position 
         * */       
        void setSolverQyx_low(int xIndex);    
        
        /**
         * Set xIndex=0 lower terms of tridiagonal system of the Qy projector step A*Qyy=f, A=[up d low]        
         * */ 
        void setSolverQyx0_low();     

        /**
         * Set Solver  Offshore and Onshore Boundary conditions in Qy procector step 
         * */
        void setSolverQyx_BC();   
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    END Solver Methods Qy projector
        ////
        /////////////////////////////////////////////////////////////////////////////////////
            
        

        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    After Solver Methods: Get Zs, U,V and offshore and onshore values
        ////
        /////////////////////////////////////////////////////////////////////////////////////
        /**
         * Computes Qx in offhore position after Qx projector step
         * */
        void setQx_BC();
        
        /**
         * Computes Qy in Offhore and Onshore position after Qy projector step
         * */
        void setQy_BC();          

        /**
         * Obtain free surfe guess solving explictily mass conservation equation before Qy step
         * @param xIndex index of the grid x position 
         * */
        void solveZs_guess(int xIndex);
        
        /**
         * Obtain free surfe solving explictily mass conservation equation after Qy step
         * @param xIndex index of the grid x position 
         * */
        void solveZs(int xIndex);     

        /**
         * Get the flow crosshore velocity (U=Qx/D) in the boundary offshore position x=-dx/2         
         * */
        void getU0();             
        
        /**
         * Get the flow crosshore velocity (U=Qx/D) 
         * @param xIndex index of the grid x position 
         * */
        void getU(int xIndex);                     
        
        /**
         * Get the flow longshore velocity (U=Qx/D) in the boundary offshore position x=-dx/2         
         * */
        void getV0();                      
        
        /**
         * Get the flow longshore velocity (U=Qx/D)
         * @param xIndex index of the grid x position 
         * */
        void getV(int xIndex);                     
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    END After Solver Methods: Get Zs, U,V and offshore and onshore values
        ////
        /////////////////////////////////////////////////////////////////////////////////////


        //Metodos para las derivadas temporales
        // void calculateTurbulence_x();
        // void getdQx_dt(int xIndex);

        
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    Mesh Operations, interpolation to c,x and y nodes, and derivatives
        ////    TO REMOVE
        ////
        /////////////////////////////////////////////////////////////////////////////////////                
        void C2X(double* var_x, double *   var_c,int nx, int ny);
        void X2C(double* var_c, const double * var_x,const double* x2c_d,const double* x2c_up,int nx, int ny);
        void C2Y(double* var_y, const double * var_c,int nx, int ny);
        void C2Y(double* var_y, const double * var_c, int ny);
        void Y2C(double* var_c, const double * var_y,int nx, int ny);
        void Y2C(double* var_c, const double * var_y, int ny);
        void difC2X(double* df_dx, const double * var_c,const double* dx,int nx, int ny); 
        void difC2Y(double* _df_dy, const double* var_c, double dy, int nx, int ny);             
        void difY2C(double* _df_dy, const double* var_y, double dy, int nx, int ny);     
        void difY2C(double* _df_dy, const double* var_y, double dy, int ny);     
        
        /** 
         * x derivatives betwen center nodes. For Offshore use continuity -> F0=2*F(1)-F(2).
         * For Onshore the source function is set to 0 -> F(Lx)=0
         * */
        void difC2C_x_addBoundary(double * _df_dx, const double* var_c, const double* difx_c_d,const double* difx_c_up,const double* difx_c_low,int nx, int ny);        
        void difC2C_y(double* df_dy, const double * var_c,double dy,int nx, int ny);     
        void inv(double * inv_f, const double * f, int nx, int ny);
        /** 
         * Performs 1/F^2
         * */
        void inv2(double * inv2_f, const double * f, int nx, int ny);
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    END Mesh Operations
        ////
        /////////////////////////////////////////////////////////////////////////////////////
        

        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    Parametrizations
        ////
        /////////////////////////////////////////////////////////////////////////////////////                
        /**
         * Calculate the general drag term -> tau / (rho*D*Q)
         * */        
        void calculateDragTerm(double  __restrict_arr * drag, double  __restrict_arr * cd,double  __restrict_arr * urms, double  __restrict_arr * u, double  __restrict_arr * v, double  __restrict_arr * inv_d,int ny);
        void calculateDragAndViscosity();
        void (Hydro::*dragCoefficient)( double*, int,double *)=NULL;
        void (Hydro::*viscosity)( double*, double*, double*, double*, int, double *)=NULL;        
        void dragCoefficient_Constant( double* D, int n,double* cD);
        void dragCoefficient_LogDepth( double* D,int n,double* cD);
        void dragCoefficient_ManningStrickler( double* D, int n,double* cD);
        void viscosity_Constant( double*D,double*Dw, double*Dr, double*H, int n, double * nuD );
        void viscosity_BreakingRollers( double*D, double*Dw, double*Dr, double*H, int n, double * nuD );
        void viscosity_BreakingDepth( double*D, double*Dw, double*Dr, double*H, int n, double * nuD );
        void viscosity_BreakingHrms( double*D, double*Dw, double*Dr, double*H, int n, double * nuD );
        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////    END Parametrizations
        ////
        /////////////////////////////////////////////////////////////////////////////////////
};

#endif


