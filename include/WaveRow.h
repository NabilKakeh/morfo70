/*********************************************************************************
 * 
 *                  WaveRow (Morfo70)
 * 
 *  WaveRow class stores wave and rollers data of a domain row. (x=x_row, y=(0,Ly))
 *  Obtains reftacted wave ray and  propagate rollers from the x-upper row.
 *  Implementation in:
 *      - WaveRow.cpp
 *      - WaveDir.cpp
 *      - WaveEnergy.cpp
 *      - WaveParametrizations.cpp
 * 
 * *******************************************************************************/

#ifndef WAVEROW_H_
#define WAVEROW_H_


#include "Bathymetry.h"
#include <cmath>
#include "OpenMP.h"


class WaveRow{
    private:
        int ny;
        double dy;
        double dx;
        double inv_2dy;
        double * inv_S_i;
        double * SinTheta_i;
        double * CosTheta_i;
        double * Theta_i;
        double * Y_i;             
        int* indexUpY;
        int* indexLowY;                         
        double * ddw_dh;
        double * ddf_dh;
        double * G;
        double * Gp;        
        double *WC; //Wave Current Interaction
        double *WC_r; //Wave Rollers Interaction

         
    public:              
        double *D;   
        double *lK;
        double *U;
        double *V;
        double *U_w;
        double *TanTheta;
        double *CosTheta;        
        double *SinTheta;
        double *Theta;
        double *Sigma;
        double *C;
        double *Cg;
        double * dw;        
        double * df;
        double * sxx;
        double * sxy;
        double * syy;
        double * sxx_r;
        double * sxy_r;
        double * syy_r;
        double * lE;
        double * E;
        double * Dw;
        double * H;        
        double * E_r;
        double * D_r;
        double * d_r;        
         double *Cx;
        double *Cy;
         double *Cgx;
        double *Cgy;
        double *Urms_H;
        
    /**
    *  Constructor: Allocate memory and Initalizes all wave and rollers variables
    *  @param ny number of Y nodes (size of the row)
    *  @param dy y-width of row cells
    *  @param indexUpY row positions of the nodes y+1 nodes
    *  @param indexLowY row positions of the nodes y-1 nodes
    * */
    WaveRow(int ny,double dy,  int* indexUpY,int* indexLowY);    

    /**
    *  Destructor, free memory, it is only called by Wave.cpp when propagations ends 
    * */
    ~WaveRow();    

        
    
    /**
     *  Init row values before compute WITHOUT interpolation
     *  OUT: set water depth , currents and x-cell-lengh
     *  @param row_m1 previus row (for initalising H)
     *  @param dx cell x-width 
     *  @param D Water Depth
     *  @param U x-current (for wave current interaction)
     *  @param V y-current (for wave current interaction)
     *  @param U_w Current projected in wave direction -> Uw=U*cos(Theta)+V*sin(Theta) (for wave current interaction)
     * */
    void initRow(const WaveRow* rm1,const double dx, const double* D,const double * U, const double* V,const double* U_w);
    
    
    /**
     *  Init row values before compute WITH interpolation (it is a sub-grid cell)
     *  OUT: set water depth , currents and x-cell-lengh
     *  @param row_m1 previus row (for initalising H)
     *  @param dx sub-cell x-width 
     *  @param dx_grid grid-cell x-width 
     *  @param D_p1 Water Depth of the upper-x grid node 
     *  @param U_p1 x-current of the upper-x grid node 
     *  @param V_p1 y-current of the upper-x grid node 
     *  @param U_w_p1 Current projected in wave direction  of the upper-x grid node 
     *  @param D_m1 Water Depth of the lower-x grid node 
     *  @param U_m1 x-current of the lower-x grid node 
     *  @param V_m1 y-current of the lower-x grid node 
     *  @param U_w_m1 Current projected in wave direction  of the lower-x grid node 
     * */
    
    
    void initRow(const WaveRow* rm1,const double dx, const double dx_grid, const double *D_p1, const double *U_p1, const double *V_p1, const double *U_w_p1, const double *D_m1,const double *U_m1, const double *V_m1, const double *U_w_m1);
     /**
     *  Init first row (x=0) with along-y homogeneous  values for each node.
     *  OUT: set offshore values and calculate related parametrizations  (K,C,Cg,Dw,Df)
     *  @param hOff Ofsshore Wave Height
     *  @param thetaOff Offshore wave Directions
     *  @param T Wave period
     *  @param D Water Depth
     *  @param U x-current (for wave current interaction)
     *  @param V y-current (for wave current interaction)
     *  @param U_w Current projected in wave direction -> Uw=U*cos(Theta)+V*sin(Theta) (for wave current interaction)
     * */
    void initFirstRow(const double hOff, const double thetaOff, const double T, const double *D, const double *U, const double *V, const double *U_w);    
    
    
    /**
    *  Calculate the parallel dispersion relationship with doppler effect
    *  OUT: for a waveRow set log(K), Sigma, Urms/H, C and Cg
    *  @param T wave period  
    *  @param np number of processor (for serial mode np=1) 
    */
    void dispersion_parallel(const double T,const int* _indexIni, const int* _chunk,int np);

    
    /**
    *  Computes in a parallel way solveDir method for a row
    *  @param np number of processor (for serial mode np=1) 
    * */
    void solveDir_parallel(const WaveRow* rm1,const int* _indexIni, const int* _chunk,int np);   
    
    
    /**
    *  Computes in a parallel way solveEnergy for a row
    *  @param np number of processor (for serial mode np=1) 
    * */
    void solveEnergy_parallel(const int* _indexIni, const int* _chunk,int np);
    
    /**
    *  Set the constant terms (no change in NR iterations)  for the wave energy objective function and its derivative
    *  OUT: sets G0 (obcjetive function) and Gp0 (derivative of the objective functio)
    *  @param rm1 previus row (for compute derivatives and incident ray)
    * */
    void energyPreprocessor(const WaveRow *rm1);

    /**
    *  Solves implicitly the WEEBEQ for a piece of the row, using Newton-Raphson
    *  OUT: sets lEr->log(Er) and D_r/E (Dr->Rollers disipatio)
    *  @param rm1 previus row (for compute derivatives and incident ray) 
    * */    
    void solveRollers(const WaveRow* rm1);    
   
private:
    
    /**
    *  Calculate dispersion relationship with doppler effect, using Newton Raphson method 
    *  OUT: for a waveRow set the wave number-> log(K),   and  Sigma, Urms/H, C and Cg
    *  @param w angular frequency 2*pi/T
    */
    void dispersion(int indexIni,int chunk,double w);    
    
    /**
    *  Calculate dispersion relationship with doppler effect, using Ondina 2018 aproximation 
    *  OUT: for a waveRow set the wave number-> log(K),   and  Sigma, Urms/H, C and Cg
    *  @param w angular frequency 2*pi/T   
    */
    void dispersion_fast(int indexIni,int chunk,double w); 

    /**
    *  Computes the wave ray refraccion calculating rot(K)=0 in a ray form: dTheta_ds=d(logK)_dst : s-> incident ray-dir, st-> transverse ray-dir
    *  OUT: For a chunk of a wave row, sets Theta, sin, cos and tan. Wace and rollers Radiaton stresses Sij/E, and projected celerities
    *  @param _tanTheta_m1 tangent of the previous row
    *  @param _lK_m1  log of the wave number -> log(K), for the previous row
    *  @param _lK  log of the wave number
    * */   
    void solveDir(const double* tanTheta_m1,const double* lK_m1,const double* lK,int indexIni,int chunk);      

    /**
    *  Solves implicitly the WEBEQ for a piece of the row, using Neton Raphson
    *  OUT: sets lE->log(E), H, Dw, Dw/E and Df/E. 
    * */  
    void solveEnergy(int indexIni, int chunk);
    
    

/**
 * 
 *  Wave Dissipation Parametrizations
 *
 * */
private:
    /**
     *  Virtual method to point wave breaking dissipation parametrizations
     * */
    void (WaveRow::*waveBreakingDissipation)(int, const double*,const double * ,const double * ,const int*,int,double*,double*)=nullptr;                
    

    /**
    *  Calculate the Wave Breaking Dissipation according Church And Thornton 1996
    *  OUT: Dw/E and d(Dw/E)_dH (derivative)
    *  @param _H wave height
    *  @param _Sigma wave angular frequency
    *  @param _D water depth  
    */
    void waveBreakingDissipationCT( int chunk,const double* _H,const double * _Sigma, const double *_D,const int* _nodesToSolve,int nNodesToSolve,double * _dw, double * _ddw_dh);
    
    /**
     *  Calculate the Wave Breaking Dissipation according Thornton and Guza 1985
     *  OUT: Dw/E and d(Dw/E)_dH (derivative)
     *  @param _H wave height
     *  @param _Sigma wave angular frequency
     *  @param _D water depth  
     */
    void waveBreakingDissipationTG( int chunk,const double* _H,const double * _Sigma, const double *_D,const int* _nodesToSolve,int nNodesToSolve,double * _dw, double * _ddw_dh);
    
    /**
     *  Calculate the Wave Breaking Dissipation  for regular waves (Laboratory)
     *  OUT: Dw/E and d(Dw/E)_dH (derivative)
     *  @param _H wave height
     *  @param _Sigma wave angular frequency
     *  @param _D water depth  
     */
    void waveBreakingDissipationRegular( int chunk,const double* _H,const double * _Sigma, const double *_D,const int* _nodesToSolve,int nNodesToSolve,double * _dw, double * _ddw_dh);
    
    /**
     *  Virtual method to point bottom Friction coeffienction parametrizations
     * */
    void (WaveRow::*BottomFriction)(int, const double*,const double  ,const double * ,const int*,int,double*,double*)=nullptr;                
    
    /**
     *  Calculate the Wave Bottom Friction Dissipation  (Df) for constant wave drag friction coefficient (fw)
     *  OUT: Df/E and d(Df/E)_dH (derivative)
     *  @param _H wave height
     *  @param _sigma wave height
     *  @param _Urms_H Urbital velocity at bottom (z0) / H
     */
    void bottomFrictionConstant(int chunk, const double *_H,const double  _sigma,const double * _Urms_H,const int* _nodesToSolve,int nNodesToSolve,double * _df, double * _ddf_dh);

    /**
     *  Calculate the Wave Bottom Friction Dissipation  for wave drag friction coeffienct calculated with Soulsby (1997) formula
     *  OUT: Df/E and d(Df/E)_dH (derivative)
     *  @param _H wave height
     *  @param _sigma wave height
     *  @param _Urms_H Urbital velocity at bottom (z0) / H
     */ 
    void bottomFrictionSoulsby97(int chunk, const double *_H,const double  _sigma,const double * _Urms_H,const int* _nodesToSolve,int nNodesToSolve,double * _df, double * _ddf_dh);

    
    
    
    

};

#endif