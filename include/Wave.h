/*********************************************************************************
 * 
 *                  WaveRow (Morfo70)
 * 
 *  WaveRow class stores wave and rollers data of a domain row. (x=x_row, y=(0,Ly))
 *  Obtains reftacted wave ray and  propagate rollers from the x-upper row
 *  
 * 
 * *******************************************************************************/
#ifndef WAVE_H_
#define WAVE_H_


#include "WaveRow.h"

class Wave{
public:    
    
private:
    int ny;
    int nx;
    int size;
    double* x;
    int* indexUpY;
    int* indexLowY;
    double* dx;         
    double dy;
    double *D;    
    double *U;
    double *V;
    double *U_w;
    double *CosTheta;
    double *SinTheta;
    
public:
    double *K;    //NABIL TMP
    double *Theta;
    double *Cg;
    double *H;
    double *E;
    double *Urms;    
    double *Sxx_w;
    double *Sxy_w;
    double *Syy_w;
    double *Dw;
    double *E_r;
    double *D_r;
    double *Sxx_r;
    double *Sxy_r;
    double *Syy_r;
    double *Sxx;
    double *Sxy;
    double *Syy;
   
    
   
public:

    /**
    *  Constructor, allocate memory and initalizes field and dummie variables.
    *  @param grid morfo70 grid for referencing dimensions and cell data
    *  @param grid morfo70 batimetrix for referencingwater depth
    * */
    Wave(Grid* grid,Bathymetry* bathy);


    /**
     *  Destructor, free memory. Is called by main.cpp when morfo70 ends.
     * */
    ~Wave();

    
    /**
     * Computes wave and rollers propagation, row by row. 
     * @param hOff Offshore wave Height
     * @param thetaOff Offshore wave direction
     * @param thetaOff Wave period
     * */
    void calculateWaveField(double hOff, double thetaOff, double T);

    
    /**
     *  set currents in order to apply wave-current interactions
     *  @param _U Crosshore flow 
     *  @param _V longshore flow 
     * */
    void setUV(double __restrict_arr const  * U,double __restrict_arr const * V);  

    
    /**
     *  set currents from resume file 
     *  @param resumeFile morfo70 netCDF file with time-step results 
     * */
    void setUVFromResumeFile(string resumeFile); 


    /**
     *  set wave parameters from resume file, NO NECESSARY
     *  @param resumeFile morfo70 netCDF file with time-step results 
     * */
     void setWaveEnergyFromResumeFile(string resumeFile); 


private:
    /**
     * Stores in Wave object data from a wave row computation
     * @param xIndex mesh (grid) x-index of the row
     * @param r WaveRow object
     * */
    void storageRow(const int xIndex, const WaveRow* r);
    
    
    /**
     * Get the maximum cell x-witdh of the next row taking acount wave incidence. 
     * if this value is higher than the corresponding mesh cell, return cell x-witdh 
     * 
     * OUT: return true if the row would be equal to de mesh (dx_row=dx_mesh), 
     * false otherwise
     * @param x0 x-position of wave incident row
     * @param xIndex mesh (grid) x-index of the nearest previous  mesh row
     * @param rm1   WaveRow object of the incidente  raw
     * @param dx_row[out]  computed x-width of the row.
     * */
    bool getDx(const double dx,const int xIndex,const WaveRow* rm1,double &dx_row);

};


#endif 