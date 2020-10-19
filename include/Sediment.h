/*********************************************************************************
 * 
 *                  Sediment.h (Morfo70)
 * 
 *  Sediment class computes and stores sedimet fluxes. 
 *  Also solves sand volume conservation equation (dZb/dt = 1/(1-p) * div(q) ) 
 *  in order to calculate bottom evolution
 *
 *  Implementation in:
 *      - Sediment.cpp 
 * 
 * *******************************************************************************/
#ifndef SEDIMENT_H_
#define SEDIMENT_H_
#include "Hydro.h"

class Sediment{
    public:
       
        double * Alpha;  //Stirring 
        double * gamma;  //Downslope
        double * Gamma;  //Swash
        double * qx;  //crosshore transport
        double * qy;  //longshore transport
        double * dh_dt;  //bed evolution
    
    private:
        Hydro* hydro;
        Wave* wave;
        Bathymetry* bathy;
        Grid* grid;
        double* maxDw;
        double* Zs_sh;

    public:
       /**
        * Initialize variables
        * */ 
       Sediment(Hydro* hydro,Wave* wave,Bathymetry* bathy,Grid* grid);

       /**
        * Destructor, only called at the end of morfo70
        * */ 
       ~Sediment();

       /**
        *  Compute sedimetn fluxes (qx and qy) using CWS parametrization
        * */
       void calculateSedimentFluxes();

       /**
        *  Update bottom (bathy->h), solving mass conservation equation
        * */
       void updateBottom();

       /**
        *  set fluxes from resume file
        *  @param resumeFile, morfo 70 resume file from previous computation
        * */
       void setSedimentFluxesFromResumeFile(string resumeFile);

    private:       
        /**
         *  Calculate crosshore sediment flux -> qx=alpha*(u-gamma*dh/dx) -Gamma*dh_dx 
         *  for a x postion (colummn)
         *  @param column index of the mesh x-position
         * */
        void getqx(const int column);
       
        /**
        *  Calculate crosshore sediment flux -> qy=alpha*(v-gamma*dh/dy) -Gamma*dh_dy 
        *  for a x postion (colummn)
        *  @param column index of the mesh x-position
        * */
        void getqy(const int column);
        
        /**
        *  Update bottom (bathy->h), solving mass conservation equation for a x position (column)
        *  dh/dt=- 1/(1-p) * (dqx/dx + dqy/dy)
        *  @param column index of the mesh x-position
        * */
        void updateBottomColumn(const int column);

       
        /**
        *  Get Wave Stirring (Alpha) from Constant Wave Stirring parametrization       
        *  @param column index of the mesh x-position
        * */
        void getWaveStirringCWS(const int column);
        
        /**
        *  Get downslope coefficient (gamma) from Constant Wave Stirring parametrization       
        *  @param column index of the mesh x-position
        * */
        void getGammaDownslopeCWS(const int column);
        
        /**
        *  Get Swash transport coefficient
        *  @param column index of the mesh x-position
        * */
        void getGammaSwash(const int column);
       
        /*
         *  Get maximum breaking wave dissipaion along the coast (y) 
         * */
        void getMaxWaveBreakingDissipation();       
        
        /*
         *  Get the free surface at the shoreline
         * */
        void getShorelineZs();       
       

       
};
#endif