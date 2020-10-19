/*********************************************************************************
 * 
 *                  Bathymetry.h (Morfo70)
 * 
 *  Bathymetry class stores bathymetry and water Depth data
 *  Methods init vars and load data from Morfo70 resume file
 *  Implementation in:
 *      - Bathymetry.cpp 
 * 
 * *******************************************************************************/

#ifndef BATHYMETRY_H_
#define BATHYMETRY_H_
#include "Grid.h"
#include "OptionsParametersConstants.h"






class Bathymetry {
private:
	//Mesh info
	int nx;
	int ny;
	int size;	
public:
	double* Z;  //Constant Profile bathymetry (Z=0 in MSL Z>0 in water Z<0 in dry region)
	double* Zb; //Beach bathymetry Zb(x,y)=Z(x,y)-h(x,y)
	double* h;  //Bottom  perturbation h(x,y), from alongshore homogeneus bathy (Z)
	double* D;  //Water Depth D(x,y)=Zb(x,y)+ Zs(x,y)->Free Surface	

	/**
	 * Constructor, init vars (Z and h) from bathy ncFile. Water Depth -> D=Zb+zs0
	 * @param bathyFile morfo70 nc file
	 * @param zs0 ini free surface for the all domain
	 * @param grid mesh
	 * */
	Bathymetry(string bathyFile,double zs0, Grid* grid);	

	/**
	 * Destructor, free memory. Only called when morfo70 ends.
	 * */
	~Bathymetry();

	/**
	*  load data from previous computation file (Z, D and h) and set necesary variables in order to
	*  resume simulation
	*  @param resumeFile morfo70 netCDF file with time-step results 
	* */
	void setBathyFromResumeFile(string resumeFile);
};

#endif /* BATHYMETRY_H_ */
