/*********************************************************************************
 * 
 *                  Bathymetry.cpp (Morfo70)
 * 
 *  Bathy constructor , destructor and resume implementation
 *
 *  Methods and parameters description in "Bathymetry.h"  
 * 
 * *******************************************************************************/

#include "Bathymetry.h"

Bathymetry::Bathymetry(string bathyFile,double zs0, Grid* grid) {


	//Initialize	
	this->Z=new double[grid->size];
	this->Zb=new double[grid->size];
	this->h=new double[grid->size];
	this->D=new double[grid->size];
	this->nx=grid->nx;
	this->ny=grid->ny;
	this->size=grid->size;	

	 MorfoIO::Netcdf ncBathy=MorfoIO::Netcdf();
	 ncBathy.Open(bathyFile,MorfoIO::Netcdf::OpenMode::Read);
	 ncBathy.ReadVar("Z",this->Z);	 
	 ncBathy.ReadVar("h",this->h);
	 ncBathy.Close();


	//Set Zb and D. D is forced to be largerd than Dmin
	for (int i=0;i<grid->size;i++){		
		Zb[i]=Z[i]-h[i];				
		double d=Zb[i]+zs0;
		if (d<Parameters::Hydro::Dmin) 
			d=Parameters::Hydro::Dmin;		
		D[i]=d;
	}

}

Bathymetry::~Bathymetry() {
	delete[] this->D;
	delete[] this->Z;
	delete[] this->Zb;
	delete[] this->h;
}


void Bathymetry::setBathyFromResumeFile(string resumeFile){
	 MorfoIO::Netcdf ncBathy=MorfoIO::Netcdf();
	 ncBathy.Open(resumeFile,MorfoIO::Netcdf::OpenMode::Read);
	 ncBathy.ReadVar("Z",this->Z);
	 //ncBathy.ReadVar("Zb",this->Zb);
	 ncBathy.ReadVar("h",this->h);
	 //ncBathy.ReadVar("Zs",this->D);	 
	 ncBathy.ReadVar("D",this->D);	 
	 ncBathy.Close();

	 for (int i=0;i<this->size;i++){
		Zb[i]=Z[i]-h[i];		
		//D[i]=D[i]+Zb[i];		
	}
}

