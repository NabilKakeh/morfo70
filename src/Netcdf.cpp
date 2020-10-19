/*********************************************************************************
 * 
 *                  Netcdf.cpp (Morfo70)
 * 
 *  Netcdf wrapper for morfo70 files, implementation
 *
 *  Methods and parameters description in "Netcdf.h"  
 * 
 * *******************************************************************************/


#include "Netcdf.h"
namespace MorfoIO{


Netcdf::Netcdf() {
}

Netcdf::~Netcdf() {
	this->Close();
}

void Netcdf::Open(string file, OpenMode openMode){
	switch (openMode){
	case Netcdf::OpenMode::Read:
		nc_open(file.c_str(),NC_NOWRITE | NC_NETCDF4,&this->ncId);
		this->ReadDims();
		break;
	case Netcdf::OpenMode::Write:
		int succes=nc_create(file.c_str(),NC_WRITE | NC_NETCDF4,&this->ncId);
		if (succes!=0) //Netcdf3 (no compression)		
			succes=nc_create(file.c_str(),NC_WRITE,&this->ncId);					
		break;
	}
}

void Netcdf::Close(){
	if (this->ncId>0)
			nc_close(this->ncId);
}

void Netcdf::ReadDims(){
	int xDimId,yDimId;
	char foo[32] = "";
	size_t xLen,yLen;
	nc_inq_dimid(ncId,"x",&xDimId);
	nc_inq_dim(ncId,xDimId,foo,&xLen);
	nc_inq_dimid(ncId,"y",&yDimId);
	nc_inq_dim(ncId,yDimId,foo,&yLen);
	this->nx=xLen;
	this->ny=yLen;
	this->size=nx*ny;
}

void Netcdf::ReadXY(double* x,double* y){
	nc_inq_varid(ncId,"x",&xId);
	nc_inq_varid(ncId,"y",&yId);
	nc_get_var_double(ncId,xId,x);
	nc_get_var_double(ncId,yId,y);

}

void Netcdf::ReadVar(string varName, double* var){
	int varId;
	//var=new double[this->size];
	double tmpVar[nx][ny];
	nc_inq_varid(ncId,varName.c_str(),&varId);
	nc_get_var_double(ncId,varId,&tmpVar[0][0]);

	int counter=0;

		for (int j = 0; j < nx; ++j) {
			for (int i = 0; i < ny; ++i) {
				counter=i+j*ny;
				var[counter]=tmpVar[j][i];
				}
		}

}

void Netcdf::ReadAtt(string attName, double * att){
	int status=nc_get_att_double(ncId,NC_GLOBAL,attName.c_str(),att);

}

void Netcdf::ReadVarY(string varName, double* var){
	int varId;	
	nc_inq_varid(ncId,varName.c_str(),&varId);	
	nc_get_var_double(ncId,varId,var);
	
		
}

void Netcdf::WriteXY(double* x,double* y,int nx, int ny){
	nc_redef(ncId);
	this->nx=nx;
	this->ny=ny;
	this->size=nx*ny;
	nc_def_dim(ncId,"x",nx,&xdimId);
	nc_def_var(ncId,"x",NC_DOUBLE,1,&xdimId,&this->xId);
	nc_def_dim(ncId,"y",ny,&ydimId);
	nc_def_var(ncId,"y",NC_DOUBLE,1,&ydimId,&this->yId);
	nc_enddef(ncId);
	nc_put_var_double(ncId,this->xId,x);
	nc_put_var_double(ncId,this->yId,y);


}

void Netcdf::WriteVar(string varName, double* var){
		nc_redef(ncId);
		int varId;
		int  dims[2]={xdimId,ydimId};
		nc_def_var(ncId,varName.c_str(),NC_DOUBLE,2,dims,&varId);
		nc_def_var_deflate(ncId,varId,1,1,9);
		nc_enddef(ncId);
		nc_put_var_double(ncId,varId,var);
}

void Netcdf::WriteVarY(string varName, double* var)
{
	nc_redef(ncId);
	int varId;
	//int  dims[1]={ydimId};
	nc_def_var(ncId,varName.c_str(),NC_DOUBLE,1,&ydimId,&varId);
	nc_def_var_deflate(ncId,varId,1,1,9);
	nc_enddef(ncId);
	nc_put_var_double(ncId,varId,var);
}

void Netcdf::WriteAtt(string attName, double  att){
	int succes=nc_put_att_double(ncId,NC_GLOBAL,attName.c_str(),NC_DOUBLE,1,&att);

}

}






