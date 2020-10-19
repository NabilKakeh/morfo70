/*********************************************************************************
 * 
 *                  Netcdf.h (Morfo70)
 * 
 *  Netcdf Class, wrapper to netcdf_c function to read and write bathymetry and 
 *  ouptut Morfo70 files
 *  Implementation in:
 *      - Netcdf.cpp 
 * 
 * *******************************************************************************/

#ifndef NETCDF_H_
#define NETCDF_H_

#include <netcdf.h>
#include <iostream>
using namespace std;


namespace MorfoIO{
class Netcdf {

public:
	
	enum class OpenMode{
		Read,
		Write
	};

public:
	int nx;
	int ny;
	int size; //nx*ny

private:
	int ncId;  
	int xId;
	int yId;
	int xdimId;
	int ydimId;

public:
	
	Netcdf();
	~Netcdf();

	/**
	 *  Open file, in read or write mode 
	 *  TODO: check netcdf3 and netcdf4 compatiblity; check incorrect files and error manager
	 *  @param file morfo70 file
	 *  @param mode read or write
	 * */
	void Open(string file, OpenMode mode);
	/**
	 * Close file, call netcdf_c free function
	* */
	void Close();

	/**
	 * Get nx, and ny from bathy file
	 * */
	void ReadDims();

	/**
	*  read  x and y mesh position from bathy file
	* @param x[out]  x(nx) positions 
	* @param y[out]  y(ny) positions 
	* */
	void ReadXY(double* x,double* y);

	/**
	*  read morfo70 grid var (double)
	* @param varNanme var name
	* @param var[out]  double(nx*ny) vector for filling values
	* */
	void ReadVar(string varName, double* var);

	/**
	* read morfo70 x=cte var (double)
	* @param varNanme var name
	* @param var[out]  double(ny) vector  
	* */
	void ReadVarY(string varName, double* var);
	
	/**
	* read morfo70 global attribute (time?)
	* TODO: check netcdf3 and netcdf4 errors
	* @param attNanme attribute name
	* @param  att[out]  attribute value
	* */
	void ReadAtt(string attName, double * att);

	/**
	*  write x and y mesh values 
	* TODO: DEPRECATED
	* @param x  x(nx) positions 
	* @param y  y(ny) positions 
	* */
	void WriteXY(double* x,double* y,int nx,int ny);
	
	/**
	*  write morfo70 grid var	
	* @param varName  var name
	* @param var  var(nx,ny) -> double(nx*ny) var values
	* */
	void WriteVar(string varName, double* var);

	/**
	*  write morfo70 x=cte  var	
	* @param varName  var name
	* @param var  var(ny) -> double(ny)  values
	* */
	void WriteVarY(string varName, double* var);
	
	/**
	* write morfo70 global attribute 
	* TODO: check netcdf3 and netcdf4 errors
	* @param attNanme attribute name
	* @param  att  attribute value
	* */
	void WriteAtt(string attName, double  att);
};
}

#endif /* NETCDF_H_ */
