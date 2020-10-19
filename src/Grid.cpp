/*********************************************************************************
 * 
 *                  Grid.cpp (Morfo70)
 * 
 *  Mesh methods implementation 
 *
 *  Methods and parameters description in "Grid.h"  
 * 
 * *******************************************************************************/

#include "Grid.h"




Grid::Grid(string bathyFile){

    //Initialice Dimensions
     MorfoIO::Netcdf ncBathy=MorfoIO::Netcdf();
	 ncBathy.Open(bathyFile,MorfoIO::Netcdf::OpenMode::Read);
	 ncBathy.ReadDims();
	 this->nx=ncBathy.nx;
	 this->ny=ncBathy.ny;
	 this->x=new double[nx];
	 this->y=new double[ny];
	 ncBathy.ReadXY(this->x,this->y);	
     this->dy=y[1]-y[0];
     this->size=this->nx*this->ny;
     ncBathy.Close();    

    //dx,dy
	this->dx=new double[nx-1]; //nx-1, no tiene espaciado al final
	for(int i=0;i<this->nx-1;i++)
		this->dx[i]=x[i+1]-x[i];	

	this->dx_x=new double[nx]; //nx, tiene espaciado con valor el valor offshore
    this->dx_x[0]=dx[0];
	for(int i=1;i<this->nx-1;i++)
			this->dx_x[i]=(this->dx[i]+this->dx[i-1])/2;
	this->dx_x[nx-1]=this->dx[nx-2];
    
    //dif_x
    this->difx_x_up=new double[nx-1]();
    this->difx_x_d=new double[nx-1]();
    this->difx_x_low=new double[nx-1]();
    this->dif2x_x_up=new double[nx-1]();
    this->dif2x_x_d=new double[nx-1]();
    this->dif2x_x_low=new double[nx-1]();
    for (int i=0;i<nx-1;i++)
    {
        this->calculateDifCoefficients(dx_x[i],dx_x[i+1],difx_x_up[i],difx_x_d[i],difx_x_low[i]);
        this->calculateDif2Coefficients(dx_x[i],dx_x[i+1],dif2x_x_up[i],dif2x_x_d[i],dif2x_x_low[i]);
    }

    //dif_x en el centro, [nx] pueden derivar todo el dominio
    this->difx_c_up=new double[nx]();
    this->difx_c_d=new double[nx]();
    this->difx_c_low=new double[nx]();
    this->dif2x_c_up=new double[nx]();
    this->dif2x_c_d=new double[nx]();
    this->dif2x_c_low=new double[nx]();
    
    //Primera columna
     this->calculateDifCoefficients(dx[0],dx[1],difx_c_up[0],difx_c_d[0],difx_c_low[0]);
     this->calculateDif2Coefficients(dx[0],dx[1],dif2x_c_up[0],dif2x_c_d[0],dif2x_c_low[0]);
    for (int i=1;i<nx-1;i++)
    {
        this->calculateDifCoefficients(dx[i-1],dx[i],difx_c_up[i],difx_c_d[i],difx_c_low[i]);
        this->calculateDif2Coefficients(dx[i-1],dx[i],dif2x_c_up[i],dif2x_c_d[i],dif2x_c_low[i]);
    }
    //ultima columna, se crea un nodo central en la mimas posicon que el staggered, dx(end)=dx_x(end)/2=dx(end-1)/2
    this->calculateDifCoefficients(dx[nx-2],dx[nx-2]/2,difx_c_up[nx-1],difx_c_d[nx-1],difx_c_low[nx-1]);
    this->calculateDif2Coefficients(dx[nx-2],dx[nx-2]/2,dif2x_c_up[nx-1],dif2x_c_d[nx-1],dif2x_c_low[nx-1]);

    //X2C
    this->x2c_up=new double[nx];
    this->x2c_d=new double[nx];
    double dx0=dx[0];
    double dx1=dx0;
    double inv_dx1pdx0=1/(dx1+dx0);
    x2c_up[0]=dx0*inv_dx1pdx0;
    x2c_d[0]=dx1*inv_dx1pdx0;
    for (int i=1;i<nx-1;i++)
    {
        dx0=dx1;
        dx1=dx[i];
        inv_dx1pdx0=1/(dx1+dx0);
        x2c_up[i]=dx0*inv_dx1pdx0;
        x2c_d[i]=dx1*inv_dx1pdx0;        
    }
    x2c_up[nx-1]=0.5;
    x2c_d[nx-1]=0.5; 
    
    
    
  

    //Lateral Indexes  TO IMPROVE
    this->indexUpY=new int[ny];
    this->indexLowY=new int[ny];
    for (int i=1;i<ny-1;i++)
    {
        indexUpY[i]=i+1;
        indexLowY[i]=i-1;
    }
    indexUpY[0]=1;
    indexUpY[ny-1]=0;
    indexLowY[0]=ny-1;
    indexLowY[ny-1]=ny-2;
}

Grid::~Grid(){
    
    delete[] this->x;
    delete[] this->dx;
    delete[] this->dx_x;
    delete[] this->difx_x_d;delete[] this->difx_x_up;delete[] this->difx_x_low;
    delete[] this->difx_c_d;delete[] this->difx_c_up;delete[] this->difx_c_low;
    delete[] this->dif2x_x_d;delete[] this->dif2x_x_up;delete[] this->dif2x_x_low;    
    delete[] this->dif2x_c_d;delete[] this->dif2x_c_up;delete[] this->dif2x_c_low;    
    delete[] this->y;
    delete[] this->indexUpY;
    delete[] this->indexLowY;
    delete[] this->x2c_d;delete[] this->x2c_up;
}

void Grid::calculateDifCoefficients(double dx0,double dx1,double & up, double & d, double & low){    
    double dx02=dx0*dx0;
	double dx12=dx1*dx1;
	double dx0dx1=dx0*dx1;
	double dx0_p_dx1=dx0+dx1;
	d=(dx12-dx02)/(dx0_p_dx1*dx0dx1);
	up=dx02/(dx0_p_dx1*dx0dx1);
	low=-dx12/(dx0_p_dx1*dx0dx1);
	
}

void Grid::calculateDif2Coefficients(double dx0,double dx1,double & up, double & d, double & low){
		double dx0dx1=dx0*dx1;
		double dx0_p_dx1=dx0+dx1;
		d=-2/(dx0dx1);
		up=2/(dx0_p_dx1*dx1);
		low=2/(dx0_p_dx1*dx0);
		
}

