/*******************************************************************************************
*  Wave Class:  stores wave field and computes wave and rollers propagation of the entire 
*  domain from offshore to onshore
*
*   Methods and parameters description in "Wave.h" 
*
*******************************************************************************************/

#include "Wave.h"
#include <omp.h>


Wave::Wave(Grid *grid, Bathymetry *bathy)
{
    //Reference Grid and Bathymetry vars
    this->ny = grid->ny;
    this->nx = grid->nx;
    this->x = grid->x;
    this->dx = grid->dx;
    this->dy = grid->dy;
    this->D = bathy->D;
    this->size=grid->size;
    this->indexUpY=grid->indexUpY;
    this->indexLowY=grid->indexLowY;

    //Initialice wave vars
    this->K = new double[grid->size]; //NABIL TMP -> this is to storage log(K)
    this->Cg = new double[grid->size]();
    this->E = new double[grid->size]();
    this->H = new double[grid->size]();
    this->Theta = new double[grid->size]();
    this->Sxx_w = new double[grid->size]();
    this->Sxy_w = new double[grid->size]();
    this->Syy_w = new double[grid->size]();
    this->Dw = new double[grid->size]();
    this->U = new double[grid->size]();
    this->V = new double[grid->size]();
    this->U_w = new double[grid->size]();

    //Initialice wave rollers
    this->E_r = new double[grid->size]();
    this->D_r = new double[grid->size]();
    this->Sxx_r = new double[grid->size]();
    this->Sxy_r = new double[grid->size]();
    this->Syy_r = new double[grid->size]();
    
    this->Sxx = new double[grid->size]();
    this->Sxy = new double[grid->size]();
    this->Syy = new double[grid->size]();
    this->Urms = new double[grid->size]();

    this->CosTheta=new double[grid->size]();
    this->SinTheta=new double[grid->size]();

   
}


Wave::~Wave()
{
    delete[] K; //NABIL TMP    
    delete[] Cg;
    delete[] E;
    delete[] H;
    delete[] Theta;
    delete[] Sxx_w;delete[] Sxy_w; delete[] Syy_w;
    delete[] Dw;
    delete[] U;
    delete[] V;
    delete[] U_w;
    delete[] E_r;delete[] D_r;
    delete[] Sxx_r;delete[] Sxy_r; delete[] Syy_r;
    
    delete[] Sxx;delete[] Sxy; delete[] Syy;
    delete[] CosTheta;delete[] SinTheta;
    
}


void Wave::calculateWaveField(double hOff, double thetaOff, double T)
{

    const int ny = this->ny;
    const int nx = this->nx;

    ///////////////////////////////////////////////////////
    ////////OPEN MP (only works if it was compiled with DUSE_OMP=1)
    ///////////////////////////////////////////////////////
    //Number of processors if OMP flag was used, otherwise (serial) returns 1
    int np=Omp::get_num_processors();  
    
    int chunk_base = ny / np;
    int remain = ny - chunk_base * np;
    int *chunk = new int[np];
    int *indexIni = new int[np]();
    for (int i = 0; i < np; i++)
    {
        chunk[i] = chunk_base;

        if (remain > 0)
            chunk[i] += 1;
        remain -= 1;
    }
    for (int i = 1; i < np; i++)
        indexIni[i] += indexIni[i - 1] + chunk[i - 1];
    ///////////////////////////////////////////////////////
    ////////End OPEN MP
    ///////////////////////////////////////////////////////


    //x=0 
    double x0 = x[0]; //offshore x-position
    int xIndex = 0; //offshore x-index-position
    WaveRow* row = new WaveRow(ny, dy, indexUpY, indexLowY); //WaveRow Object to compute
    WaveRow* row_m1 = new WaveRow(ny, dy, indexUpY, indexLowY); //WaveRow Object of the x-upper row
    WaveRow* row_tmp=nullptr; //Dummy variable  for allow switch rows

    //Init and store data of the first row
    row_m1->initFirstRow(hOff, thetaOff, T, D, U, V, U_w); 
    storageRow(0, row_m1);

    double dx_row;
    bool storeResults;    
    
    while (xIndex < nx - 1)
    {
      
        storeResults = getDx(x0, xIndex, row_m1, dx_row);
       

        int index = (xIndex + 1) * ny;
        int index_m1 = xIndex * ny;
        if (storeResults)  //row is a mesh-row
            row->initRow(row_m1,dx_row, D + index, U + index, V + index, U_w + index);
        else //row is located in a intermediate position between mesh-rows. Interpolation needed
            row->initRow(row_m1,dx_row, dx[xIndex], D + index, U + index, V + index, U_w + index, D + index_m1, U + index_m1, V + index_m1, U_w + index_m1);

        //solve dispersion relationship
        row->dispersion_parallel(T, indexIni, chunk, np);                
        
        //get the ray refraction
        row->solveDir_parallel(row_m1, indexIni, chunk, np);
        
        //compute constan termo for energy solver
        row->energyPreprocessor(row_m1);
        
        //Solver Energy
        row->solveEnergy_parallel(indexIni, chunk, np);
        
        if (Options::Rollers::CalculateRollers)
             row->solveRollers(row_m1);        

        if (storeResults)
        {
            xIndex = xIndex + 1;
            storageRow(xIndex, row);        
        }
        else                

        //WaveRow pointint to next position
        row_tmp = row_m1;
        row_m1 = row;
        row = row_tmp;
        x0 = x0 + dx_row; //update x position        
        
    }

    delete(row_m1);
    delete(row);
    delete[] chunk;
    delete[] indexIni;
}


void Wave::storageRow(const int xIndex, const WaveRow *r)
{
    int ny = this->ny;
    int iniIndex = xIndex * ny;
    std::copy(r->Theta, r->Theta + ny, this->Theta + iniIndex);
    std::copy(r->CosTheta, r->CosTheta + ny, this->CosTheta + iniIndex);
    std::copy(r->SinTheta, r->SinTheta + ny, this->SinTheta + iniIndex);
    std::copy(r->Cg, r->Cg + ny, this->Cg + iniIndex);
    std::copy(r->H, r->H + ny, this->H + iniIndex);
    std::copy(r->Dw, r->Dw + ny, this->Dw + iniIndex);

    std::copy(r->lK, r->lK + ny, this->K + iniIndex); // NABIL TMP

    double *E = new double[ny];
    std::copy(r->E, r->E + ny, E);
    double *Sxx_w = new double[ny];
    std::copy(r->sxx, r->sxx + ny, Sxx_w);
    double *Sxy_w = new double[ny];
    std::copy(r->sxy, r->sxy + ny, Sxy_w);
    double *Syy_w = new double[ny];
    std::copy(r->syy, r->syy + ny, Syy_w);

    double *E_r = new double[ny];
    std::copy(r->E_r, r->E_r + ny, E_r);
    double *D_r = new double[ny];
    std::copy(r->d_r, r->d_r + ny, D_r);
    double *Sxx_r = new double[ny];
    std::copy(r->sxx_r, r->sxx_r + ny, Sxx_r);
    double *Sxy_r = new double[ny];
    std::copy(r->sxy_r, r->sxy_r + ny, Sxy_r);
    double *Syy_r = new double[ny];
    std::copy(r->syy_r, r->syy_r + ny, Syy_r);
    double* Sxx=new double[ny];
    double* Sxy=new double[ny];
    double* Syy=new double[ny];
    double* Urms=new double[ny];
    std::copy(r->Urms_H, r->Urms_H + ny, Urms);
    

    for (int i = 0; i < ny; i++)
    {
        
        double e = E[i];
        Sxx_w[i] *= e;
        Sxy_w[i] *= e;
        Syy_w[i] *= e;

        double e_r = E_r[i];
        Sxx_r[i] *= e_r;
        Sxy_r[i] *= e_r;
        Syy_r[i] *= e_r;
        D_r[i] *= e_r;

        Sxx[i]=Sxx_w[i]+Sxx_r[i];
        Sxy[i]=Sxy_w[i]+Sxy_r[i];
        Syy[i]=Syy_w[i]+Syy_r[i];
        Urms[i]*=r->H[i];
        

    
    }

    std::copy(E, E + ny, this->E + iniIndex);
    std::copy(Sxx_w, Sxx_w + ny, this->Sxx_w + iniIndex);
    std::copy(Sxy_w, Sxy_w + ny, this->Sxy_w + iniIndex);
    std::copy(Syy_w, Syy_w + ny, this->Syy_w + iniIndex);

    std::copy(E_r, E_r + ny, this->E_r + iniIndex);
    std::copy(D_r, D_r + ny, this->D_r + iniIndex);
    std::copy(Sxx_r, Sxx_r + ny, this->Sxx_r + iniIndex);
    std::copy(Sxy_r, Sxy_r + ny, this->Sxy_r + iniIndex);
    std::copy(Syy_r, Syy_r + ny, this->Syy_r + iniIndex);

    std::copy(Sxx, Sxx + ny, this->Sxx + iniIndex);
    std::copy(Sxy, Sxy + ny, this->Sxy + iniIndex);
    std::copy(Syy, Syy + ny, this->Syy + iniIndex);

    std::copy(Urms, Urms + ny, this->Urms + iniIndex);

    
    

    delete[] E;
    delete[] Sxx_w;delete[] Sxy_w;delete[] Syy_w;
    delete[] E_r;
    delete[] Sxx_r; delete[] Sxy_r; delete[] Syy_r;
    delete[] D_r;
    delete[] Sxx;delete[] Sxy;delete[] Syy;
    delete[] Urms;
}


bool Wave::getDx(const double x0, const int xIndex, const WaveRow *rm1, double &dx_row)
{

    int ny = this->ny;
    double dy = this->dy;
    double *tanTheta_rm1 = new double[ny];
    double _x = this->x[xIndex+1];
    std::copy(rm1->TanTheta, rm1->TanTheta + ny, tanTheta_rm1);
    

    bool storeResults = true;
    double tanThetaMax = -999;
    double tanTheta;

    //Get the maximun value of tan(Theta) of the incdente
    for (int i = 0; i < ny; i++)
    {
        tanTheta = abs(tanTheta_rm1[i]);
        if (tanTheta > tanThetaMax)
            tanThetaMax = tanTheta;
    }

    
    double dx_remain = _x - x0; //distance vs the actual postion and next mesh row
    if (tanThetaMax == 0)
    {
        dx_row = dx_remain;
        delete[] tanTheta_rm1;     
        return true;
    }

    dx_row = dy / tanThetaMax; //minumun distante requires to solve wave dir

    if (dx_row < dx_remain)
    {
        storeResults = false;
        dx_row = dx_row * 0.9;
    }
    else
    {
        storeResults = true;
        dx_row = dx_remain;
    }
    delete[] tanTheta_rm1;
    

    return storeResults;
}


void Wave::setUV(double __restrict_arr const  * _U,double __restrict_arr const * _V)
{
    
    std::copy(_U,_U+size,this->U);
    std::copy(_V,_V+size,this->V);
    for (int i=0;i<nx;i++)    
        for (int j=0;j<ny;j++)        
            U_w[i*ny+j]=CosTheta[i*ny+j]*U[i*ny+j]+SinTheta[i*ny+j]*V[i*ny+j];        
            
    
}


void Wave::setUVFromResumeFile(string resumeFile)
{
     MorfoIO::Netcdf ncResume=MorfoIO::Netcdf();
	 ncResume.Open(resumeFile,MorfoIO::Netcdf::OpenMode::Read);
	 ncResume.ReadVar("Theta",this->Theta);
	 ncResume.ReadVar("U",this->U);
     ncResume.ReadVar("V",this->V);	 
	 ncResume.Close();

	 for (int i=0;i<nx;i++)    
        for (int j=0;j<ny;j++)        
                U_w[i*ny+j]=cos(Theta[i*ny+j])*U[i*ny+j]+sin(Theta[i*ny+j])*V[i*ny+j];        
	
}


void Wave::setWaveEnergyFromResumeFile(string resumeFile){
     MorfoIO::Netcdf ncResume=MorfoIO::Netcdf();
	 ncResume.Open(resumeFile,MorfoIO::Netcdf::OpenMode::Read);
     ncResume.ReadVar("E",this->E);	 	 
     ncResume.ReadVar("H",this->H);	 	 
	 ncResume.ReadVar("Theta",this->Theta);	 	 
	 ncResume.Close();
}
