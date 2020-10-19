/*********************************************************************************
 * 
 *               HydroBoundaryConditions.cpp (Morfo70) 
 * 
 *      General Methods por Hydro Class
 *          - Constructor and Destructor
 *          - Main methos for solving n iterations/time steps with the same wave forcing ->Run()
 *          - Resume state from resume file
 *          - Free surface computation
 * 
 *   Methods and parameters description in "Hydro.h"  
 * 
 * *******************************************************************************/

#include "Hydro.h"

Hydro::Hydro(Grid* grid,Bathymetry* bathy,Wave* wave){
    this->grid=grid;
    this->bathy=bathy;
    this->wave=wave;

    //Init vars
    int size=grid->size;
    int ny=grid->ny;
    this->U=new double[size+ny]();    
    this->Qx=new double[size+ny]();
    this->U_c=new double[size]();
    this->V=new double[size+ny]();    
    this->Qy=new double[size+ny]();
    this->V_c=new double[size]();
    this->dQx_dx=new double[size]();
    this->dQy_dy=new double[size]();
    this->Zs=new double[size]();
    
    //Get Zs from Zs=D-Zb
    for (int i=0;i<size;i++)
        Zs[i]=bathy->D[i]-bathy->Z[i]+bathy->h[i];
    

    /////////////////////////////////////////////////////////////////////////////
    ////  Set Boundary conditions parameters
    /////////////////////////////////////////////////////////////////////////////
    //Sponge
    const double lSponge=grid->dx_x[0]*5;
    const double exp_dx_2ls=exp(-0.5*grid->dx_x[0]/lSponge);
    const double exp_dx_ls=exp(-grid->dx_x[0]/lSponge);

    switch (Options::Hydro::OffshoreBC)
    {
        case Options::Hydro::OffshoreBCType::Sponge:
            this->bcOffshoreQx_f1=0.5*exp_dx_2ls/(1-0.5*exp_dx_2ls);
            this->bcOffshoreQx_f0=0;
            break;
        case Options::Hydro::OffshoreBCType::ConstantFreeSurface:
            this->bcOffshoreQx_f1=1;
            this->bcOffshoreQx_f0=0;
            break;
    }
    //For Qy projector
    this->bcOffshoreQy=exp_dx_ls;
    /////////////////////////////////////////////////////////////////////////////
    ////  END Set Boundary conditions parameters
    /////////////////////////////////////////////////////////////////////////////

    
    //Common vars one run (n iterations)
    this->cD=new double[size]();
    this->nuD=new double[size]();
    this->cD_x=new double[size-ny]();
    this->cD_y=new double[size]();
    this->Urms_x=new double[size-ny]();
    this->Urms_y=new double[size]();
    this->inv_Urms2_x=new double[size-ny]();
    this->inv_Urms2_y=new double[size]();
    this->nuD_x=new double[size-ny]();    
    this->nuD_y=new double[size]();    
    this->dnuD_dx=new double[size-ny]();
    this->dnuD_dy=new double[size]();
    this->dnuDx_dy=new double[size-ny]();
    this->dnuDy_dx=new double[size]();
    this->Forcing_x=new double[size]();
    this->Forcing_y=new double[size]();

    //Specific vars for each iteration / hydro time step
    this->guess_D=new double[size+ny]();
    this->guess_Zs=new double[size+ny]();
    this->Di=new double[size+ny]();  //D_x, D_y
    this->inv_Di=new double[size+ny]();    
    this->Ui_j=new double[size]();    //U_y, V_x
    this->DragTerm=new double[size]();   
    this->IsWater=new bool[size];

    
    //Solver matrix elements A=[low, d , up]
    this->d=new double[size]();
    this->f=new double[size]();
    this->up=new double[size]();
    this->low=new double[size]();
    this->v_tmp=new double[size]();


    //Init parametrizations methods  (Drag)      
    switch (Options::Hydro::DragModel) 
    {
        case (Options::Hydro::DragModelType::ConstantDrag):
            this->dragCoefficient=&Hydro::dragCoefficient_Constant; 
            break;
        case (Options::Hydro::DragModelType::LogDepth):
            this->dragCoefficient=&Hydro::dragCoefficient_LogDepth; 
            break;
        case (Options::Hydro::DragModelType::ManningStrickler):
            this->dragCoefficient=&Hydro::dragCoefficient_ManningStrickler; 
            break;
    }  

    //Init parametrizations methods  (Viscosity)      
    switch (Options::Hydro::ViscosityModel) 
    {
        case (Options::Hydro::ViscosityModelType::ConstantViscosity):
            this->viscosity=&Hydro::viscosity_Constant;
            break;
        case (Options::Hydro::ViscosityModelType::BreakingHrms):
            this->viscosity=&Hydro::viscosity_BreakingHrms;
            break;
        case (Options::Hydro::ViscosityModelType::BreakingDepth):
            this->viscosity=&Hydro::viscosity_BreakingDepth;
            break;
        case (Options::Hydro::ViscosityModelType::BreakingRollers):
            this->viscosity=&Hydro::viscosity_BreakingRollers;
            break;
    }  
    
}

Hydro::~Hydro()
{
    //State variables
    delete[] Qx;
    delete[] dQx_dx;
    delete[] Qy;
    delete[] dQy_dy;
    delete[] U;
    delete[] U_c;
    delete[] V;
    delete[] V_c;    
    delete[] Zs;
    
    


    //Common vars one run (n iterations)
    delete[] cD;delete[] nuD;
    delete[] cD_x;delete[] nuD_x; delete[] dnuD_dx;delete[] dnuDx_dy;
    delete[] cD_y;delete[] nuD_y; delete[] dnuD_dy;delete[] dnuDy_dx;
    delete[] Forcing_x; delete[] Urms_x;delete[] inv_Urms2_x;
    delete[] Forcing_y; delete[] Urms_y;delete[] inv_Urms2_y;
    
    //specific hydro time step vars
    delete[] Di; 
    delete[] guess_D;delete[] guess_Zs;
    delete[] Ui_j;
    delete[] inv_Di;
    delete[] DragTerm;
    delete[] IsWater;    
    delete[] d;delete[] up; delete[] low; delete[] f;   delete[]v_tmp;

}


void Hydro::Run(int nIter){
    int nx=this->grid->nx;
    int ny=this->grid->ny;    

    //Computes Forcing
    this->getForcing_x();
    this->getForcing_y();
    MeshOperations::C2X(this->Urms_x,this->wave->Urms, nx,ny);
    
    //Get drag coeficcient and viscosity
    this->C2Y(this->Urms_y,this->wave->Urms,nx,ny);    
    this->calculateDragAndViscosity();
    MeshOperations::C2X(this->cD_x,this->cD,nx,ny);
    
     this->C2Y(this->cD_y,this->cD,nx,ny);     
     this->difC2X(this->dnuD_dx,this->nuD,this->grid->dx,nx,ny);
     MeshOperations::C2X(this->nuD_x,this->nuD,nx,ny);
     this->C2Y(this->nuD_y,this->nuD,nx,ny);
     
     this->difC2C_y(this->dnuDx_dy,this->nuD_x,this->grid->dy,nx-1,ny);
     this->difC2Y(this->dnuD_dy,this->nuD,this->grid->dy,nx,ny);
     this->difC2C_x_addBoundary(this->dnuDy_dx,this->nuD_y,this->grid->
     difx_c_d,this->grid->difx_c_up,this->grid->difx_c_low,nx,ny);

    
     //Main loop 
     for (int k=0;k<nIter;k++)
     {        
          this->solveQx();        
          this->solveQy();
        
     }
}

void Hydro::setHydroFromResumeFile(string resumeFile){
    const int ny=this->grid->ny;
    

    MorfoIO::Netcdf ncResume=MorfoIO::Netcdf();
	 ncResume.Open(resumeFile,MorfoIO::Netcdf::OpenMode::Read);	 
	 ncResume.ReadVar("U",this->U_c);
     ncResume.ReadVarY("U0",this->U); //x=0, first column
     ncResume.ReadVarY("Qx0",this->Qx); //x=0, first column
     ncResume.ReadVar("Qx_x",this->Qx+ny);
     ncResume.ReadVar("U_x",this->U+ny);
     ncResume.ReadVar("V",this->V_c);
     ncResume.ReadVarY("V0",this->V); //x=0, first column
     ncResume.ReadVarY("Qy0",this->Qy);	 //x=0, first column
     ncResume.ReadVar("Qy_y",this->Qy+ny);
     ncResume.ReadVar("V_y",this->V+ny);
     ncResume.ReadVar("Zs",this->Zs);	      
     ncResume.ReadVar("D",this->bathy->D);	      

     ncResume.ReadVar("dQx_dx",this->dQx_dx);	      
     ncResume.ReadVar("dQy_dy",this->dQy_dy);	      
	 ncResume.Close();
    
}


void Hydro::solveZs_guess(int column){
    
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column;        
    double const  __restrict_arr * dqy_dy=this->dQy_dy+i*ny;        
    double const  __restrict_arr * dqx_dx=this->dQx_dx+i*ny;    
    double const   __restrict_arr * zs=this->Zs+i*ny;
    double const  __restrict_arr * d=this->bathy->D+i*ny;
    double   __restrict_arr * guess_zs=this->guess_Zs+i*ny;
    double   __restrict_arr * guess_d=this->guess_D+i*ny;
    for (int j=0;j<ny;j++)
    {
           
            const double incrZs=-dt*(dqx_dx[j]+dqy_dy[j]);
            guess_zs[j]=zs[j]+incrZs;
            guess_d[j]=d[j]+incrZs;
    }
    
    
}


void Hydro::solveZs(int column){
    
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column;        
    double const  __restrict_arr * dqy_dy=this->dQy_dy+i*ny;        
    double const  __restrict_arr * dqx_dx=this->dQx_dx+i*ny;    
    double   __restrict_arr * zs=this->Zs+i*ny;
    double   __restrict_arr * d=this->bathy->D+i*ny;
    for (int j=0;j<ny;j++)
    {
           
            double incrZs=-dt*(dqx_dx[j]+dqy_dy[j]);
            zs[j]+=incrZs;
            d[j]+=incrZs;
    }
    
    
}



