/*********************************************************************************
 * 
 *                  Sediment.cpp (Morfo70)
 * 
 *  Sediment transport and bottom evolution implementation
 *
 *  Methods and parameters description in "Sediment.h"  
 * 
 * *******************************************************************************/

#include "Sediment.h"

Sediment::Sediment(Hydro* hydro,Wave* wave,Bathymetry* bathy,Grid* grid)
{
    this->grid=grid;
    this->bathy=bathy;
    this->wave=wave;
    this->hydro=hydro;
    int size=grid->size;

    this->qx=new double[size]();
    this->qy=new double[size]();
    this->dh_dt=new double[size]();

    this->Alpha=new double[size]();
    this->Gamma=new double[size]();
    this->gamma=new double[size]();

    this->maxDw=new double[grid->ny]();
    this->Zs_sh=new double[grid->ny]();
    
}

Sediment::~Sediment()
{
    delete[] Alpha;
    delete[] Gamma;
    delete[] gamma;
    delete[] qx;
    delete[] qy;
    delete[] dh_dt;
    delete[] maxDw;
    delete[] Zs_sh;
    
}

void Sediment::calculateSedimentFluxes()
{
    const int nx=this->grid->nx;
    //Calculo los paranetros comunes para la zona de swash
    this->getMaxWaveBreakingDissipation();
    this->getShorelineZs();

    this->getWaveStirringCWS(0);
    this->getGammaDownslopeCWS(0);
    this->getGammaSwash(0);
    this->getqy(0);

    #pragma omp parallel for
    for (int i=1;i<nx;i++)
    {
        this->getWaveStirringCWS(i);
        this->getGammaDownslopeCWS(i);
        this->getGammaSwash(i);
        this->getqx(i-1);
        this->getqy(i);
    }
}

void Sediment::updateBottom()
{
    const int nx=this->grid->nx;

    #pragma omp parallel for
    for (int i=1;i<nx;i++)    
        this->updateBottomColumn(i);
    
}

void Sediment::setSedimentFluxesFromResumeFile(string resumeFile){
    MorfoIO::Netcdf ncResume=MorfoIO::Netcdf();
	 ncResume.Open(resumeFile,MorfoIO::Netcdf::OpenMode::Read);	 
	 //ncResume.ReadVar("dh_dt",this->dh_dt);
     ncResume.ReadVar("qx",this->qx);
     ncResume.ReadVar("qy",this->qy);
     ncResume.Close();
     
     const int ny=this->grid->ny;
     const int nx=this->grid->nx;
    const double p=Parameters::Bottom::p;    
    const double _1_m1p=1/(1-p);
    const double inv_dy=1/this->grid->dy;          

     for (int i=1;i<nx;i++)
     {
        const double inv_dx=1/this->grid->dx_x[i];      //[i+1]??
        double dqy_dy=inv_dy*(qy[i*ny]-qy[i*ny+ny-1]);        
        double dqx_dx=inv_dx*(qx[i*ny]-qx[(i-1)*ny]);        
        this->dh_dt[i*ny]=-_1_m1p*(dqx_dx+dqy_dy);          
        for (int j=1;j<ny;j++)
        {
            double dqy_dy=inv_dy*(qy[i*ny+j]-qy[i*ny+j-1]);        
            double dqx_dx=inv_dx*(qx[i*ny+j]-qx[(i-1)*ny+j]);        
            this->dh_dt[i*ny+j]=-_1_m1p*(dqx_dx+dqy_dy);   
        }
        
     }
     
    
}

void Sediment::updateBottomColumn(const int column)
{    
    const int ny=this->grid->ny;
    const int i=column;             

    const double dt_2=0.5*Parameters::Run::dt;

    const double inv_dx=1/this->grid->dx_x[i];    //[i+1]??
    const double inv_dy=1/this->grid->dy;        
    const double p=Parameters::Bottom::p;    
    const double _1_m1p=1/(1-p);

    const double  __restrict_arr *_qx=this->qx+(i-1)*ny; //El aterior
    const double  __restrict_arr *_qy=this->qy+i*ny;
    double  __restrict_arr *_h=this->bathy->h+i*ny;
    double  __restrict_arr *_Zb=this->bathy->Zb+i*ny;
    double  __restrict_arr *_Zs=this->hydro->Zs+i*ny;
    double  __restrict_arr *_dh_dt_m1=this->dh_dt+i*ny;

    double dqy_dy=inv_dy*(_qy[0]-_qy[ny-1]);        
    double dqx_dx=inv_dx*(_qx[ny]-_qx[0]);        
    double _dh_dt=-_1_m1p*(dqx_dx+dqy_dy);        
    double incr_h=dt_2*(3*_dh_dt-_dh_dt_m1[0]);  //Explicit solution Adams Bashford (2)
    _h[0]+=incr_h;
    _Zb[0]-=incr_h;
    _Zs[0]+=incr_h;
    _dh_dt_m1[0]=_dh_dt;

    for (int j=1;j<ny;j++)
    {
        dqy_dy=inv_dy*(_qy[j]-_qy[j-1]);        
        dqx_dx=inv_dx*(_qx[j+ny]-_qx[j]);        
        _dh_dt=-_1_m1p*(dqx_dx+dqy_dy);        
        incr_h=dt_2*(3*_dh_dt-_dh_dt_m1[j]);  //Explicit solution Adams Bashford (2)
        _h[j]+=incr_h;
        _Zb[j]-=incr_h;
        _Zs[j]+=incr_h;
        _dh_dt_m1[j]=_dh_dt;
    }
    


}

void Sediment::getqx(const int column){
    const int ny=this->grid->ny;
    const int i=column;             
    const double inv_dx=1/this->grid->dx[i];        
    const double  __restrict_arr *_g=this->gamma+i*ny;
    const double  __restrict_arr *_Alpha=this->Alpha+i*ny;
    const double  __restrict_arr *_G=this->Gamma+i*ny;
    const double  __restrict_arr *_h=this->bathy->h+i*ny;
    const double  __restrict_arr *_U=this->hydro->U+(i+1)*ny;
    double  __restrict_arr *_qx=this->qx+i*ny;

    for (int j=0;j<ny;j++)
    {
        double A_x=0.5*(_Alpha[j]+_Alpha[j+ny]);
        double g_x=0.5*(_g[j]+_g[j+ny]);
        double G_x=0.5*(_G[j]+_G[j+ny]);
        double dh_dx=(_h[j+ny]-_h[j])*inv_dx;
        _qx[j]=A_x*(_U[j]-g_x*dh_dx)-G_x*dh_dx;
    }
}

void Sediment::getqy(const int column){
    const int ny=this->grid->ny;
    const int i=column;             
    const double inv_dy=1/this->grid->dy;        
    const double  __restrict_arr *_g=this->gamma+i*ny;
    const double  __restrict_arr *_Alpha=this->Alpha+i*ny;
    const double  __restrict_arr *_G=this->Gamma+i*ny;
    const double  __restrict_arr *_h=this->bathy->h+i*ny;
    const double  __restrict_arr *_V=this->hydro->V+(i+1)*ny;
    double  __restrict_arr *_qy=this->qy+i*ny;

    for (int j=0;j<ny-1;j++)
    {
        double A_y=0.5*(_Alpha[j]+_Alpha[j+1]);
        double g_y=0.5*(_g[j]+_g[j+1]);
        double G_y=0.5*(_G[j]+_G[j+1]);
        double dh_dy=(_h[j+1]-_h[j])*inv_dy;
        _qy[j]=A_y*(_V[j]-g_y*dh_dy)-G_y*dh_dy;
    }
    double A_y=0.5*(_Alpha[ny-1]+_Alpha[0]);
    double g_y=0.5*(_g[ny-1]+_g[0]);
    double G_y=0.5*(_G[ny-1]+_G[0]);
    double dh_dy=(_h[0]-_h[ny-1])*inv_dy;
    _qy[ny-1]=A_y*(_V[ny-1]-g_y*dh_dy)-G_y*dh_dy;
}

void Sediment::getGammaDownslopeCWS(const int column)
{
    const int ny=this->grid->ny;
    const int i=column;             
    const double gamma0_tanPhi=Parameters::Sediment::CWS::gamma0*inv_tanPhi_i;
    double  __restrict_arr *_g=this->gamma+i*ny;
    const  double  __restrict_arr *_d=this->bathy->D+i*ny;
    const double __restrict_arr * zb=this->bathy->Zb+i*ny;
    const double __restrict_arr * zs_sh=this->Zs_sh;    
    const double Dmin=Parameters::Hydro::Dmin;
    const double A_swash=Parameters::Sediment::Swash::A_swash;
    

    //const double _2_Doff_on=2/(A_swash-Dmin);
    const double _2_Doff_on=2/(A_swash);

    for (int j=0;j<ny;j++)   
    {
        const int isDry=(_d[j]<=Dmin);
        const double d_dry=(1-isDry)*_d[j]+isDry*(zb[j]+zs_sh[j]);        
        //const double sf_drySwash=(tanh(pi*(_2_Doff_on*(_d[j]-Dmin)-1))+1)*0.5;  
        const double sf_drySwash=(tanh(pi*(_2_Doff_on*(d_dry+A_swash)-1))+1)*0.5; 
        //_g[j]=gamma0_tanPhi*sf_drySwash;

        //Segun EGU
        const double u_bsl=(tanh(pi*(2*(_d[j]-Dmin)/(A_swash-Dmin)-1))+1)/2; 
        _g[j]=u_bsl*gamma0_tanPhi;
    }
    

    
}

void Sediment::getWaveStirringCWS(const int column)
{
    const int ny=this->grid->ny;
    const int i=column;       
    const double alpha0_min=Parameters::Sediment::CWS::alpha0_min;
    const double incrAlpha=Parameters::Sediment::CWS::alpha0_max-alpha0_min;
    const double Dmax=Parameters::Sediment::CWS::D_alpha0_max;
    const double Dmax_min=Parameters::Sediment::CWS::D_alpha0_min-Dmax;

    const double __restrict_arr * d=this->bathy->D+i*ny;    

    double  __restrict_arr *_alpha=this->Alpha+i*ny;

    for (int j=0;j<ny;j++)
    {
         double x=1+2*(Dmax-d[j])/Dmax_min;
        _alpha[j]=alpha0_min+incrAlpha*0.5*(1+tanh(pi*x));

    }

}

void Sediment::getGammaSwash(const int column)
{
    
    const int i=column;
    const int ny=this->grid->ny;
    const double G0=Parameters::Sediment::Swash::Gamma0;    
    const double Dmin=Parameters::Hydro::Dmin+1E-4;
    const double A_swash=Parameters::Sediment::Swash::A_swash;
    const double inv_rho_g_sm1=1/(g*rho*(s-1));
    const double __restrict_arr * _d=this->bathy->D+i*ny;
    const double __restrict_arr * zb=this->bathy->Zb+i*ny;
    const double __restrict_arr * zs_sh=this->Zs_sh;    
    const double __restrict_arr * _maxDw=this->maxDw;    
    double __restrict_arr * _G=this->Gamma+i*ny;
    
    const double _2_Doff_on=2/(A_swash-Dmin);

    for (int j=0;j<ny;j++)
    {
        const int isDry=(_d[j]<=Dmin);
        const double d_dry=(1-isDry)*_d[j]+isDry*(abs(zb[j]+zs_sh[j]));
        const double sf_swash=1-(tanh(pi*(_2_Doff_on*(d_dry-Dmin)-1))+1)*0.5;
        _G[j]=G0*inv_rho_g_sm1*sf_swash*_maxDw[j]; //TODO
        //_G[j]=sf_swash*G0*Parameters::Sediment::CWS::alpha0_max; //TODO CTE
    }


}

void Sediment::getMaxWaveBreakingDissipation(){
    const int ny=this->grid->ny;
    const int nx=this->grid->nx;    
    const double __restrict_arr * dw=this->wave->Dw;
    double __restrict_arr * _maxDw=this->maxDw;

    for (int i=0;i<nx;i++)
        for (int j=0;j<ny;j++)
            _maxDw[j]=max(_maxDw[j],dw[i*ny+j]);

}

void Sediment::getShorelineZs(){
    const int ny=this->grid->ny;
    const int nx=this->grid->nx;    
    const double __restrict_arr * zs=this->hydro->Zs;
    const double __restrict_arr * d=this->bathy->D;
    double __restrict_arr * zs_sh=this->Zs_sh;
    const double Dmin=Parameters::Hydro::Dmin+1E-4;
    int * isSh=new int[ny]();
    for (int i=1;i<nx;i++)
        for (int j=0;j<ny;j++)
        {            
            const int isDry=(d[i*ny+j]<=Dmin);
            const int isDryLow=(d[(i-1)*ny+j]<=Dmin);
            isSh[j]=!(isDry*isDryLow);
            zs_sh[j]=(1-isSh[j])*zs_sh[j]+zs[i*ny+j]*isSh[j];

        }
    
    delete[] isSh;        

}