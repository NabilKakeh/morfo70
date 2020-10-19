/*******************************************************************************************
*  WaveRow Class   stores and computes wave refraction in a row (Xrow, Y). 
*  IndexUp and IndexLow variables are used to compute faster peridic boundary conditions.
*
*  Methods and parameters description in "WaveRow.h" 
*
*  WaverRow.ccp contains constructor, desctructor and initialize methods. Other methods are 
*  implemented in:
    - WaveEnergy.cpp 
    - WaveDir.cpp 
    - WaveRollers.cpp 
    - WaveParametrizations.cpp
*
*******************************************************************************************/

#include "WaveRow.h"


WaveRow::WaveRow(int ny, double dy, int *indexUpY, int *indexLowY)
{
    this->ny = ny;
    D = new double[ny];   
    lK = new double[ny];    
    Urms_H=new double[ny];
    Theta = new double[ny];
    CosTheta = new double[ny];
    SinTheta = new double[ny];
    TanTheta = new double[ny];
    Sigma = new double[ny];
    U = new double[ny];
    V = new double[ny];
    U_w = new double[ny];
    C = new double[ny];
    Cg = new double[ny];
    inv_S_i = new double[ny];
    Theta_i = new double[ny];
    SinTheta_i = new double[ny];
    CosTheta_i = new double[ny];
    Y_i = new double[ny];       
    df=new double[ny]();
    dw=new double[ny]();
    ddf_dh=new double[ny]();
    ddw_dh=new double[ny]();
    sxx=new double[ny]();
    sxy=new double[ny]();
    syy=new double[ny]();
    sxx_r=new double[ny]();
    sxy_r=new double[ny]();
    syy_r=new double[ny]();
    lE=new double[ny]();
    H=new double[ny]();
    E=new double[ny]();
    Dw=new double[ny]();
    G=new double[ny]();
    Gp=new double[ny]();    
    WC=new double[ny]();
    WC_r=new double[ny]();
    E_r=new double[ny]();
    D_r=new double[ny]();
    d_r=new double[ny]();    
     Cx=new double[ny]();
    Cy=new double[ny]();
     Cgx=new double[ny]();
    Cgy=new double[ny]();    
    


    this->dy = dy;
    this->inv_2dy = 0.5 / dy;

    this->indexUpY = indexUpY;
    this->indexLowY = indexLowY;
    
    // Set the Wave Breking model
    switch (Options::Wave::BreakingModel)
     {
         case (Options::Wave::BreakingModelType::ChurchThornton):
             this->waveBreakingDissipation=&WaveRow::waveBreakingDissipationCT;        
             break;
        case (Options::Wave::BreakingModelType::ThorntonGuza):
             this->waveBreakingDissipation=&WaveRow::waveBreakingDissipationTG;        
             break;
        case (Options::Wave::BreakingModelType::RegularWaves):
             this->waveBreakingDissipation=&WaveRow::waveBreakingDissipationRegular;        
             break;
        case (Options::Wave::BreakingModelType::None):
            break;
     }

     // Set the model for the wave skin friction
     switch (Options::Wave::fwModel)
    {
         case (Options::Wave::FricCoefficientModelType::Constant):
             this->BottomFriction=&WaveRow::bottomFrictionConstant;
             break;
        case (Options::Wave::FricCoefficientModelType::Soulsby97):
             this->BottomFriction=&WaveRow::bottomFrictionSoulsby97;
             break;
     }
}


WaveRow::~WaveRow()
{
    delete[] Theta; delete[] CosTheta;delete[]SinTheta;delete[] TanTheta;
    delete[] D;delete[]lK;delete[] Urms_H;    
    delete[]U;delete[]V;delete[]U_w;
    delete[]Sigma;delete[]C;delete[]Cg;        
    delete[] inv_S_i;delete[]SinTheta_i;delete[]CosTheta_i;delete[]Theta_i;delete[]Y_i;   
    
    
    delete[] G;    delete[]  Gp;    
    delete[] WC;delete[] WC_r;
    delete[] dw;delete[] ddw_dh;
    delete[] df;     delete[] ddf_dh;
    delete[] sxx;delete[] sxy;delete[] syy;    
    delete[] lE;delete[] H;delete[] E;  delete[] Dw;
    delete[] sxx_r;delete[] sxy_r;delete[] syy_r;
    delete[] E_r;delete[] D_r;delete[] d_r;    
    delete[] Cx;delete[] Cy;
    delete[] Cgx;delete[] Cgy;

    
    
    
 
}



void WaveRow::initFirstRow(const double hOff, const double thetaOff, const double T, const double *D, const double *U, const double *V, const double *U_w)
{

    int ny=this->ny;    
    std::copy(D,D+ny,this->D);
    std::copy(U,U+ny,this->U);
    std::copy(V,V+ny,this->V);
    std::copy(U_w,U_w+ny,this->U_w);        
    
    //Calulate parametrizations for first node
    const double w=2*pi/T;     
    dispersion(0,ny,w);
    int nodes=0;
    double dw0=0;
    double ddw_dh0=0;
    double df0=0;
    double ddf_dh0=0;    
    
    
    (this->*waveBreakingDissipation)(1,&hOff, &w,D,&nodes,1,&dw0,&ddw_dh0);
     if (Options::Wave::BottomFriction)
         (this->*BottomFriction)(1,&hOff,w,Urms_H,&nodes,1,&df0,&ddf_dh0);
         

    //Initialize vars
    double* Theta=new double[ny];    
    double* CosTheta=new double[ny];    
    double* SinTheta=new double[ny];    
    double* TanTheta=new double[ny];        
    double* Cgx=new double[ny];    
    double* Cgy=new double[ny];    
    double* Cx=new double[ny];    
    double* Cy=new double[ny];        
    double* H=new double[ny];
    double* E=new double[ny];
    double* Dw=new double[ny];
    double* lE=new double[ny];
    double* dw=new double[ny];
    double* df=new double[ny];
    double* sxx=new double[ny];
    double* sxy=new double[ny];
    double* syy=new double[ny];
    double* sxx_r=new double[ny];
    double* sxy_r=new double[ny];
    double* syy_r=new double[ny];

    

    
    //Set vars for all nodes
    double EOff=hOff*hOff*cteEnergy;        
    double lEOff=log(EOff);        
    #pragma omp parallel for    
    for (int i = 0; i < ny; i++)
    {
        
        Theta[i] = thetaOff;
        double cosTheta=cos(thetaOff);
        double sinTheta=sin(thetaOff);
        CosTheta[i] = cosTheta;
        SinTheta[i] = sinTheta;
        TanTheta[i] = sinTheta/ cosTheta;            
        Cgx[i]=Cg[i]*cosTheta;
        Cgy[i]=Cg[i]*sinTheta;
        Cx[i]=C[i]*cosTheta;
        Cy[i]=C[i]*sinTheta;
        double cg_c=Cg[i]/C[i];        
        double cosTheta2=cosTheta*cosTheta;
        double cossinTheta=sinTheta*cosTheta;
        sxx[i]=((1+cosTheta2)*cg_c-0.5);
        syy[i]=((2-cosTheta2)*cg_c-0.5);
        sxy[i]=(cossinTheta*cg_c);
        sxx_r[i]=2*cosTheta2;
        syy_r[i]=2*(1-cosTheta2);
        sxy_r[i]=2*cossinTheta;
        H[i]=hOff;
        E[i]=EOff;        
        Dw[i]=dw0*EOff;
        lE[i]=lEOff;
        dw[i]=dw0;
        df[i]=df0;
    }
    
    std::copy(Theta,Theta+ny,this->Theta);
    std::copy(CosTheta,CosTheta+ny,this->CosTheta);
    std::copy(SinTheta,SinTheta+ny,this->SinTheta);
    std::copy(TanTheta,TanTheta+ny,this->TanTheta);
    std::copy(sxx,sxx+ny,this->sxx);
    std::copy(syy,syy+ny,this->syy);
    std::copy(sxy,sxy+ny,this->sxy);
    std::copy(sxx_r,sxx_r+ny,this->sxx_r);
    std::copy(syy_r,syy_r+ny,this->syy_r);
    std::copy(sxy_r,sxy_r+ny,this->sxy_r); 
     std::copy(Cgx,Cgx+ny,this->Cgx);
     std::copy(Cgy,Cgy+ny,this->Cgy);
     std::copy(Cx,Cx+ny,this->Cx);
     std::copy(Cy,Cy+ny,this->Cy);
    std::copy(H,H+ny,this->H);
    std::copy(E,E+ny,this->E);
    std::copy(Dw,Dw+ny,this->Dw);
    std::copy(lE,lE+ny,this->lE);
    std::copy(dw,dw+ny,this->dw);
    std::copy(df,df+ny,this->df);
    delete[] Theta;delete[] CosTheta;delete[] SinTheta;delete[] TanTheta;
    delete[] Cgx;delete[] Cgy;
    delete[] Cx;delete[] Cy;
    delete[] H;delete[] lE;delete[] E;delete[] dw;delete[] df;
    delete[] Dw;
    delete[]sxx;delete[]sxy;delete[]syy;
    delete[]sxx_r;delete[]sxy_r;delete[]syy_r;
}


void WaveRow::initRow(const WaveRow* row_m1,const double dx, const double *D, const double *U, const double *V,  const double *U_w)
{
    int ny=this->ny;
    std::copy(D,D+ny,this->D);
    std::copy(U,U+ny,this->U);
    std::copy(V,V+ny,this->V);
    std::copy(U_w,U_w+ny,this->U_w);
    std::copy(row_m1->H,row_m1->H+ny,this->H);
    std::copy(row_m1->lE,row_m1->lE+ny,this->lE);
    this->dx = dx;
}

void WaveRow::initRow(const WaveRow* row_m1,const double dx, const double dx_grid, const double *D_p1, const double *U_p1, const double *V_p1,  const double *U_w_p1,const double *D_m1,const double *U_m1, const double *V_m1, const double *U_w_m1)
{
    this->dx = dx;
    double dx_dx = dx / dx_grid;
    int ny=this->ny;
    double* d_p1=new double[ny];
    double* d_m1=new double[ny];
    double* u_p1=new double[ny];
    double* u_m1=new double[ny];
    double* v_p1=new double[ny];
    double* v_m1=new double[ny];
    double* u_w_p1=new double[ny];
    double* u_w_m1=new double[ny];    
    std::copy(D_p1,D_p1+ny,d_p1);
    std::copy(D_m1,D_m1+ny,d_m1);    
    std::copy(U_p1,U_p1+ny,u_p1);
    std::copy(U_m1,U_m1+ny,u_m1);    
    std::copy(V_p1,V_p1+ny,v_p1);
    std::copy(V_m1,V_m1+ny,v_m1);    
    std::copy(U_w_p1,U_w_p1+ny,u_w_p1);
    std::copy(U_w_m1,U_w_m1+ny,u_w_m1);    
    
    std::copy(row_m1->H,row_m1->H+ny,this->H);
    std::copy(row_m1->lE,row_m1->lE+ny,this->lE);

    #pragma omp parallel for
    for (int i = 0; i < ny; i++)
    {
        d_m1[i] += (d_p1[i] - d_m1[i]) * dx_dx;
        u_m1[i] += (u_p1[i] - u_m1[i]) * dx_dx;
        v_m1[i] += (v_p1[i] - v_m1[i]) * dx_dx;
        u_w_m1[i] += (u_w_p1[i] - u_w_m1[i]) * dx_dx;
    }
    std::copy(d_m1,d_m1+ny,this->D);
    std::copy(u_m1,u_m1+ny,this->U);
    std::copy(v_m1,v_m1+ny,this->V);
    std::copy(u_w_m1,u_w_m1+ny,this->U_w);

    delete [] d_m1;
    delete [] d_p1;
    delete [] u_m1;
    delete [] u_p1;
    delete [] v_m1;
    delete [] v_p1;
    delete [] u_w_m1;
    delete [] u_w_p1;
    //,u_m1,v_m1,u_w_m1;
    
}














