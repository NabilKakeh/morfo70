/*******************************************************************************************
  Wave Parametrizations:
     - Parallel (OpenMP) dispersion relationship with doppler (Defatult->Newton-Raphson and fast-> explicit Ondina solution)
     - Wave Breaking dissipation: Church-Thornton, Thornton-Guza and Regular waves formulations
     - Wave Bottom Friction dissipation: computes integral of bottom shear stress (Horikawa 1988) using constanta drag coeffieient or Soulsby 1997

    Methods and parameters description in "WaveRow.h" 

*******************************************************************************************/

#include "WaveRow.h"


void WaveRow::dispersion(int indexIni,int chunk,double w){
  
  const double z0=Parameters::Bottom::z0;
  const int maxIterations=10;
  const double tolerance=1E-14;
  const double  __restrict_arr *_D=this->D+indexIni;
  const double  __restrict_arr *_U_w=this->U_w+indexIni;
  double  __restrict_arr *_lK=this->lK+indexIni;
  double __restrict_arr *_Sigma=this->Sigma+indexIni;
  double __restrict_arr *_C=this->C+indexIni;
  double __restrict_arr *_Cg=this->Cg+indexIni;
  double __restrict_arr *_Urms_H=this->Urms_H+indexIni;

  for (int i=0;i<chunk;i++)
  {      
      const double d=_D[i];      
      const double inv_sqrtgh=1/sqrt(g*d);
      double alpha=_U_w[i]*inv_sqrtgh;      
      //alpha=0;
      const double mu=w*w*d/g;
      const double sqrtmu=sqrt(mu);
      
      double k=w*inv_sqrtgh;
      double kD=k*d;
      double tanhkD=tanh(kD);

      int iter=1;
      double err=9999;
    while (iter<maxIterations && err>tolerance)      
    {        
        const double mu_=(sqrtmu-alpha*kD);
        const double f=kD*tanhkD-mu_*mu_;
        const double df_dkD=tanhkD+kD*(1-tanhkD*tanhkD)+2*alpha*mu_;
        const double incrkD=f/df_dkD;
        kD+=-incrkD;
        tanhkD=tanh(kD);
        
        err=abs(f);
        iter+=1;


    }

    k=kD/d;
    const double sigma = sqrt(g * k * tanhkD);
    _Sigma[i] = sigma;
    _lK[i] = log(k);
    //C and Cg
    const double cc = sigma / k;
    _C[i] = cc;
    _Cg[i] = 0.5 * (1 + k * d * (1.0 / tanhkD - tanhkD)) * cc;    
    //Urms_H
    const double coshkd_2=1/(1-tanhkD*tanhkD);                          
    _Urms_H[i]=0.5*g*cosh(k*z0)/(cc*sqrt(coshkd_2));
  }
}


void WaveRow::dispersion_fast(int indexIni,int chunk,double w)
{
       double *D = new double[chunk];
        double *U_w = new double[chunk];
        double *lK = new double[chunk];
        //double *SinhKD = new double[chunk];
        double *Urms_H = new double[chunk];
        double *Sigma = new double[chunk];
        double *C = new double[chunk];
        double *Cg = new double[chunk];
    
        std::copy(this->D+indexIni, this->D +indexIni +chunk, D);
        std::copy(this->U_w+indexIni, this->U_w+indexIni + chunk, U_w);

    const double w2=w*w;
    double nu_D = w2 / g;
    double a = 14.979441698850072;
    double b = 1.1491762457673198;
    double c = 89.70572188644893;
    double d = -4.073909516735278;
    double e = 5.028261191874598;
    double z0=Parameters::Bottom::z0;

    for ( int i = 0; i < chunk; i++)
    {
        
        double dd = D[i];
        double u_w = U_w[i];
        double sqrt_gD = sqrt(g * dd);
        double nu = nu_D * dd;
        double nu2 = nu * nu;
        double nu3 = nu2 * nu;
        double sqrt_nu = sqrt(nu);
        double alpha = u_w / sqrt_gD;        
        double _2alpha = 2 * alpha;

        //Inital guest (Ondina 2019)
        double f = -(a * nu + b * nu2 + nu3) / (c + d * nu + e * nu2 + nu3) - 1;

        //mu0
        //double f_aux = (sqrt(-expm1(nu * f)) + _2alpha * sqrt_nu);
        double f_aux = (sqrt(1-exp(nu * f)) + _2alpha * sqrt_nu);
        double mu0 = nu / f_aux;
        //Newton-Raphson

        double mu = mu0;
        double alphamu, tanh_mu,mu_tanh_mu, G, dG_dmu, aux;
        for (int j = 0; j < 3; j++)
        {
            alphamu = alpha * mu;
            //alphamu=0;
            tanh_mu = tanh(mu);
            mu_tanh_mu = mu * tanh_mu;
            aux = sqrt_nu - alphamu;
            G = mu_tanh_mu - aux * aux;
            dG_dmu = tanh_mu + mu - tanh_mu * mu_tanh_mu + _2alpha * aux;
            mu = mu - G / dG_dmu;
        }

        double k = mu / dd;
        double sigma = sqrt(g * k * tanh_mu);
        Sigma[i] = sigma;
        lK[i] = log(k);
        //C y Cg
        double cc = sigma / k;
        C[i] = cc;
        Cg[i] = 0.5 * (1 + k * dd * (1.0 / tanh_mu - tanh_mu)) * cc;

        
        //URMS_H
         double coshkd_2=1/(1-tanh_mu*tanh_mu);                                    
         Urms_H[i]=0.5*g*cosh(k*z0)/(cc*sqrt(coshkd_2));


    }


    std::copy(lK, lK + chunk, this->lK+indexIni);
    std::copy(C, C + chunk, this->C+indexIni);
    std::copy(Cg, Cg + chunk, this->Cg+indexIni);
    std::copy(Sigma, Sigma + chunk, this->Sigma+indexIni);
    //std::copy(SinhKD, SinhKD + chunk, this->SinhKD+indexIni);
    std::copy(Urms_H, Urms_H + chunk, this->Urms_H+indexIni);
    
    delete[] D;
    delete[] lK;
    delete[] Cg;
    delete[] C;
    delete[] Sigma;
    delete[] U_w;
    //delete[] SinhKD;
    delete[] Urms_H;
   
}


void WaveRow::dispersion_parallel(const double T,const int* _indexIni, const int* _chunk,int np)
{    
   int* chunk=new int[np];
   int* indexIni=new int[np]();
   std::copy(_chunk,_chunk+np, chunk);
    std::copy(_indexIni,_indexIni+np, indexIni);    
    double w = 2 * pi / T;    
    //double w2 = w * w;
    

    #pragma omp parallel shared(indexIni,chunk,w)  
    {
        int pid=Omp::get_pid();        
        dispersion(indexIni[pid],chunk[pid],w);
    }

    delete[] chunk;delete[] indexIni;
    
}


void WaveRow::waveBreakingDissipationCT( int chunk,const double* _H,const double * _Sigma, const double *_D,const int* _nodesToSolve,int nNodesToSolve,double * _dw, double * _ddw_dh){        
    double * H=new double[chunk];    
    std::copy(_H,_H+chunk,H);
    double * Sigma=new double[chunk];
    std::copy(_Sigma,_Sigma+chunk,Sigma);    
    double * D=new double[chunk];
    std::copy(_D,_D+chunk,D);   
    int * nodesToSolve=new int[chunk];
    std::copy(_nodesToSolve,_nodesToSolve+chunk,nodesToSolve);    

    double * dw=new double[chunk];
    std::copy(_dw,_dw+chunk,dw);
    double * ddw_dh=new double[chunk];      
    double inv_gamma_b=1/Parameters::Wave::gamma_b;  // 1/gamma_b   
    double inv_4sqrtpi=0.25/sqrt_pi; // 1/(4*sqrt(pi))
    double B33=3*Parameters::Wave::B3;  // 3*B^3      
    
    for (int i=0;i<nNodesToSolve;i++)
    {           
        int nodeIndex=nodesToSolve[i];
        double h=H[nodeIndex];
        double inv_D=1/D[nodeIndex]; // 1/D
        double f1=B33*Sigma[nodeIndex]*inv_4sqrtpi*inv_D; // 3*B^3*sigma/( 4*sqrt(pi)*D)
        double dcc_dh=inv_D*inv_gamma_b;  // 1/( gamma_b*D)
        double cc=h*dcc_dh;
        double cc_2=cc*cc;
        double inv_1pcc_2=1/(1+cc_2);
        double pow_inv_1pcc_2_2d5=inv_1pcc_2*inv_1pcc_2*sqrt(inv_1pcc_2);
        double pow_inv_1pcc_2_3d5=inv_1pcc_2*pow_inv_1pcc_2_2d5;
        double f2=1-pow_inv_1pcc_2_2d5;
        double df2_h=-2.5*pow_inv_1pcc_2_3d5*2*cc_2;        
        double tanh8ccm1=tanh(8*(cc-1));        

        double f3=1+tanh8ccm1;
        double tanh8ccm1_2=tanh8ccm1*tanh8ccm1;
        double df3_h=(1-tanh8ccm1_2)*8.*dcc_dh;
        dw[nodeIndex]=h*f1*f2*f3;
        ddw_dh[nodeIndex]=f1*(f2*f3+h*(df2_h*f3+df3_h*f2));  
    }
     std::copy(dw,dw+chunk,_dw);
    std::copy(ddw_dh,ddw_dh+chunk,_ddw_dh);
    delete[]D;
    delete[] Sigma;
     delete[] H;
     delete[] nodesToSolve;
    
     delete[] ddw_dh;delete[] dw;
    
}


void WaveRow::waveBreakingDissipationTG( int chunk,const double* _H,const double * _Sigma, const double *_D,const int* _nodesToSolve,int nNodesToSolve,double * _dw, double * _ddw_dh){        
    double * H=new double[chunk];    
    std::copy(_H,_H+chunk,H);
    double * Sigma=new double[chunk];
    std::copy(_Sigma,_Sigma+chunk,Sigma);    
    double * D=new double[chunk];
    std::copy(_D,_D+chunk,D);   
    int * nodesToSolve=new int[chunk];
    std::copy(_nodesToSolve,_nodesToSolve+chunk,nodesToSolve);    

    double * dw=new double[chunk]; //Dw/E
    std::copy(_dw,_dw+chunk,dw);
    double * ddw_dh=new double[chunk];  //d(Dw/E)_dH

    
    
    const double inv_gamma_b=1/Parameters::Wave::gamma_b;   // 1/gamma_b 
    const double kappa_sigma=3*Parameters::Wave::B3/(4*sqrt_pi*Parameters::Wave::gamma_b*Parameters::Wave::gamma_b);  //B^3/(4*sqrt(pi)*gamma_b^2)
    
    for (int i=0;i<nNodesToSolve;i++)
    {           
        const int nodeIndex=nodesToSolve[i];        
        const double h=H[nodeIndex];
        const double kappa=kappa_sigma*Sigma[nodeIndex];
        const double inv_D=1/D[nodeIndex];
        const double gamma=h*inv_D;
        const double g2=gamma*gamma;
        const double g3=g2*gamma;
        const double dcc_dh=inv_D*inv_gamma_b;
        const double cc=h*dcc_dh;
        const double ff=(1+cc*cc);
        const double ff1=pow(ff,-2.5);                
        const double ff2=(1-ff1);
        const double ff3=ff1/ff;        
        dw[nodeIndex]=kappa*g3*ff2;                                                
        ddw_dh[nodeIndex]=kappa*g2*(3*ff2*inv_D+2.5*gamma*ff3*2*cc*dcc_dh);      
    }
     std::copy(dw,dw+chunk,_dw);
    std::copy(ddw_dh,ddw_dh+chunk,_ddw_dh);
    delete[]D;
    delete[] Sigma;
     delete[] H;
     delete[] nodesToSolve;
    
     delete[] ddw_dh;delete[] dw;
    
}


void WaveRow::waveBreakingDissipationRegular( int chunk,const double* _H,const double * _Sigma, const double *_D,const int* _nodesToSolve,int nNodesToSolve,double * _dw, double * _ddw_dh){        
    double * H=new double[chunk];    
    std::copy(_H,_H+chunk,H);
    double * Sigma=new double[chunk];
    std::copy(_Sigma,_Sigma+chunk,Sigma);    
    double * D=new double[chunk];
    std::copy(_D,_D+chunk,D);   
    int * nodesToSolve=new int[chunk];
    std::copy(_nodesToSolve,_nodesToSolve+chunk,nodesToSolve);    

    double * dw=new double[chunk];  //Dw/E
    std::copy(_dw,_dw+chunk,dw);
    double * ddw_dh=new double[chunk]; //d(Dw/E)_dH

    
    
    const double inv_gamma_b=1/Parameters::Wave::gamma_b;    //  1/gamma_b
    const double B3_pi=Parameters::Wave::B3/pi; // B^3/pi 
    const double m_mon=Parameters::Wave::m_mon;
    const double n_mon=Parameters::Wave::n_mon;

    
    for (int i=0;i<nNodesToSolve;i++)
    {   const int nodeIndex=nodesToSolve[i];        
        const double h=H[nodeIndex];        
        const double inv_D=1/D[nodeIndex];        
        const double dcc_dh=inv_D*inv_gamma_b;  //1/ (gamma_b*D)
        const double cc=h*dcc_dh;    //H/(gamma_b*D)

        const double cc_to_n=pow(cc,n_mon);
        const double cc_to_m=pow(cc,m_mon);
        const double cc_to_mm1=(max(0.0,cc_to_m/cc));
        const double cc_to_nm1=(max(0.0,cc_to_n/cc));
        //const double cc_to_mm1=pow(cc,m_mon-1);
        //const double cc_to_nm1=pow(cc,n_mon-1);
        const double f1=B3_pi*inv_D*Sigma[nodeIndex];        
        const double f2=cc_to_m;        
        const double f3=-expm1(-cc_to_n);
        
        const double df2_dh=m_mon*(cc_to_mm1)*dcc_dh;
        const double df3_dh=(f3-1)*(-n_mon*(cc_to_nm1))*dcc_dh;        
        dw[nodeIndex]=f1*h*f2*f3;
        ddw_dh[nodeIndex]=f1*(f2*f3+h*df2_dh*f3+h*df3_dh*f2);                 
    }
     std::copy(dw,dw+chunk,_dw);
    std::copy(ddw_dh,ddw_dh+chunk,_ddw_dh);
    delete[]D;
    delete[] Sigma;
     delete[] H;
     delete[] nodesToSolve;
    
     delete[] ddw_dh;delete[] dw;
    
}


void WaveRow::bottomFrictionConstant(int chunk, const double *_H,const double _sigma,const double * _Urms_H,const int* _nodesToSolve,int nNodesToSolve,double * _df, double * _ddf_dh){
    
    
    double * H=new double[chunk];
    std::copy(_H,_H+chunk,H);
    int * nodesToSolve=new int[chunk];
     std::copy(_nodesToSolve,_nodesToSolve+chunk,nodesToSolve);    
    double * Urms_H=new double[chunk];
    std::copy(_Urms_H,_Urms_H+chunk,Urms_H);   
    
    
     double * df=new double[chunk];   //Df/E
     double * ddf_dh=new double[chunk];    //d(Df/E)_dH
     std::copy(_df,_df+chunk,df);
     std::copy(_ddf_dh,_ddf_dh+chunk,ddf_dh);  
    double _64fw_6pig=21.333333333*Parameters::Wave::fw0*inv_2pi/g;    //Cte -> 64*fw/(6*pi*g)
    
    
    for (int i=0;i<nNodesToSolve;i++)
    {        
        double urms_h_cube=Urms_H[i];
        urms_h_cube=urms_h_cube*urms_h_cube*urms_h_cube; //(Urms/H)^3
        
        ddf_dh[i]=_64fw_6pig*urms_h_cube;
        df[i]=ddf_dh[i]*H[i];

    }
    std::copy(df,df+chunk,_df);
    std::copy(ddf_dh,ddf_dh+chunk,_ddf_dh);
    //delete[] SinhKD;
    delete[] Urms_H;
    delete[] H;
    delete[] nodesToSolve;
    delete[] ddf_dh;delete[] df;

}


void WaveRow::bottomFrictionSoulsby97(int chunk, const double *_H,const double _sigma,const double * _Urms_H,const int* _nodesToSolve,int nNodesToSolve,double * _df, double * _ddf_dh){
    
    
    double * H=new double[chunk];
    std::copy(_H,_H+chunk,H);
    int * nodesToSolve=new int[chunk];
     std::copy(_nodesToSolve,_nodesToSolve+chunk,nodesToSolve);    
    double * Urms_H=new double[chunk];
    std::copy(_Urms_H,_Urms_H+chunk,Urms_H);   
    
    
    double * df=new double[chunk];    //Df/E
    double * ddf_dh=new double[chunk];    //d(Df/E)_dH
    double z0=Parameters::Bottom::z0;    //Bed apparent Rougnnes
    double _64_6pig=21.3333333*inv_2pi/g; // Cte->64/(6*pi*g)
    
    std::copy(_df,_df+chunk,df);
    std::copy(_ddf_dh,_ddf_dh+chunk,ddf_dh);      
    
    
    
    for (int i=0;i<nNodesToSolve;i++)
    {    

        //Computes wave skin friction  coefficient (fw) according Souslby 1997
        double urms_h=Urms_H[i];
        double h=H[i];
        double urms=urms_h*h;
        double fw=1.39*sqrt(z0*_sigma/urms); //fast fomula x^0.5 ~ x^0.52
        //double fw=1.39*pow(z0*_sigma/urms,0.52); //real formula, slower
        double urms_h_cube=urms_h*urms_h*urms_h;
        
        ddf_dh[i]=_64_6pig*fw*urms_h_cube;
        df[i]=ddf_dh[i]*h;

    }
    std::copy(df,df+chunk,_df);
    std::copy(ddf_dh,ddf_dh+chunk,_ddf_dh);
    
    delete[] Urms_H;
    delete[] H;
    delete[] nodesToSolve;
    delete[] ddf_dh;delete[] df;

}




