/*********************************************************************************
 * 
 *                  HydroParametrizations.cpp (Morfo70)
 * 
 *  Implements parametrizations for solving hydro field 
 * 
 *  Main methods: 
 *      - calculateDragTerm: get the related bottom shear stress quantity tau/(rho*D*U), 
 *        following Federsen formulation
 *      - calculateDragAndViscosity: Computes the corresponding paramatrizations for 
 *        get drag coefficient (cD) and the horizontal momentum difusivity term (nu)
 * 
 *   Methods and parameters description in "Hydro.h"  
 * 
 * *******************************************************************************/
#include "Hydro.h"
    

    void Hydro::calculateDragTerm(double  __restrict_arr * drag, double  __restrict_arr * cd,double __restrict_arr * urms, double  __restrict_arr * u, double  __restrict_arr * v, double  __restrict_arr * inv_d,int ny){        
        const double alpha_f2=(Parameters::Hydro::alpha_f*Parameters::Hydro::alpha_f);

        const double inv_sqrt2=1/sqrt_2;        
            for (int j=0;j<ny;j++)
            {   
                drag[j]=cd[j]*urms[j]*inv_sqrt2*sqrt(alpha_f2+2*(v[j]*v[j]+u[j]*u[j])/(urms[j]*urms[j]))*inv_d[j];
            }
            
    }
    
    void Hydro::calculateDragAndViscosity(){        
        int ny=this->grid->ny;
        int nx=this->grid->nx;
        int n=ny;
        #pragma omp parallel for shared(n)
        for (int i=0;i<nx;i++){

            double * D=this->bathy->D+i*ny;
            double * cD=this->cD+i*ny;
            (this->*dragCoefficient)(D,n,cD);
            double * Dw=this->wave->Dw+i*ny;
            double * D_r=this->wave->D_r+i*ny;
            double * H=this->wave->H+i*ny;
            double * nuD=this->nuD+i*ny;
            (this->*viscosity)(D,Dw,D_r,H,n,  nuD );
        }
        
    }

    void  Hydro::dragCoefficient_Constant( double* D, int n,double* cD)   {
            std::fill_n(cD,n,Parameters::Hydro::cD0);
    }
    
    void  Hydro::dragCoefficient_LogDepth( double* D, int n,double * cD)   {                
        const double z0=Parameters::Bottom::z0;            
        
        for (int i=0;i<n;i++){           

            double d=D[i];                        
            d=max(d,z0+0.001);
            const double cd=0.4/(log(d/z0)-1.0);                
            cD[i]=cd*cd;            
        }
    }
   
    void  Hydro::dragCoefficient_ManningStrickler( double* D, int n,double * cD)   {
        const double ka=Parameters::Hydro::ka;
        
        for (int i=0;i<n;i++)
            cD[i]=0.015*cbrt(ka/D[i]);
        
    }

    void  Hydro::viscosity_Constant(double*D, double*Dw,double*Dr, double*H, int n, double * nuD ){
        std::copy(D,D+n,nuD);        
        for (int i=0;i<n;i++)
            nuD[i]*=Parameters::Hydro::nu0;
    }
    
    void Hydro::viscosity_BreakingRollers( double*D, double*Dw, double*Dr, double*H, int n, double * nuD ){
        std::copy(D,D+n,nuD);        
        double M=Parameters::Hydro::M;
        double nu0=Parameters::Hydro::nu0;        
        double alpha_r=Parameters::Rollers::alpha_r;
        double inv_rho=1/rho;
        for (int i=0;i<n;i++)
        {

            nuD[i]*=nu0+M*D[i]*cbrt(((1-alpha_r)*Dw[i]+Dr[i])*inv_rho);         
        }
    }

    void Hydro::viscosity_BreakingDepth( double*D, double*Dw,double*Dr, double*H, int n, double * nuD ){
        std::copy(D,D+n,nuD);        
        double M=Parameters::Hydro::M;
        double nu0=Parameters::Hydro::nu0;                
        double inv_rho=1/rho;
        for (int i=0;i<n;i++)
        {
            nuD[i]*=nu0+M*D[i]*cbrt(Dw[i]*inv_rho);
        }
    }

    void Hydro::viscosity_BreakingHrms( double*D, double*Dw, double*Dr, double*H, int n, double * nuD ){
        std::copy(D,D+n,nuD);        
        double M=Parameters::Hydro::M;
        double nu0=Parameters::Hydro::nu0;                
        float inv_rho=1/rho;
        for (int i=0;i<n;i++)
        {            
            nuD[i]*=nu0+M*H[i]*cbrt(Dw[i]*inv_rho);         
            
        }
    }

