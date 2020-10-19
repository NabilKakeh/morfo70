/*********************************************************************************
 * 
 *                  HydroQx.cpp (Morfo70)
 * 
 *  Implementation of Methods to solve x momentum equation of NLSW.
 *  It is solved in tws steps: 
 *      1) Predictor d_dy terms ar ommied and equation is integrated chrosshore (x)
 *      2) Projector, predictor solution is corrected integratind d_dy terms 
 * 
 *  Main method: solveQx
 * 
 *   Methods and parameters description in "Hydro.h"  
 * 
 * *******************************************************************************/
#include "Hydro.h"


void Hydro::solveQx()
{
    int nx=this->grid->nx;
    int ny=this->grid->ny;

    this->C2X(this->Ui_j,this->V_c,nx,ny);     //Get V in x-cells

    //////////////////////////////////////////////
    //Prediction Step A*Qxx=f, A=[low, d, up]
    //////////////////////////////////////////////
    //Set prediction solver 
    this->getD_x0();
    this->getD_x(0);
    this->setWaterNodes_x(0);
    this->calculateDragTerm(this->DragTerm, this->cD_x,this->Urms_x,this->U+ny, this->Ui_j,this->inv_Di+ny,ny);
    this->setSolverQxx_f(0);
    this->setSolverQxx_d(0);
    this->setSolverQxx0_low(); //Caso especial que no esta en tierra seca nunca    
    #pragma omp parallel for
    for (int i=1;i<nx-1;i++)
    {     
        this->getD_x(i);          
        this->setWaterNodes_x(i);
        this->setSolverQxx_up(i-1);
        this->calculateDragTerm(this->DragTerm+i*ny, this->cD_x+i*ny,this->Urms_x+i*ny,this->U+(i+1)*ny, this->Ui_j+i*ny,this->inv_Di+(i+1)*ny,ny);
        this->setSolverQxx_f(i);
        this->setSolverQxx_d(i);           
        this->setSolverQxx_low(i);         
    }       
    this->setSolverQxx_up(nx-2);   //Debe ser igual a cero
    this->setSolverQxx_BC();    //Set boundary conditions in matrix solver

     //Solve Qxx
     Math::TridiagonalSolver(this->Qx+ny,this->v_tmp,this->d,this->low,this->up,this->f,nx-1,ny);
    //////////////////////////////////////////////
    //END Prediction Step A*Qxx=f, A=[low, d, up]
    //////////////////////////////////////////////

    //////////////////////////////////////////////
    //Projector Step A*Qxy=f, A=[low, d, up]
    //////////////////////////////////////////////
    #pragma omp parallel for
    for (int i=0;i<nx-1;i++)
    {
         //Set solver
         this->setSolverQxy_f(i);
         this->setSolverQxy_d(i);
         this->setSolverQxy_up(i);
         this->setSolverQxy_low(i);
         //Solve periodic tridiagonal system
         Math::TridiagonalSolverPeriodic(this->d+i*ny,this->low+i*ny,this->up+i*ny,this->f+i*ny,ny,this->Qx+(i+1)*ny);        
         
    }
    //////////////////////////////////////////////
    //END Projector Step A*Qxy=f, A=[low, d, up]
    //////////////////////////////////////////////    

    //Set offshore boundary conditio to solution
     this->setQx_BC();


    //Get Zs
    this->difC2X(this->dQx_dx,this->Qx,this->grid->dx_x,nx+1,ny);        
    #pragma omp parallel for
    for (int i=0;i<nx;i++)    
        this->solveZs_guess(i);

    //Get U
     this->getU0();
      #pragma omp parallel for
      for (int i=0;i<nx-1;i++)    
          this->getU(i);
    
    //Get U_c (in center nodes)
    MeshOperations::X2C(this->U_c,this->U,this->grid->x2c_d,this->grid->x2c_up,nx,ny);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          Qx Prediction Solver Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Hydro::setSolverQxx_f(int i){
     
        const int ny=this->grid->ny;    
        const double dt=Parameters::Hydro::dt;    
        const double inv_dt=1/dt;
        const double inv_dy=1/this->grid->dy;
        const double inv_dx=1/this->grid->dx[i];        
        double  __restrict_arr *zs=this->Zs+i*ny;
        double  __restrict_arr *dqy_dy=this->dQy_dy+i*ny;
        double  __restrict_arr *v=this->V+(i+1)*ny;
        double  __restrict_arr *qx=this->Qx+(i+1)*ny;
        double  __restrict_arr *forcing_x=this->Forcing_x+i*ny;
        double  __restrict_arr *d_x=this->Di+(i+1)*ny;        
        double  __restrict_arr *nud_x=this->nuD_x+i*ny;
        double  __restrict_arr *dnud_dy=this->dnuDx_dy+i*ny;
        bool  __restrict_arr *isWater=this->IsWater+i*ny;
        double   *_f=this->f+i*ny;

        for (int j=1;j<ny;j++)
        {            
            const double dzs_dx=inv_dx*(zs[ny+j]-zs[j]);                             
            const double d2qy_dxy=inv_dx*(dqy_dy[ny+j]-dqy_dy[j]);                    
            const double dV_dx=0.5*inv_dx*((v[ny+j]+v[ny+j-1])-(v[j]+v[j-1]));        
            const double d2V_dxy=inv_dx*inv_dy*((v[ny+j]-v[ny+j-1])-(v[j]-v[j-1]));                            
            _f[j]=isWater[j]*(forcing_x[j]+qx[j]*inv_dt+g*d_x[j]*(dt*d2qy_dxy -dzs_dx)+dnud_dy[j]*dV_dx+nud_x[j]*d2V_dxy);               
    
        }
        //BC
         const double dzs_dx=inv_dx*(zs[ny]-zs[0]);        
         const double d2qy_dxy=inv_dx*(dqy_dy[ny]-dqy_dy[0]);        
         const double dV_dx=0.5*inv_dx*((v[ny]+v[2*ny-1])-(v[0]+v[ny-1]));        
         const double d2V_dxy=inv_dx*inv_dy*((v[ny]-v[2*ny-1])-(v[0]-v[ny-1]));        
        
         _f[0]=isWater[0]*(forcing_x[0]+qx[0]*inv_dt+g*d_x[0]*(dt*d2qy_dxy -dzs_dx)+dnud_dy[0]*dV_dx+nud_x[0]*d2V_dxy);                        
}

void Hydro::setSolverQxx_d(int column){
    
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;    
    const double inv_dt=1/dt;    
    const int i=column;
    const double difx=this->grid->difx_x_d[i];
    const double dif2x=this->grid->dif2x_x_d[i];     
    double  __restrict_arr *d_x=this->Di+(i+1)*ny; //????
    double  __restrict_arr *inv_d_x=this->inv_Di+(i+1)*ny; //????    
    double  __restrict_arr *nud_x=this->nuD_x+i*ny;
    double  __restrict_arr *dnud_dx=this->dnuD_dx+i*ny;
    double  __restrict_arr *dragTerm=this->DragTerm+i*ny;    
    double  __restrict_arr *u=this->U+(i+1)*ny;    
    //double const __restrict_arr *dnud_dy=this->dnuDx_dy+i*ny;
    bool  __restrict_arr *isWater=this->IsWater+i*ny;
    double  __restrict_arr *_d=this->d+i*ny;

        for (int j=0;j<ny;j++)        
        {
            bool isl=!isWater[j]; //is land
            //const double inv_d_x=1/d_x[j];
            _d[j]=(inv_dt+dragTerm[j]+difx*(u[j]-2*dnud_dx[j]*inv_d_x[j])-dif2x*(g*d_x[j]*dt+2*nud_x[j]*inv_d_x[j]))*isWater[j]+isl;

        }
    
    
}

void Hydro::setSolverQxx_up(int column){    
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;    
    const int i=column;
    const double difx=this->grid->difx_x_up[i];
    const double dif2x=this->grid->dif2x_x_up[i];
    double  __restrict_arr *uup=this->U+(i+2)*ny;    
    double  __restrict_arr *d_x=this->Di+(i+1)*ny; //????
    double  __restrict_arr *inv_d_xup=this->inv_Di+(i+2)*ny; //????    
    double  __restrict_arr *nud_x=this->nuD_x+i*ny;
    double  __restrict_arr *dnud_dx=this->dnuD_dx+i*ny;
    bool    __restrict_arr *isWater=this->IsWater+i*ny;
    double  __restrict_arr *_up=this->up+i*ny;
    

    
        
        for (int j=0;j<ny;j++)        
        {    
             double _isWater=isWater[j]&&isWater[j+ny];         
            _up[j]=(difx*(uup[j]-2*dnud_dx[j]*inv_d_xup[j])-dif2x*(g*d_x[j]*dt+2*nud_x[j]*inv_d_xup[j]))*_isWater; 
        }
 
 
}

void Hydro::setSolverQxx_low(int column){    
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;    
    const int i=column;
    const double difx=this->grid->difx_x_low[i];
    const double dif2x=this->grid->dif2x_x_low[i];
    double __restrict_arr *ulow=this->U+i*ny;    
    double __restrict_arr *d_x=this->Di+(i+1)*ny; //????
    double __restrict_arr *inv_d_xlow=this->inv_Di+i*ny; //????    
    double __restrict_arr *nud_x=this->nuD_x+i*ny;
    double __restrict_arr *dnud_dx=this->dnuD_dx+i*ny;
    bool   __restrict_arr *isWater=this->IsWater+i*ny;
    double  __restrict_arr *_low=this->low+i*ny;
        
    for (int j=0;j<ny;j++)        
    {            
        double _isWater=isWater[j]&&isWater[j-ny];         
        _low[j]=(difx*(ulow[j]-2*dnud_dx[j]*inv_d_xlow[j])-dif2x*(g*d_x[j]*dt+2*nud_x[j]*inv_d_xlow[j]))*_isWater;         

    }
 
 
}

void Hydro::setSolverQxx0_low(){    
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;    
    const int i=0;
    const double difx=this->grid->difx_x_low[i];
    const double dif2x=this->grid->dif2x_x_low[i];
    double __restrict_arr *ulow=this->U+i*ny;    
    double __restrict_arr *d_x=this->Di+(i+1)*ny; //????
    double __restrict_arr *inv_d_xlow=this->inv_Di+i*ny; //????    
    double __restrict_arr *nud_x=this->nuD_x+i*ny;
    double __restrict_arr *dnud_dx=this->dnuD_dx+i*ny;    
    double  __restrict_arr *_low=this->low+i*ny;
        
    for (int j=0;j<ny;j++)        
        _low[j]=(difx*(ulow[j]-2*dnud_dx[j]*inv_d_xlow[j])-dif2x*(g*d_x[j]*dt+2*nud_x[j]*inv_d_xlow[j]));         
}

void Hydro::setWaterNodes_x(int column)
{
    
    const int ny=this->grid->ny;
    const int i=column;
    const double Dmin=Parameters::Hydro::Dmin+1E-4;
    double   __restrict_arr * d=this->bathy->D+i*ny;    
    double   __restrict_arr * zs=this->Zs+i*ny;    
    double   __restrict_arr * dqx_dx=this->dQx_dx+i*ny;
    
    bool   __restrict_arr * isWater=this->IsWater+i*ny;
    
        for (int j=0;j<ny;j++)
        {
            bool dup_gt_min=d[j+ny]>Dmin;
            bool d_gt_min=d[j]>Dmin;
            bool zsup_gt_zs=zs[j+ny]>zs[j];
            bool zs_gt_zsup=zs[j]>zs[j+ny];            
            bool sink=dqx_dx[j]>0;
            bool sinkup=dqx_dx[j+ny]>0;

            
            isWater[j]=d_gt_min*dup_gt_min+(dup_gt_min*(!d_gt_min)*zsup_gt_zs+d_gt_min*(!dup_gt_min)*zs_gt_zsup)*!(sink+sinkup);           
        }
        
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          Qx Proyector Solver Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Hydro::setSolverQxy_f(int column){
     
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column;    
    double __restrict_arr  * dragTerm=this->DragTerm+ i*ny;
    double __restrict_arr  * qx=this->Qx+(i+1)*ny;    
    bool __restrict_arr  * isWater=this->IsWater+i*ny;
    double __restrict_arr  * _f=this->f+i*ny;            

        for (int j=0;j<ny;j++)      
            _f[j]=isWater[j]*qx[j]*(1+dt*dragTerm[j]);
            //_f[j]=isWater[j]*(qx[j]/dt+dragTerm[j]*qx[j]);
    
}

void Hydro::setSolverQxy_d(int column){
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column;    
    const double _2inv_dy2=2/(this->grid->dy*this->grid->dy);
    double __restrict_arr * dragTerm=this->DragTerm+ i*ny;
    double __restrict_arr  * nud_x=this->nuD_x+i*ny;    
    double __restrict_arr  * inv_d_x=this->inv_Di+(i+1)*ny;    
    bool __restrict_arr  * isWater=this->IsWater+i*ny;    
    double __restrict_arr  * _d=this->d+i*ny;            
        
        for (int j=0;j<ny;j++)      
        {   
            bool isl=!isWater[j];
            _d[j]=isWater[j]*((1+dt*dragTerm[j])+dt*inv_d_x[j]*nud_x[j]*_2inv_dy2)+isl;
            //_d[j]=isWater[j]*(1/dt+dragTerm[j]+inv_d_x[j]*nud_x[j]*_2inv_dy2)+isl;
        }
}

void Hydro::setSolverQxy_up(int column){         
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column;    
    const double inv_2dy=0.5/this->grid->dy;
    const double inv_dy2=4*inv_2dy*inv_2dy;

    double __restrict_arr  * nud_x=this->nuD_x+i*ny;    
    double __restrict_arr  * dnudx_dy=this->dnuDx_dy+i*ny;    
    double __restrict_arr  * inv_d_x=this->inv_Di+(i+1)*ny;    
    double  __restrict_arr  *v=this->Ui_j+i*ny;
    bool __restrict_arr  * isWater=this->IsWater+i*ny;  
    double __restrict_arr  * _up=this->up+i*ny;                  
        
    for (int j=0;j<ny-1;j++)      
    {            
        double _isWater=isWater[j]&&isWater[j+1];       
        _up[j]=dt*(v[j+1]*inv_2dy-inv_d_x[j+1]*(nud_x[j]*inv_dy2+dnudx_dy[j]*inv_2dy))*_isWater;
        //_up[j]=(v[j+1]*inv_2dy-inv_d_x[j+1]*(nud_x[j]*inv_dy2+dnudx_dy[j]*inv_2dy))*isWater[j];
    }
    
    //BC periodicas
    double _isWater=isWater[ny-1]&&isWater[0];
    _up[ny-1]=dt*(v[0]*inv_2dy-inv_d_x[0]*(nud_x[ny-1]*inv_dy2+dnudx_dy[ny-1]*inv_2dy))*_isWater;    
    //_up[ny-1]=(v[0]*inv_2dy-inv_d_x[0]*(nud_x[ny-1]*inv_dy2+dnudx_dy[ny-1]*inv_2dy))*isWater[ny-1];    
}

void Hydro::setSolverQxy_low(int column){         
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column;    
    const double inv_2dy=0.5/this->grid->dy;
    const double inv_dy2=4*inv_2dy*inv_2dy;

    double __restrict_arr  * nud_x=this->nuD_x+i*ny;    
    double __restrict_arr  * dnudx_dy=this->dnuDx_dy+i*ny;    
    double __restrict_arr  * inv_d_x=this->inv_Di+(i+1)*ny;    
    double  __restrict_arr *v=this->Ui_j+i*ny;
    bool __restrict_arr  * isWater=this->IsWater+i*ny;
    double __restrict_arr  * _low=this->low+i*ny;              
        
        for (int j=1;j<ny;j++)      
        {           
            double _isWater=isWater[j]&&isWater[j-1];         
             _low[j]=dt*(-v[j-1]*inv_2dy-inv_d_x[j-1]*(nud_x[j]*inv_dy2-dnudx_dy[j]*inv_2dy))*_isWater;
            //_low[j]=(-v[j-1]*inv_2dy-inv_d_x[j-1]*(nud_x[j]*inv_dy2-dnudx_dy[j]*inv_2dy))*isWater[j];
        }
        
        //BC periodicas
        double _isWater=isWater[0]&&isWater[ny-1];         
        _low[0]=dt*(-v[ny-1]*inv_2dy-inv_d_x[ny-1]*(nud_x[0]*inv_dy2-dnudx_dy[0]*inv_2dy))*_isWater;
        //_low[0]=(-v[ny-1]*inv_2dy-inv_d_x[ny-1]*(nud_x[0]*inv_dy2-dnudx_dy[0]*inv_2dy))*isWater[0];
}




/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          Pre-Solver and Post-Solver Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Hydro::getForcing_x(){
    int nx=this->grid->nx;    
    int ny=this->grid->ny;
    double inv_rho=1/rho;
    double inv_4dy=0.25/this->grid->dy;
    inv_4dy*=inv_rho;
    //double * dx=this->grid->dx;

     double * forcing_x;
     double* Sxy=this->wave->Sxy;
     double* Sxy_up=this->wave->Sxy+ny;
     double* Sxx=this->wave->Sxx;
     double* Sxx_up=this->wave->Sxx+ny;

    for (int i=0;i<nx-1;i++){
        double inv_rhodx=inv_rho/this->grid->dx[i];        
        forcing_x=this->Forcing_x+i*ny;
        for (int j=1;j<ny-1;j++)    
            forcing_x[j]=-(inv_4dy*(Sxy[j+1]-Sxy[j-1]+Sxy_up[j+1]-Sxy_up[j-1])+(Sxx_up[j]-Sxx[j])*inv_rhodx);
            //Forcing_x[i*ny+j]=-(inv_4dy*(wave->Sxy[i*nx+j+1]-wave->Sxy[i*nx+j-1]+wave->Sxy[(i+1)*nx +j+1]-wave->Sxy[(i+1)*nx +j-1])+(wave->Sxx[(i+1)*nx +j]-wave->Sxx[i*nx+j])*inv_rhodx);
           
         forcing_x[0]=-(inv_4dy*(Sxy[1]-Sxy[ny-1]+Sxy_up[1]-Sxy_up[ny-1])+(Sxx_up[0]-Sxx[0])*inv_rhodx);
         forcing_x[ny-1]=-(inv_4dy*(Sxy[0]-Sxy[ny-2]+Sxy_up[0]-Sxy_up[ny-2])+(Sxx_up[ny-1]-Sxx[ny-1])*inv_rhodx);       
        //Forcing_x[i*nx]=-(inv_4dy*(wave->Sxy[i*nx+1]-wave->Sxy[i*nx+ny-1]+wave->Sxy[(i+1)*nx+1]-wave->Sxy[(i+1)*nx+ny-1])+(wave->Sxx[(i+1)*nx]-wave->Sxx[i*nx])*inv_rhodx);
        //Forcing_x[i*nx+ny-1]=-(inv_4dy*(wave->Sxy[i*nx]-wave->Sxy[i*nx+ny-2]+wave->Sxy[(i+1)*nx]-wave->Sxy[(i+1)*nx+ny-2])+(wave->Sxx[(i+1)*nx+ny-1]-wave->Sxx[i*nx+ny-1])*inv_rhodx);       
         Sxy+=ny;
         Sxy_up+=ny;
         Sxx+=ny;
         Sxx_up+=ny;
    }
    
}

void Hydro::getD_x0()
{    
    const int ny=this->grid->ny;
    double   __restrict_arr * d=this->bathy->D;    
    double   __restrict_arr * d_x=this->Di;
    double   __restrict_arr * inv_d_x=this->inv_Di;
    //double* inv_d=this->inv_D;
    for (int j=0;j<ny;j++)    {
        d_x[j]=1.5*d[j]-0.5*d[j+ny];    
        inv_d_x[j]=1.0/d_x[j];
    }
}

void Hydro::getD_x(int column)
{    
    const int ny=this->grid->ny;
    const int i=column;
    double   __restrict_arr * d=this->bathy->D+i*ny;    
    double   __restrict_arr * d_x=this->Di+(i+1)*ny;
    double   __restrict_arr * inv_d_x=this->inv_Di+(i+1)*ny;
  
    //inv_d+=ny;
    for (int j=0;j<ny;j++)
    {
        d_x[j]=0.5*(d[j]+d[j+ny]);        
        inv_d_x[j]=1/d_x[j];
    }
    
}

void Hydro::getU0()
{
    const int ny=this->grid->ny;
    double   __restrict_arr * d=this->guess_D;    
    double   __restrict_arr * qx=this->Qx;        
    double   __restrict_arr * u=this->U;        
    
    for (int j=0;j<ny;j++)    {
        const double d_x0=1.5*d[j]-0.5*d[j+ny];    
        u[j]=qx[j]/d_x0;        
    }    
}

void Hydro::getU(int column)
{

    const int ny=this->grid->ny;
    const int i=column;
    //const double xd=this->grid->x2c_d[i];
    //double xup=this->grid->x2c_up[i];
    double   __restrict_arr * d=this->guess_D+i*ny;    
    double   __restrict_arr * qx=this->Qx+(i+1)*ny;        
    double   __restrict_arr * u=this->U+(i+1)*ny;        
    //double   __restrict_arr * uc=this->U_c+i*ny;        
    
    for (int j=0;j<ny;j++)    {
        const double d_x=0.5*(d[j]+d[j+ny]);    
        u[j]=qx[j]/d_x;  
        //uc[j]=xd*u[j-ny] +xup*u[j];      
    }    

}



