/*********************************************************************************
 * 
 *                  HydroQy.cpp (Morfo70)
 * 
 *  Implementation of Methods to solve y momentum equation of NLSW.
 *  It is solved in tws steps: 
 *      1) Predictor d_dx terms ar ommited and equation is integrated alongshore (y)
 *      2) Projector, predictor solution is corrected integratind d_dx terms along x
 * 
 *  Main method: solveQy
 * 
 *   Methods and parameters description in "Hydro.h"  
 * 
 * *******************************************************************************/
#include "Hydro.h"

void Hydro::solveQy()
{
    int nx=this->grid->nx;
    int ny=this->grid->ny; 
    

    //////////////////////////////////////////////
    //Prediction Step A*Qyy=f, A=[low, d, up]
    //////////////////////////////////////////////    
    this->getD_y0();
    #pragma omp parallel for
    for (int i=0;i<nx-1;i++)
    { 
        //Set prediction solver 
        this->setWaterNodes_y(i);      
        this->getD_y(i);
        this->setSolverQyy_f(i);
        this->C2Y(this->Ui_j+i*ny,this->U_c+i*ny,ny); //U_y
        this->calculateDragTerm(DragTerm+i*ny,cD_y+i*ny,Urms_y+i*ny,Ui_j+i*ny, V+(i+1)*ny,inv_Di+(i+1)*ny,ny);
        this->setSolverQyy_d(i);
        this->setSolverQyy_low(i);
        this->setSolverQyy_up(i);
        //Solve System
         Math::TridiagonalSolverPeriodic(this->d+i*ny,this->low+i*ny,this->up+i*ny,this->f+i*ny,ny,this->Qy+(i+1)*ny);
    }
    this->getD_y(nx-1);
    //////////////////////////////////////////////
    // END Prediction Step A*Qyy=f, A=[low, d, up]
    //////////////////////////////////////////////    

    //////////////////////////////////////////////
    // Projector Step A*Qyx=f, A=[low, d, up]
    //////////////////////////////////////////////       
    //Set projector solver    
    //First i=0 outside the loop, in order to set low(i-1)=0
    this->setSolverQyx_f(0);
    this->setSolverQyx_d(0);
    this->setSolverQyx0_low();
    this->setSolverQyx_up(0);
    #pragma omp parallel for
     for (int i=1;i<nx-1;i++)
     {         
         this->setSolverQyx_f(i);
         this->setSolverQyx_d(i);
         this->setSolverQyx_low(i);
         this->setSolverQyx_up(i);
     }
    this->setSolverQyx_BC();
    
    //Solve System
    Math::TridiagonalSolver(this->Qy+ny,this->v_tmp,this->d,this->low,this->up,this->f,nx-1,ny);
    //////////////////////////////////////////////
    // END Projector Step 
    //////////////////////////////////////////////       
    
    //Apply Onshore and Offshore Boundary condition to solution
    this->setQy_BC();

    //Zs ans V    
    #pragma omp parallel for
    for (int i=0;i<nx;i++)
    {        
        this->difY2C(this->dQy_dy+i*ny,this->Qy+(i+1)*ny,this->grid->dy,ny); 
        this->solveZs(i);
        this->getV(i);
        this->Y2C(this->V_c+i*ny,this->V+(i+1)*ny,ny);
    }
    this->getV0();
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          Qy Predictor Solver Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Hydro::setWaterNodes_y(int i)
{
    const int ny=this->grid->ny;
    const double Dmin=Parameters::Hydro::Dmin+1E-4;
    double   __restrict_arr * d=this->guess_D+i*ny;    
    double   __restrict_arr * zs=this->guess_Zs+i*ny;    
    double   __restrict_arr * dqy_dy=this->dQy_dy+i*ny;
    bool  __restrict_arr * isWater=this->IsWater+i*ny;
    
        for (int j=0;j<ny-1;j++)
        {
            bool dup_gt_min=d[j+1]>Dmin;
            bool d_gt_min=d[j]>Dmin;
            bool zsup_gt_zs=zs[j+1]>zs[j];
            bool zs_gt_zsup=zs[j]>zs[j+1];            
            bool sink=dqy_dy[j]>0;
            bool sinkup=dqy_dy[j+1]>0;            
            isWater[j]=d_gt_min*dup_gt_min+(dup_gt_min*(!d_gt_min)*zsup_gt_zs+d_gt_min*(!dup_gt_min)*zs_gt_zsup)*!(sink+sinkup);
            //isWater[j]=d_gt_min*dup_gt_min;
        }
        bool dup_gt_min=d[0]>Dmin;
        bool d_gt_min=d[ny-1]>Dmin;
        bool zsup_gt_zs=zs[0]>zs[ny-1];
        bool zs_gt_zsup=zs[ny-1]>zs[0];            
        bool sink=dqy_dy[ny-1]>0;
        bool sinkup=dqy_dy[0]>0;            
        isWater[ny-1]=d_gt_min*dup_gt_min+(dup_gt_min*(!d_gt_min)*zsup_gt_zs+d_gt_min*(!dup_gt_min)*zs_gt_zsup)*!(sink+sinkup);        
        //isWater[ny-1]=d_gt_min*dup_gt_min;
    
   
}

void Hydro::setSolverQyy_f(int i){
    const int ny=this->grid->ny;    
    const double dt=Parameters::Hydro::dt;    
    const double inv_dt=1/dt;    
    const double inv_dy=1/this->grid->dy;
    const double inv_dx=1/this->grid->dx_x[i];        
    double const x2c_d=this->grid->x2c_d[i];
    double const x2c_up=this->grid->x2c_up[i];
    double  __restrict_arr *zs=this->Zs+i*ny;
    double  __restrict_arr *dqx_dx=this->dQx_dx+i*ny;
    double  __restrict_arr *forcing_y=this->Forcing_y+i*ny;    
    double  __restrict_arr *qy=this->Qy+(i+1)*ny;    
    double  __restrict_arr * d_y=this->Di+(i+1)*ny;     
    double  __restrict_arr * u=this->U+i*ny;  
    double  __restrict_arr * nud_y=this->nuD_y+i*ny;
    double  __restrict_arr * dnudy_dx=this->dnuDy_dx+i*ny;;      
    bool  __restrict_arr * isWater=this->IsWater;  
    
    double __restrict_arr* _f=this->f+i*ny;  
    for (int j=0;j<ny-1;j++)
    {   const double dzs_dy=inv_dy*(zs[j+1]-zs[j]);  
        const double d2qx_dxy=inv_dy*(dqx_dx[j+1]-dqx_dx[j]);                                   
        const double du_dy=inv_dy*(x2c_up*(u[ny+j+1]-u[ny+j])+x2c_d*(u[j+1]-u[j]));        
        const double d2u_dxy=inv_dx*inv_dy*((u[ny+j+1]-u[ny+j])-(u[j+1]-u[j]));   
        _f[j]=isWater[j]*(forcing_y[j]+qy[j]*inv_dt+g*d_y[j]*(dt*d2qx_dxy -dzs_dy)+dnudy_dx[j]*du_dy+nud_y[j]*d2u_dxy);
    }

    //BC (ny-1)
    const double dzs_dy=inv_dy*(zs[0]-zs[ny-1]);  
    const double d2qx_dxy=inv_dy*(dqx_dx[0]-dqx_dx[ny-1]);                                   
    const double du_dy=inv_dy*(x2c_up*(u[ny]-u[2*ny-1])+x2c_d*(u[0]-u[ny-1]));        
    const double d2u_dxy=inv_dx*inv_dy*((u[ny]-u[2*ny-1])-(u[0]-u[ny-1]));        
    _f[ny-1]=isWater[ny-1]*(forcing_y[ny-1]+qy[ny-1]*inv_dt+g*d_y[ny-1]*(dt*d2qx_dxy -dzs_dy)+dnudy_dx[ny-1]*du_dy+nud_y[ny-1]*d2u_dxy);
    
}

void Hydro::setSolverQyy_d(int i){
    
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;    
    const double inv_dt=1/dt;        
    const double _2inv_dy2=2/(this->grid->dy*this->grid->dy);    
    double  __restrict_arr * d_y=this->Di+(i+1)*ny;    
    double  __restrict_arr * inv_d_y=this->inv_Di+(i+1)*ny;
    double  __restrict_arr * dragTerm=this->DragTerm +i*ny;                    
    double  __restrict_arr * nud_y=this->nuD_y+i*ny;                
    bool  __restrict_arr * isWater=this->IsWater;
    double __restrict_arr* _d=this->d+i*ny;  
    for (int j=0;j<ny;j++)        
    {
        bool isl=!isWater[j];
        _d[j]=(inv_dt+dragTerm[j]+g*d_y[j]*dt*_2inv_dy2+2*nud_y[j]*inv_d_y[j]*_2inv_dy2)*isWater[j]+isl;            
    }


    
    
}

void Hydro::setSolverQyy_up(int i){    
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const double inv_dy=1/this->grid->dy;
    const double inv_dy2=inv_dy*inv_dy;            
    double  __restrict_arr * d_y=this->Di+(i+1)*ny;
    double  __restrict_arr * inv_d_y=this->inv_Di+(i+1)*ny;    
    double  __restrict_arr * v=this->V+(i+1)*ny;            
    double  __restrict_arr * nud_y=this->nuD_y+i*ny;            
    double  __restrict_arr * dnud_dy=this->dnuD_dy+ i*ny;            
    bool  __restrict_arr * isWater=this->IsWater+i*ny;
    double __restrict_arr* _up=this->up+i*ny;  

        for (int j=0;j<ny-1;j++) 
        {
            double _isWater=isWater[j]&&isWater[j+1];                                    
            _up[j]=(0.5*inv_dy*(v[j+1]-2*dnud_dy[j]*inv_d_y[j+1])-inv_dy2*(g*d_y[j]*dt+2*nud_y[j]*inv_d_y[j+1]))*_isWater;             
        }
        double _isWater=isWater[ny-1]&&isWater[0];                                    
        _up[ny-1]=(0.5*inv_dy*(v[0]-2*dnud_dy[ny-1]*inv_d_y[0])-inv_dy2*(g*d_y[ny-1]*dt+2*nud_y[ny-1]*inv_d_y[0]))*_isWater;         
    
}

void Hydro::setSolverQyy_low(int i){
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const double inv_dy=1/this->grid->dy;
    const double inv_dy2=inv_dy*inv_dy;            
    double  __restrict_arr * d_y=this->Di+(i+1)*ny;
    double  __restrict_arr * inv_d_y=this->inv_Di+(i+1)*ny;    
    double  __restrict_arr * v=this->V+(i+1)*ny;            
    double  __restrict_arr * nud_y=this->nuD_y+i*ny;            
    double  __restrict_arr * dnud_dy=this->dnuD_dy+ i*ny;            
    bool  __restrict_arr * isWater=this->IsWater+i*ny;
    double __restrict_arr* _low=this->low+i*ny;  
    for (int j=1;j<ny;j++) 
    {
        double _isWater=isWater[j]&&isWater[j-1];                                                   
        _low[j]=(0.5*inv_dy*(-v[j-1]+2*dnud_dy[j]*inv_d_y[j-1])-inv_dy2*(g*d_y[j]*dt+2*nud_y[j]*inv_d_y[j-1]))*_isWater;         
    }
    double _isWater=isWater[0]&&isWater[ny-1];                                                   
    _low[0]=(0.5*inv_dy*(-v[ny-1]+2*dnud_dy[0]*inv_d_y[ny-1])-inv_dy2*(g*d_y[0]*dt+2*nud_y[0]*inv_d_y[ny-1]))*_isWater;     
              
        
        

        
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          Qy Proyector Solver Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Hydro::setSolverQyx_f(int column){  
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column;
    
    double  __restrict_arr * dragTerm=this->DragTerm+i*ny;
    double  __restrict_arr * Qy=this->Qy+(i+1)*ny;    
    bool  __restrict_arr * isWater=this->IsWater+ i*ny;
    double  __restrict_arr * _f=this->f+i*ny;    
    for (int j=0;j<ny;j++)              
        _f[j]=isWater[j]*Qy[j]*(1+dt*dragTerm[j]);
        //_f[j]=isWater[j]*Qy[j]*(1/dt+dragTerm[j]);
    
}

void Hydro::setSolverQyx_d(int column){
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column; 
    double const x2c_d=this->grid->x2c_d[i];
    double const x2c_up=this->grid->x2c_up[i];
    double  __restrict_arr *  dragTerm=this->DragTerm +i*ny;    
    double  __restrict_arr *  nud=this->nuD_y+i*ny;    
    double  __restrict_arr *  dnud_dx=this->dnuDy_dx+i*ny;    
    double  __restrict_arr *  inv_d=this->inv_Di+(i+1)*ny;    
    bool  __restrict_arr *  isWater=this->IsWater+i*ny;    
    double  __restrict_arr *  u =this->U+i*ny;
    double  __restrict_arr *  _d=this->d+i*ny;
    
    const double dx_d=this->grid->difx_c_d[i];
    const double d2x_d=this->grid->dif2x_c_d[i];
    const double inv_dx=1/this->grid->dx_x[i];

        
        for (int j=0;j<ny-1;j++)      
        {   
            const bool isl=!isWater[j];
            const double du_dx=0.5*inv_dx*((u[j+ny+1]+u[j+ny])-(u[j+1]+u[j]));            
            const double u_y=0.5*(x2c_up*(u[ny+j+1]+u[ny+j])+x2c_d*(u[j+1]+u[j]));                    
            _d[j]=isWater[j]*((1+dt*dragTerm[j])+dt*(u_y*dx_d+du_dx)-dt*(inv_d[j]*nud[j]*d2x_d+inv_d[j]*dnud_dx[j]*dx_d))+isl;
            //_d[j]=isWater[j]*((1/dt+dragTerm[j])+(u_y*dx_d+du_dx)-(inv_d[j]*nud[j]*d2x_d+inv_d[j]*dnud_dx[j]*dx_d))+isl;
        }
        const bool isl=!isWater[ny-1];
        const double du_dx=0.5*inv_dx*((u[ny]+u[2*ny-1])-(u[0]+u[ny-1]));            
        const double u_y=0.5*(x2c_up*(u[ny]+u[2*ny-1])+x2c_d*(u[0]+u[ny-1]));                            
        _d[ny-1]=isWater[ny-1]*((1+dt*dragTerm[ny-1])+dt*(u_y*dx_d+du_dx)-dt*(inv_d[ny-1]*nud[ny-1]*d2x_d+inv_d[ny-1]*dnud_dx[ny-1]*dx_d))+isl;
        //_d[ny-1]=isWater[ny-1]*((1/dt+dragTerm[ny-1])+(u_y*dx_d+du_dx)-(inv_d[ny-1]*nud[ny-1]*d2x_d+inv_d[ny-1]*dnud_dx[ny-1]*dx_d))+isl;
        
        
    
}

void Hydro::setSolverQyx_up(int column){
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column; 
    
    double  __restrict_arr * nud=this->nuD_y+i*ny;    
    double  __restrict_arr * dnud_dx=this->dnuDy_dx+i*ny;    
    double  __restrict_arr * inv_d_y_up=this->inv_Di+(2+i)*ny;        
    bool  __restrict_arr * isWater=this->IsWater+i*ny;    
    double  __restrict_arr * u=this->Ui_j+i*ny;
    double  __restrict_arr * _up=this->up+i*ny;    

    const double dx_up=this->grid->difx_c_up[i];
    const double d2x_up=this->grid->dif2x_c_up[i];
        
    for (int j=0;j<ny;j++)    
    {
        double _isWater=isWater[j]&&isWater[j+ny];                                                                     
        _up[j]=_isWater*dt*(u[j]*dx_up-nud[j]*d2x_up*inv_d_y_up[j]-dnud_dx[j]*dx_up*inv_d_y_up[j]);    
        //_up[j]=isWater[j]*(u[j]*dx_up-nud[j]*d2x_up*inv_d_y_up[j]-dnud_dx[j]*dx_up*inv_d_y_up[j]);    
    }
               
        
        
}

void Hydro::setSolverQyx_low(int column){
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=column; 


    double  __restrict_arr * nud=this->nuD_y+i*ny;    
    double  __restrict_arr * dnud_dx=this->dnuDy_dx+i*ny;    
    double  __restrict_arr * inv_d_y_low=this->inv_Di+i*ny;        
    bool    __restrict_arr * isWater=this->IsWater+i*ny;    
    double  __restrict_arr * u=this->Ui_j+i*ny;
    double  __restrict_arr * _low=this->low+i*ny;    
    
    const double dx_low=this->grid->difx_c_low[i];
    const double d2x_low=this->grid->dif2x_c_low[i];
        
    for (int j=0;j<ny;j++)  
    {                   
         double _isWater=isWater[j]&&isWater[j-ny];      
        _low[j]=_isWater*dt*(u[j]*dx_low-nud[j]*d2x_low*inv_d_y_low[j]-dnud_dx[j]*dx_low*inv_d_y_low[j]);        
        //low[j]=isWater[j]*(u[j]*dx_low-nud[j]*d2x_low*inv_d_y_low[j]-dnud_dx[j]*dx_low*inv_d_y_low[j]);        
    }
        
    
    
}

 void Hydro::setSolverQyx0_low(){
    const int ny=this->grid->ny;
    const double dt=Parameters::Hydro::dt;
    const int i=0; 


    double  __restrict_arr * nud=this->nuD_y+i*ny;    
    double  __restrict_arr * dnud_dx=this->dnuDy_dx+i*ny;    
    double  __restrict_arr * inv_d_y_low=this->inv_Di+i*ny;        
    bool    __restrict_arr * isWater=this->IsWater+i*ny;    
    double  __restrict_arr * u=this->Ui_j+i*ny;
    double  __restrict_arr * _low=this->low+i*ny;    
    
    const double dx_low=this->grid->difx_c_low[i];
    const double d2x_low=this->grid->dif2x_c_low[i];
        
    for (int j=0;j<ny;j++)                      
        _low[j]=isWater[j]*dt*(u[j]*dx_low-nud[j]*d2x_low*inv_d_y_low[j]-dnud_dx[j]*dx_low*inv_d_y_low[j]);        
        
        
    
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          Pre-Solver and Post-Solver Methods 
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Hydro::getForcing_y()
{
    int nx=this->grid->nx;    
    int ny=this->grid->ny;
    double inv_rho=1/rho;
    double inv_dy=1/this->grid->dy;
    inv_dy*=inv_rho;

    double * forcing_y=this->Forcing_y;
    double* Syy=this->wave->Syy;    
    double* Sxy_low=this->wave->Sxy;
    double* Sxy=this->wave->Sxy+ny;
    double* Sxy_up=this->wave->Sxy+2*ny;
    double* difx_c_d=this->grid->difx_c_d;
    double* difx_c_up=this->grid->difx_c_up;
    double* difx_c_low=this->grid->difx_c_low;
    

    //Primera Fila
    double d_dx_d=difx_c_d[0]*inv_rho;    
    double d_dx_up=difx_c_up[0]*inv_rho;    
    double d_dx_low=difx_c_low[0]*inv_rho;    
    for (int j=0;j<ny-1;j++)  
    {          
        forcing_y[j]=-(inv_dy*(Syy[j+1]-Syy[j])+0.5*((Sxy_up[j]+Sxy_up[j+1])*d_dx_up+(Sxy[j]+Sxy[j+1])*d_dx_d+(Sxy_low[j]+Sxy_low[j+1])*d_dx_low));           
    }
    forcing_y[ny-1]=-(inv_dy*(Syy[0]-Syy[ny-1])+0.5*((Sxy_up[0]+Sxy_up[ny-1])*d_dx_up+(Sxy[0]+Sxy[ny-1])*d_dx_d+(Sxy_low[0]+Sxy_low[ny-1])*d_dx_low));

    for (int i=1;i<nx-1;i++){
        d_dx_d=difx_c_d[i]*inv_rho;    
        d_dx_up=difx_c_up[i]*inv_rho;    
        d_dx_low=difx_c_low[i]*inv_rho;          

        Syy+=ny;
        forcing_y+=ny;        
        for (int j=0;j<ny-1;j++)    
            forcing_y[j]=-(inv_dy*(Syy[j+1]-Syy[j])+0.5*((Sxy_up[j]+Sxy_up[j+1])*d_dx_up+(Sxy[j]+Sxy[j+1])*d_dx_d+(Sxy_low[j]+Sxy_low[j+1])*d_dx_low));           
           
        forcing_y[ny-1]=-(inv_dy*(Syy[0]-Syy[ny-1])+0.5*((Sxy_up[0]+Sxy_up[ny-1])*d_dx_up+(Sxy[0]+Sxy[ny-1])*d_dx_d+(Sxy_low[0]+Sxy_low[ny-1])*d_dx_low));
        
        Sxy_low=Sxy;
        Sxy=Sxy_up;
        Sxy_up+=ny;
    }

    //Ultima Fila
    Syy+=ny;    
    forcing_y+=ny;        
    d_dx_d=difx_c_d[nx-1]*inv_rho;        
    d_dx_low=difx_c_low[nx-1]*inv_rho;
    
    for (int j=0;j<ny-1;j++)    
            forcing_y[j]=-(inv_dy*(Syy[j+1]-Syy[j])+0.5*((Sxy[j]+Sxy[j+1])*d_dx_d+(Sxy_low[j]+Sxy_low[j+1])*d_dx_low));
           
    forcing_y[ny-1]=-(inv_dy*(Syy[0]-Syy[ny-1])+0.5*((Sxy[0]+Sxy[ny-1])*d_dx_d+(Sxy_low[0]+Sxy_low[ny-1])*d_dx_low));        

    
}

void Hydro::getD_y0()
{
    const int ny=this->grid->ny;
    double const __restrict_arr * d=this->guess_D;      
    double  __restrict_arr * d_y= this->Di;
    double  __restrict_arr * inv_d_y= this->inv_Di;
    //Primera linea
     for (int j=0;j<ny-1;j++)   
     { 
        d_y[j]=1*(d[j+1]+d[j])-0.5*(d[ny+j+1]+d[ny+j]); 
        inv_d_y[j]=1/d_y[j];
     }
     d_y[ny-1]=1*(d[0]+d[ny-1])-0.5*(d[0]+d[2*ny-1]); 
     inv_d_y[ny-1]=1/d_y[ny-1];
}

void Hydro::getD_y(int i)
{
     const int ny=this->grid->ny;
     double  __restrict_arr * d=this->guess_D+i*ny;      
     double  __restrict_arr * d_y= this->Di+(i+1)*ny;
     double  __restrict_arr * inv_d_y= this->inv_Di+(i+1)*ny;

      for (int j=0;j<ny-1;j++)   
      { 
         d_y[j]=0.5*(d[j+1]+d[j]); 
         inv_d_y[j]=1/d_y[j];
      }
        d_y[ny-1]=0.5*(d[0]+d[ny-1]); 
         inv_d_y[ny-1]=1/d_y[ny-1];
}

void Hydro::getV0()
{
    const int ny=this->grid->ny;
    double   __restrict_arr * d=this->bathy->D;    
    double   __restrict_arr * qy=this->Qy;        
    double   __restrict_arr * v=this->V;        
    
    for (int j=0;j<ny;j++)    {
        const double d_y0=2*d[j]-d[j+ny];    
        v[j]=qy[j]/d_y0;        
    }    
}

 void Hydro::getV(int column)
{

    const int ny=this->grid->ny;
    const int i=column;
    double   __restrict_arr * d=this->bathy->D+i*ny;    
    double   __restrict_arr * qy=this->Qy+(i+1)*ny;        
    double   __restrict_arr * v=this->V+(i+1)*ny;            
    
    for (int j=0;j<ny-1;j++)    {
        const double d_y=0.5*(d[j]+d[j+1]);    
        v[j]=qy[j]/d_y;          
    }    
    const double d_y=0.5*(d[ny-1]+d[0]);    
    v[ny-1]=qy[ny-1]/d_y;      

}