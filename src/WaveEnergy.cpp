/*******************************************************************************************
*  Wave Energy,   solves the Wave Energy balance Equation (WEBEQ) with wave using ray-trace method:
*        Cg*grad(lE) + div(Cg)+Dw/E+Df/E=0
*
*  Methods and parameters description in "WaveRow.h" 
*
*******************************************************************************************/

#include "WaveRow.h"



void WaveRow::solveEnergy_parallel(const int* _indexIni, const int* _chunk,int np)
{
    
     int* chunk=new int[np];
     int* indexIni=new int[np]();
     std::copy(_chunk,_chunk+np, chunk);
     std::copy(_indexIni,_indexIni+np, indexIni);
    
    

     #pragma omp parallel shared(indexIni,chunk)  
    {
        int pid=Omp::get_pid();
    //  int pid=0;
        
       solveEnergy(indexIni[pid],chunk[pid]);
    }
    delete[] chunk;delete[] indexIni;
}


void WaveRow::solveEnergy(int indexIni, int chunk)
{
     
    double * G =new double[chunk];
    std::copy(this->G,this->G+chunk,G);
    double * Gp =new double[chunk];
    std::copy(this->Gp,this->Gp+chunk,Gp);    
    
    double * H =new double[chunk];
    double * E =new double[chunk]();    
    double * Dw =new double[chunk]();        
    double * lE =new double[chunk];
    std::copy(this->lE+indexIni,this->lE+indexIni+chunk,lE);
    std::copy(this->H+indexIni,this->H+indexIni+chunk,H);
    double * D =new double[chunk];
    std::copy(this->D+indexIni,this->D+indexIni+chunk,D);
    double * Sigma =new double[chunk];    
    std::copy(this->Sigma+indexIni,this->Sigma+indexIni+chunk,Sigma);
    
    double * Urms_H =new double[chunk];
    std::copy(this->Urms_H+indexIni,this->Urms_H+indexIni+chunk,Urms_H);
    int solverMaxIterations=Parameters::Wave::SolverMaxIterations;
    double solverTolerance=Parameters::Wave::SolverTolerance;    
    
    
    
    int * solveNode=new int[chunk];
    int * nodesToSolve=new int[chunk];    

    double  logCteEnergy=log(cteEnergy);    
    

     int j=0;
     for (j=0;j<chunk;j++)
     { 
         solveNode[j]=1;    
     }
      int nNodesToSolve=chunk;
      double * dw =new double[chunk]();    
      double * ddw_dh =new double[chunk]();    
      double * df =new double[chunk]();    
      double * ddf_dh =new double[chunk]();    

    // Main Loop. Each iteration, checks if the node converges  and it is removed from the vector notesToSolve
    for (int k=0;k<solverMaxIterations;k++)
    {
         nNodesToSolve=0;
         for (j=0;j<chunk;j++){
             if (solveNode[j])
             {
                 nodesToSolve[nNodesToSolve]=j;
                 nNodesToSolve+=1;
             }

         }
         if (nNodesToSolve==0) break;
       
        //Wave disipation 
         (this->*waveBreakingDissipation)(chunk,H,Sigma,D,nodesToSolve,nNodesToSolve,dw,ddw_dh);
          if (Options::Wave::BottomFriction)
            (this->*BottomFriction)(chunk,H,Sigma[0],Urms_H,nodesToSolve,nNodesToSolve,df,ddf_dh);        
          
          //Solves iteration with Newton raphson gg-> objective function ; gp-> d(gg)_d(le) derivative
          for (int i=0;i<nNodesToSolve;i++)
          {            
             
             int nodeIndex=nodesToSolve[i];
             
              double gg=G[nodeIndex];
              double gp=Gp[nodeIndex];
              ;
              double le=lE[nodeIndex];
              
              gg+=gp*le+0.5*(dw[nodeIndex]+df[nodeIndex]);                           
              gp+=gp+0.25*(ddf_dh[nodeIndex]+ddw_dh[nodeIndex])*H[nodeIndex];                      

              le-=gg/gp;              
              double h=exp((le-logCteEnergy)*0.5);                                        
              H[nodeIndex]=h;                            
              E[nodeIndex]=cteEnergy*h*h;
              
              //Check convergence
              solveNode[nodeIndex]=abs(le-lE[nodeIndex])>solverTolerance;
              lE[nodeIndex]=le;
              Dw[nodeIndex]=dw[nodeIndex]*E[nodeIndex];
          }

    }

    //Copy and remove dummies variables 
       std::copy(lE,lE+chunk,this->lE+indexIni);
       std::copy(E,E+chunk,this->E+indexIni);
       std::copy(H,H+chunk,this->H+indexIni);
       std::copy(dw,dw+chunk,this->dw+indexIni);
       std::copy(df,df+chunk,this->df+indexIni);
       std::copy(Dw,Dw+chunk,this->Dw+indexIni);
       delete[]lE;delete[] H;delete[] D;delete[] Sigma;delete[] E;
       delete[]G;
       delete[] Gp;
       delete[] Urms_H;
       delete[] dw;       delete[] ddw_dh;
       delete[] ddf_dh;   delete[] df;       
       delete[] Dw;
       delete[] nodesToSolve;delete[] solveNode;
}


void WaveRow::energyPreprocessor(const WaveRow* rm1)
{
    
     double inv_dx=1/this->dx;
     double inv_2dy=this->inv_2dy;
     int ny=this->ny;
     int* indexUpY=new int[ny];
     int* indexLowY=new int[ny];
     std::copy(this->indexUpY,this->indexUpY+ny,indexUpY);
     std::copy(this->indexLowY,this->indexLowY+ny,indexLowY);
     double *Y_i=new double[ny];
     std::copy(this->Y_i,this->Y_i+ny,Y_i);
     double *inv_S_i=new double[ny];
     std::copy(this->inv_S_i,this->inv_S_i+ny,inv_S_i);
    double *Cg=new double[ny];
    std::copy(this->Cg,this->Cg+ny,Cg);
     double *Cg_m1=new double[ny];     
     std::copy(rm1->Cg,rm1->Cg+ny,Cg_m1);
     double *dw_m1=new double[ny];
     std::copy(rm1->dw,rm1->dw+ny,dw_m1);
     double *df_m1=new double[ny];
     std::copy(rm1->df,rm1->df+ny,df_m1);
     double *lE_m1=new double[ny];
     std::copy(rm1->lE,rm1->lE+ny,lE_m1);    
     double *Cgx_m1=new double[ny];
     std::copy(rm1->Cgx,rm1->Cgx+ny,Cgx_m1);
     double *Cgy_m1=new double[ny];
     std::copy(rm1->Cgy,rm1->Cgy+ny,Cgy_m1);
     double *Cgx=new double[ny];
     std::copy(this->Cgx,this->Cgx+ny,Cgx);
     double *Cgy=new double[ny];
     std::copy(this->Cgy,this->Cgy+ny,Cgy);
     double *U=new double[ny];
     std::copy(this->U,this->U+ny,U);
     double *U_m1=new double[ny];
     std::copy(rm1->U,rm1->U+ny,U_m1);
     double *V=new double[ny];
     std::copy(this->V,this->V+ny,V);;
     double *V_m1=new double[ny];
     std::copy(rm1->V,rm1->V+ny,V_m1);
     double *sxx=new double[ny];
     double *sxy=new double[ny];
     double *syy=new double[ny];
     std::copy(this->sxx,this->sxx+ny,sxx);
     std::copy(this->sxy,this->sxy+ny,sxy);
     std::copy(this->syy,this->syy+ny,syy);
    double *sxx_r=new double[ny];
     double *sxy_r=new double[ny];
     double *syy_r=new double[ny];
     std::copy(this->sxx_r,this->sxx_r+ny,sxx_r);
     std::copy(this->sxy_r,this->sxy_r+ny,sxy_r);
     std::copy(this->syy_r,this->syy_r+ny,syy_r);
     double *WC_m1=new double[ny];
     std::copy(rm1->WC,rm1->WC+ny,WC_m1);

     double *G=new double[ny];
     double *Gp=new double[ny];     
     double *WC=new double[ny];
     double *WC_r=new double[ny];

     

    #pragma omp parallel for
     for (int i = 0; i < ny; i++)
     {
        
         int indexUp = indexUpY[i];
         int indexLow = indexLowY[i];
         double y_i_dy = Y_i[i] * inv_2dy;
         double inv_s_i=inv_S_i[i];

        
        double cg_c=0.5*(Cg[i] +(Cg_m1[indexUp] + Cg_m1[indexLow]) * 0.5 + (Cg_m1[indexUp] - Cg_m1[indexLow])*y_i_dy);
        double u_c=0.5*(U[i] +(U_m1[indexUp] + U_m1[indexLow]) * 0.5 + (U_m1[indexUp] - U_m1[indexLow])*y_i_dy);
        double v_c=0.5*(V[i] +(V_m1[indexUp] + V_m1[indexLow]) * 0.5 + (V_m1[indexUp] - V_m1[indexLow])*y_i_dy);
        double dw_i=(dw_m1[indexUp] + dw_m1[indexLow]) * 0.5 + (dw_m1[indexUp] - dw_m1[indexLow])*y_i_dy;
        double df_i=(df_m1[indexUp] + df_m1[indexLow]) * 0.5 + (df_m1[indexUp] - df_m1[indexLow])*y_i_dy;
        double le_i=(lE_m1[indexUp] + lE_m1[indexLow]) * 0.5 + (lE_m1[indexUp] - lE_m1[indexLow])*y_i_dy;
        
        
        

         //WC interaction  
         double du_dx=(U[i]-U_m1[i])*inv_dx;
         double du_dy=(U[indexUp]-U[indexLow])*inv_2dy;
         double dv_dx=(V[i]-V_m1[i])*inv_dx;
         double dv_dy=(V[indexUp]-V[indexLow])*inv_2dy;
         double wc=2*sxx[i]*du_dx+sxy[i]*(du_dy+dv_dx)+2*syy[i]*dv_dy; 
         double wc_c=0.5*(wc +(WC_m1[indexUp] + WC_m1[indexLow]) * 0.5 + (WC_m1[indexUp] - WC_m1[indexLow])*y_i_dy); 

         ////////////////////////////////////////////////
         //////////// Quitar WC
        //  du_dx=0; //NABIL Quitar WC
        //  du_dy=0;  //NABIL Quitar WC
        //  dv_dx=0;  //NABIL Quitar WC
        //  dv_dy=0;  //NABIL Quitar WC
        //  wc=0;  //NABIL Quitar WC
        //  u_c=0; //NABIL Quitar WC
        //  v_c=0; //NABIL Quitar WC         
        //  wc_c=0; //NABIL Quitar WC
         //////////// END Quitar WC
         ////////////////////////////////////////////////
         
         WC[i]=wc;

         

         //WC interaction (ROLLERS)         
         WC_r[i]=sxx_r[i]*du_dx+sxy_r[i]*(du_dy+dv_dx)+syy_r[i]*dv_dy; 
         

         double dcgxu_dx=(Cgx[i]-Cgx_m1[i])*inv_dx+du_dx;         

         double dcgyv_dy_c=0.5*(Cgy[indexUp]+Cgy_m1[indexUp]-Cgy[indexLow]-Cgy_m1[indexLow])*inv_2dy+dv_dy;
         
         double div_cgv_c=dcgxu_dx+dcgyv_dy_c;     
         //double dle_dy_c=0.5*(lE[indexUp]+lE_m1[indexUp]-lE[indexLow]-lE_m1[indexLow])*inv_2dy;
         double dle_dy_c=(lE_m1[indexUp]-lE_m1[indexLow])*inv_2dy;
         
         G[i]=v_c*dle_dy_c+wc_c+div_cgv_c+0.5*(dw_i+df_i)-cg_c*inv_s_i*le_i-u_c*lE_m1[i]*inv_dx;    
         Gp[i]=cg_c*inv_s_i+u_c*inv_dx;
         
     }
   

    std::copy(G,G+ny,this->G);
    std::copy(Gp,Gp+ny,this->Gp);
    // std::copy(Gp0,Gp0+ny,this->Gp0);
    // std::copy(Gp1,Gp1+ny,this->Gp1);
    std::copy(WC,WC+ny,this->WC);
    std::copy(WC_r,WC_r+ny,this->WC_r);
    std::copy(lE_m1,lE_m1+ny,this->lE);
    delete[] indexUpY;delete[] indexLowY;
    delete[] Y_i;delete[] inv_S_i;
    delete[] Cg;delete[] Cgx;delete[] Cgy;
    delete[] Cg_m1; delete[] Cgx_m1;delete[] Cgy_m1;
    delete[] dw_m1;    delete[] df_m1;
    delete[] G;    delete[] Gp;
    //delete[] Gp0;delete[] Gp1;
    delete[] lE_m1;
    //delete[] lE;
    delete[] U_m1; delete[] V_m1;
    delete[] U; delete[] V;
    delete[] sxx;delete[] sxy;delete[] syy;
    delete[] sxx_r;delete[] sxy_r;delete[] syy_r;
    delete[] WC;delete[] WC_m1;    delete[] WC_r;

}

