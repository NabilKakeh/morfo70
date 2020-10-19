/*******************************************************************************************
*  Wave Rollers, solves the Wave Rollers Energy balance Equation (WREBEQ) with wave using ray-trace method:
*        2C*grad(lEr) + div(C)+Dr/E=-alpha_R*Dw
*
*  Methods and parameters description in "WaveRow.h" 
*******************************************************************************************/

#include "WaveRow.h"

void WaveRow::solveRollers(const WaveRow* rm1)
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
    double *C=new double[ny];
    std::copy(this->C,this->C+ny,C);
     double *C_m1=new double[ny];     
     std::copy(rm1->C,rm1->C+ny,C_m1);
      double *U=new double[ny];
     std::copy(this->U,this->U+ny,U);
     double *U_m1=new double[ny];
     std::copy(rm1->U,rm1->U+ny,U_m1);
     double *V=new double[ny];
     std::copy(this->V,this->V+ny,V);;
     double *V_m1=new double[ny];
     std::copy(rm1->V,rm1->V+ny,V_m1);
     double *Cx_m1=new double[ny];
     std::copy(rm1->Cx,rm1->Cx+ny,Cx_m1);
     double *Cy_m1=new double[ny];
     std::copy(rm1->Cy,rm1->Cy+ny,Cy_m1);
     double *Cx=new double[ny];
     std::copy(this->Cx,this->Cx+ny,Cx);
     double *Cy=new double[ny];
     std::copy(this->Cy,this->Cy+ny,Cy);
     double *E_r_m1=new double[ny];
     std::copy(rm1->E_r,rm1->E_r+ny,E_r_m1);    
    double *WC_r=new double[ny];
     std::copy(this->WC_r,this->WC_r+ny,WC_r);
      double *WC_r_m1=new double[ny];
     std::copy(rm1->WC_r,rm1->WC_r+ny,WC_r_m1);
     double *Dw_m1=new double[ny];
     std::copy(rm1->Dw,rm1->Dw+ny,Dw_m1); 
     double *Dw=new double[ny];
     std::copy(this->Dw,this->Dw+ny,Dw); 
     double *d_r_m1=new double[ny];
     std::copy(rm1->d_r,rm1->d_r+ny,d_r_m1);      
     
     
     double* d_r=new double[ny]; //D_r/E
     double* E_r=new double[ny];
     double d_r_div_c=2*g*Parameters::Rollers::sinBeta_r;    //D_r=2*g*sin(Beta_r)*c
     double alpha_r=Parameters::Rollers::alpha_r; //Breaking energy transfer to rollers
    
     #pragma omp parallel for
     for (int i = 0; i < ny; i++)
     {
         int indexUp = indexUpY[i];
         int indexLow = indexLowY[i];
         double y_i_dy = Y_i[i] * inv_2dy;
         double inv_s_i=inv_S_i[i];

         double c_c=0.5*(C[i] +(C_m1[indexUp] + C_m1[indexLow]) * 0.5 + (C_m1[indexUp] - C_m1[indexLow])*y_i_dy);
         double u_c=0.5*(U[i] +(U_m1[indexUp] + U_m1[indexLow]) * 0.5 + (U_m1[indexUp] - U_m1[indexLow])*y_i_dy);
         double v_c=0.5*(V[i] +(V_m1[indexUp] + V_m1[indexLow]) * 0.5 + (V_m1[indexUp] - V_m1[indexLow])*y_i_dy);

         //WC interaction  
         double du_dx=(U[i]-U_m1[i])*inv_dx;                  
         double dv_dy=(V[indexUp]-V[indexLow])*inv_2dy;         
         double wc_r_c=0.5*(WC_r[i] +(WC_r_m1[indexUp] + WC_r_m1[indexLow]) * 0.5 + (WC_r_m1[indexUp] - WC_r_m1[indexLow])*y_i_dy);

         double dcxu_dx=(Cx[i]-Cx_m1[i])*inv_dx+du_dx;
         double dcyv_dy_c=0.5*(Cy[indexUp]+Cy_m1[indexUp]-Cy[indexLow]-Cy_m1[indexLow])*inv_2dy+dv_dy;
         double div_cv_c=dcxu_dx+dcyv_dy_c;     
         double der_dy_c=(E_r_m1[indexUp]-E_r_m1[indexLow])*inv_2dy;

         d_r[i]=d_r_div_c/C[i];  
         double d_r_c=0.5*(d_r[i] +(d_r_m1[indexUp] + d_r_m1[indexLow]) * 0.5 + (d_r_m1[indexUp] - d_r_m1[indexLow])*y_i_dy);
          
         double Dw_c=0.5*(Dw[i] +(Dw_m1[indexUp] + Dw_m1[indexLow]) * 0.5 + (Dw_m1[indexUp] - Dw_m1[indexLow])*y_i_dy);
         double er_i=(E_r_m1[indexUp] + E_r_m1[indexLow]) * 0.5 + (E_r_m1[indexUp] - E_r_m1[indexLow])*y_i_dy;

         double up=c_c*inv_s_i*er_i+u_c*E_r_m1[i]*inv_dx-v_c*der_dy_c+0.5*alpha_r*Dw_c;
         double low=c_c*inv_s_i+u_c*inv_dx+div_cv_c+0.5*wc_r_c+0.5*d_r_c;

         E_r[i]=up/low;

     }
    std::copy(E_r,E_r+ny,this->E_r);
    std::copy(d_r,d_r+ny,this->d_r);    
     delete[] indexUpY;delete[] indexLowY;
    delete[] Y_i;delete[] inv_S_i;
    delete[] C;delete[] Cx;delete[] Cy;
    delete[] C_m1; delete[] Cx_m1;delete[] Cy_m1;
    delete[] Dw_m1;delete[] Dw;        
    delete[] d_r;delete[] d_r_m1;
    delete[] E_r_m1;delete[] E_r;    
    delete[] U_m1; delete[] V_m1;
    delete[] U; delete[] V;    
    delete[] WC_r;delete[] WC_r_m1;
     
}