/*******************************************************************************************
  Wave Dir, Computes the wave ray refraccion calculating rot(K)=0 in a ray form
        -solveDir: Core method
        - solvedir_paralle: call solveDir in order to compute pieces of the row for each processor.     

          Methods and parameters description in "WaveRow.h" 

*******************************************************************************************/
#include "WaveRow.h"


void WaveRow::solveDir(const double* _tanTheta_m1,const double* _lK_m1,const double* _lK,int _indexIni,int _chunk){
    int chunk=_chunk;
    int indexIni=_indexIni;
    double dx=this->dx;
    double inv_2dy = this->inv_2dy; //  1/(2*dy)
    
    double* tanTheta_m1=new double[chunk+2];
    double* lK_m1=new double[chunk+2];
    double* lK=new double[chunk+2];
    std::copy(_tanTheta_m1+indexIni, _tanTheta_m1+indexIni+chunk+2, tanTheta_m1);
    std::copy(_lK_m1+indexIni, _lK_m1+indexIni+chunk+2, lK_m1);
    std::copy(_lK+indexIni, _lK+indexIni+chunk+2, lK);
    double *Cx=new double[chunk];
    double *Cgx=new double[chunk];
    double *Cy=new double[chunk];
    double *Cgy=new double[chunk];
    std::copy(this->C+indexIni, this->C+indexIni+chunk, Cx);
    std::copy(this->Cg+indexIni, this->Cg+indexIni+chunk, Cgx);
    

    double * Y_i=new double[chunk];
    double * inv_S_i=new double[chunk];
    double * Theta=new double[chunk];
    double * TanTheta=new double[chunk];
    double * CosTheta=new double[chunk];
    double * SinTheta=new double[chunk];
    double * sxx=new double[chunk];
    double * syy=new double[chunk];
    double * sxy=new double[chunk];
    double * sxx_r=new double[chunk];
    double * syy_r=new double[chunk];
    double * sxy_r=new double[chunk];
    double dx2 = dx * dx;
    double inv_dx = 1 / dx;

    const double thetaMax=Parameters::Wave::ThetaMax;
    


    
    for (int i = 1; i < chunk+1; i++)
    {
        //Get incidente Ray 
        double tt_p1 = tanTheta_m1[i+1];
        double tt_m1 = tanTheta_m1[i-1];
        double dttheta_dy = (tt_p1 - tt_m1) * inv_2dy;
        double ttheta_med = (tt_p1 + tt_m1) * 0.5;
        double y_i = -ttheta_med / (dttheta_dy + inv_dx); //y-distance between center node and incidente position
        double s_i = sqrt(dx2 + y_i * y_i); //lenght of the ray
        double inv_s_i = 1 / s_i;
        double sinTheta_i = -y_i * inv_s_i;
        double cosTheta_i = dx * inv_s_i;
        double theta_i = asin(sinTheta_i);
        Y_i[i-1] = y_i;
        inv_S_i[i-1] = inv_s_i;

        double dlk_dx=(lK[i]-lK_m1[i])*inv_dx;

        double dlk_dy_c=0.5*inv_2dy*(lK[i+1]+lK_m1[i+1]-lK[i-1]-lK_m1[i-1]);        
        double dlk_ds = dlk_dx * cosTheta_i + dlk_dy_c * sinTheta_i;
        double dlk_dst = -dlk_dx * sinTheta_i + dlk_dy_c * cosTheta_i;
        double theta = theta_i + dlk_dst / (inv_s_i + dlk_ds * 0.5);

        theta=min(theta,thetaMax);
        theta=max(theta,-thetaMax);
        
         Theta[i-1]=theta;
         double sinTheta=sin(theta);
         double cosTheta=cos(theta);
         SinTheta[i-1]=sinTheta;
         CosTheta[i-1]=cosTheta;
         TanTheta[i-1]=sinTheta/cosTheta;

        double c=Cx[i-1];
        double cg=Cgx[i-1];
        Cx[i-1]*=cosTheta;
        Cgx[i-1]*=cosTheta;
        Cy[i-1]=c*sinTheta;
        Cgy[i-1]=cg*sinTheta;
        double cg_c=cg/c;
        double cosTheta2=cosTheta*cosTheta;
        double cossinTheta=cosTheta*sinTheta;
        sxx[i-1]=((1+cosTheta2)*cg_c-0.5);
        syy[i-1]=((2-cosTheta2)*cg_c-0.5);
        sxy[i-1]=(cossinTheta*cg_c);
        sxx_r[i-1]=2*cosTheta2;
        syy_r[i-1]=2*(1-cosTheta2);
        sxy_r[i-1]=2*cossinTheta;  
       
    }

     std::copy(SinTheta, SinTheta + chunk, this->SinTheta+indexIni);
     std::copy(CosTheta, CosTheta + chunk, this->CosTheta+indexIni);
     std::copy(TanTheta, TanTheta + chunk, this->TanTheta+indexIni);
     std::copy(Theta, Theta + chunk, this->Theta+indexIni);
     std::copy(Y_i, Y_i + chunk, this->Y_i+indexIni);
     std::copy(inv_S_i, inv_S_i + chunk, this->inv_S_i+indexIni);
    std::copy(Cx, Cx + chunk, this->Cx+indexIni);
    std::copy(Cy, Cy + chunk, this->Cy+indexIni);
    std::copy(Cgx, Cgx + chunk, this->Cgx+indexIni);
    std::copy(Cgy, Cgy + chunk, this->Cgy+indexIni);
    std::copy(sxx, sxx + chunk, this->sxx+indexIni);
    std::copy(sxy, sxy + chunk, this->sxy+indexIni);
    std::copy(syy, syy + chunk, this->syy+indexIni);
    std::copy(sxx_r, sxx_r + chunk, this->sxx_r+indexIni);
    std::copy(sxy_r, sxy_r + chunk, this->sxy_r+indexIni);
    std::copy(syy_r, syy_r + chunk, this->syy_r+indexIni);
     delete[] lK;delete[] lK_m1;delete[] tanTheta_m1;
     delete[] Theta;delete[] TanTheta;
     delete[] SinTheta;delete[] CosTheta;
     delete[] Y_i;delete[] inv_S_i;
     delete[] Cx;delete[] Cy;
     delete[] Cgx;delete[] Cgy;
    delete[] sxx;delete[] sxy;delete[] syy;
    delete[] sxx_r;delete[] sxy_r;delete[] syy_r;

  
}



void WaveRow::solveDir_parallel(const WaveRow* rm1,const int* _indexIni, const int* _chunk,int np){
     
     int* chunk=new int[np];
     int* indexIni=new int[np]();
     std::copy(_chunk,_chunk+np, chunk);
     std::copy(_indexIni,_indexIni+np, indexIni);

    double* tanTheta_per=new double[ny+2]();
    std::copy(rm1->TanTheta,rm1->TanTheta+ny,tanTheta_per+1);
    tanTheta_per[0]=tanTheta_per[ny];
    tanTheta_per[ny+1]=tanTheta_per[1];
    double* lKm1_per=new double[ny+2]();
    std::copy(rm1->lK,rm1->lK+ny,lKm1_per+1);
    lKm1_per[0]=lKm1_per[ny];
    lKm1_per[ny+1]=lKm1_per[1];
    double* lK_per=new double[ny+2]();
    std::copy(this->lK,this->lK+ny,lK_per+1);
    lK_per[0]=lK_per[ny];
    lK_per[ny+1]=lK_per[1];

     #pragma omp parallel shared(chunk,indexIni,tanTheta_per,lKm1_per,lK_per)  
    {
       int pid=Omp::get_pid();
       solveDir(tanTheta_per,lKm1_per,lK_per,indexIni[pid],chunk[pid]);
    }

        delete[] chunk;delete[] indexIni;
        delete[]tanTheta_per;delete[] lK_per;delete[] lKm1_per; 

}




