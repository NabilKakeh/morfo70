/*********************************************************************************
 * 
 *                  Math.cpp (Morfo70)
 * 
 *  Implementation of linear solvers A*y=f
 *  
 *  Descripcion in header Math.h 
 * 
 * *******************************************************************************/
#include "Math.h"

namespace Math
{
        
        void TridiagonalSolver(double * Y,double* V,double* A, double* B, double* C, double * F, int n,int m){
            const double * a=A;
            double* y=Y;
            const double * f=F;
            double * inv_w=new double[m];
            for (int j=0;j<m;j++) 
            {       
                inv_w[j]=1/a[j];            
                y[j]=f[j]*inv_w[j];           
            
            }
            double * v=V;
            double * ylow;
            const double * c=C;
            const double * b=B;        
            for (int i=1;i<n;i++)
            {
                ylow=y;
                a+=m;
                b+=m;
                y+=m;
                f+=m;
                
                for (int j=0;j<m;j++)
                {
                    v[j]=c[j]*inv_w[j];
                    inv_w[j]=1/(a[j]-b[j]*v[j]);
                    y[j]=(f[j]-b[j]*ylow[j])*inv_w[j];
                }
                c+=m;
                v+=m;
                
                
            }
            ylow=Y+(n-2)*m;
            y=ylow+m;
            v=V+(n-2)*m;
            
            for (int i=n-2;i>-1;i--)
            {            
                for (int j=0;j<m;j++)
                {
                    ylow[j]-=v[j]*y[j];
                }
                v-=m;
                y=ylow;
                ylow-=m;
                
            }
            delete[] inv_w;
        }
        
        
        void TridiagonalSolverPeriodic(double* a, double* b, double* c, double * f,  int n,double * _y){
            double t=b[0];
            double s=c[n-1];
            double* k=new double[n];
            double* inv_k=new double[n];
            double* v=new double[n];
            double* h=new double[n];       

            k[0]=a[0];
            v[0]=t;
            inv_k[0]=1/k[0];
            h[0]=s*inv_k[0];
            for (int i=1;i<n-2;i++)
            {                
                k[i]=a[i]-b[i]*c[i-1]*inv_k[i-1];
                inv_k[i]=1/k[i];
                v[i]=-b[i]*v[i-1]*inv_k[i-1];                
                h[i]=-c[i-1]*h[i-1]*inv_k[i];
            }
            
            k[n-2]=a[n-2]-b[n-2]*c[n-3]*inv_k[n-3];        
            inv_k[n-2]=1/k[n-2];
            v[n-2]=c[n-2]-b[n-2]*v[n-3]*inv_k[n-3];            
            h[n-2]=(b[n-1]-h[n-3]*c[n-3])*inv_k[n-2];
             


            double* r=new double[n];
            double* y=new double[n];
            r[0]=f[0];
            double rh=r[0]*h[0];
            double vh=v[0]*h[0];
            for (int i=1;i<n-1;i++)            {
                r[i]=f[i]-b[i]*r[i-1]*inv_k[i-1];    
                rh+=h[i]*r[i];
                vh+=h[i]*v[i];
            }
            k[n-1]=a[n-1]-vh;
            r[n-1]=f[n-1]-rh;

            y[n-1]=r[n-1]/k[n-1];
            y[n-2]=(r[n-2]-v[n-2]*y[n-1])*inv_k[n-2];

            for (int i=n-3;i>-1;i--)
                y[i]=(r[i]-c[i]*y[i+1]-v[i]*y[n-1])*inv_k[i];
            
            std::copy(y,y+n,_y);
            delete[] y;delete[] k;delete[] inv_k;delete[] r; delete[]v;delete[]h;


        }
        

        void TridiagonalSolverPeriodic(double * Y, double* A, double* B,  double* C, double * F, int n,int m){
            double* k=new double[n];
            double* inv_k=new double[n];
            double* v=new double[n];
            double* h=new double[n];
            double* r=new double[n];

            

            double* a=A;
            double* b=B;
            double* c=C;
            double* f=F;
            double* y=Y;
            for (int i=0;i<m;i++)
            {

                double t=b[0];
                double s=c[n-1];
                    

                k[0]=a[0];
                v[0]=t;
                inv_k[0]=1/k[0];
                h[0]=s*inv_k[0];
                for (int j=1;j<n-2;j++)
                {                    
                    k[j]=a[j]-b[j]*c[j-1]*inv_k[j-1];
                    inv_k[j]=1/k[j];
                    v[j]=-b[j]*v[j-1]*inv_k[j-1];                
                    h[j]=-c[j-1]*h[j-1]*inv_k[j];
                }
                
                k[n-2]=a[n-2]-b[n-2]*c[n-3]*inv_k[n-3];        
                inv_k[n-2]=1/k[n-2];
                v[n-2]=c[n-2]-b[n-2]*v[n-3]*inv_k[n-3];            
                h[n-2]=(b[n-1]-h[n-3]*c[n-3])*inv_k[n-2];


                
                
                r[0]=f[0];
                double rh=r[0]*h[0];
                double vh=v[0]*h[0];
                for (int j=1;j<n-1;j++)            {
                    r[j]=f[j]-b[j]*r[j-1]*inv_k[j-1];    
                    rh+=h[j]*r[j];
                    vh+=h[j]*v[j];
                }
                k[n-1]=a[n-1]-vh;
                r[n-1]=f[n-1]-rh;

                y[n-1]=r[n-1]/k[n-1];
                y[n-2]=(r[n-2]-v[n-2]*y[n-1])*inv_k[n-2];

                for (int j=n-3;j>-1;j--)
                    y[j]=(r[j]-c[j]*y[j+1]-v[j]*y[n-1])*inv_k[j];            

                a+=n;
                b+=n;
                c+=n;
                f+=n;
                y+=n;
            }
            delete[] k;delete[] inv_k;delete[] r; delete[]v;delete[]h;
        }








}
