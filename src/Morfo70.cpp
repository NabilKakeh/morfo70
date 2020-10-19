/*********************************************************************************
 * 
 *                  Morfo70.cpp (Morfo70)
 * 
 *  Global Morfo70 funcions in order to perform numerical simulation.
 *  Main Method: Run
 *
 *  Methods and parameters description in "Morfo70.h"  
 * 
 * *******************************************************************************/

#include "Morfo70.h"
#include <dirent.h>
#include <vector>



namespace Morfo70{

     //Offshore wave conditions
     double H_off=0.0001;
     double T_off=6;
     double Theta_off=0;
     double Zs_off=0; 

     string inputDir="";
     string outputDir="";
     list<string> iniFiles=list<string>();
     string batFile="";
     

     bool readAndCheckParameters( int argc, const char** argv )
     {
     ///Read input parameters
     for (int i=1;i<argc;i++)
     {
        string arg=argv[i];        
        
        if (arg.compare(0,5,"-ind=")==0)
            inputDir=arg.substr(5);
        if (arg.compare(0,6,"-outd=")==0)
            outputDir=arg.substr(6);
        if (arg.compare(0,6,"-inif=")==0)
        {
            string iniFiles_=arg.substr(6);
            int pos=iniFiles_.find(',');            
            while (pos>=0)                
            {
                iniFiles.push_back(inputDir+ PATH_SEPARATOR+iniFiles_.substr(0,pos));
                iniFiles_=iniFiles_.substr(pos+1);  
                pos=iniFiles_.find(',');          
            }
            iniFiles.push_back(inputDir+ PATH_SEPARATOR+iniFiles_);                           
            

            
        }
        if (arg.compare(0,6,"-batf=")==0)
            batFile=inputDir+ PATH_SEPARATOR+arg.substr(6);
    }

    ///Check Input
    if (!checkDir(outputDir))
        mkdir(outputDir.c_str(),0777);
    for (string iniFile:iniFiles)
    {
       if (!checkFile(iniFile))    {
         cout<<"Ini File "+iniFile+ " No exist"<<"\n"; 
         return false;
       }
    }
    
    if (!checkFile(batFile))  
    {  
        cout<<"Bat File "+batFile+ " No exist"<<"\n";   

        return false;
        }

        return true;

    
        
}

   
     void Run(Grid* grid,Bathymetry * bathy,Wave* wave, Hydro* hydro,Sediment* sediment, string outputDir){
          
          double time=0;
          double warmTime=Parameters::Run::warmTime;                    
          double stopTime=Parameters::Run::stopTime;     
          double bedStartTime=Parameters::Run::bedStartTime;
          double bedStopTime=Parameters::Run::bedStopTime;
          double saveInterval=Parameters::Run::saveIntervalHydro;     
          double dt=Parameters::Run::dt;  
          double H_off0=0.01; //Oleaje minimo 
          if (warmTime<=0)           
               H_off0=H_off;

          double lastSavedTime=-1;
          
          if (Options::Run::doResume)
          {
               double timeResume=resume(outputDir,bathy,wave,hydro,sediment);
               if (timeResume>0) // There is resume!!!
               {
                    time=timeResume;
                    lastSavedTime=floor(time);
               }
          }
          

          //Auxiliary vars to get for each loop the time-step variation
          // m1 -> Previous time step
          const size_t size=grid->size;          
          double* U_m1=new double[size];
          double* V_m1=new double[size];
          double* Zs_m1=new double[size];
          double* E_m1=new double[size];
          double* Er_m1=new double[size];
          double* h_m1=new double[size];
          std::copy(wave->E,wave->E+size,E_m1);
          std::copy(wave->E_r,wave->E_r+size,Er_m1);
          std::copy(hydro->U_c,hydro->U_c+size,U_m1);          
          std::copy(hydro->V_c,hydro->V_c+size,V_m1);
          std::copy(hydro->Zs,hydro->Zs+size,Zs_m1);
          std::copy(bathy->h,bathy->h+size,h_m1);
          
          
          //MAIN LOOP!!!!
          while (time<stopTime)
          {
               //Warming function, increasis slightly wave height
               double waveFactor=(1+tanh((5/warmTime)*(time -warmTime/2)))/2;  
               if (warmTime<=0)   
                    waveFactor=1;
               
               //offshore wave height
               double H_off_t=max(H_off0,waveFactor*H_off);               
               
               wave->calculateWaveField(H_off_t,Theta_off,T_off);

               bool updateBottom=(time>=bedStartTime && time<bedStopTime) && Options::Bottom::CalulateBedEvolution;
               int hydroIterations=Parameters::Hydro::runIterations;
               if (updateBottom) 
                    hydroIterations=Parameters::Hydro::runIterationsInBedEvolution;
               
               if (Options::Hydro::CalculateHydro)
               {
                     hydro->Run(hydroIterations);
                     if (Options::Wave::WaveCurrentInteraction)
                     {
                          wave->setUV(hydro->U_c,hydro->V_c);
                     }                  
               }

               if (Options::Sediment::CalculateSedimentTransport)
               {
                    sediment->calculateSedimentFluxes();
               }
               
               if (updateBottom)
               {
                    sediment->updateBottom();
                    saveInterval=Parameters::Run::saveIntervalMorfo;
               }
               else               
                    saveInterval=Parameters::Run::saveIntervalHydro;
               
               
                    

               time=time+dt;

               ///////////////////////////////////////////////////////////////////////////////
               ///Print max vars variations in shell
               ///////////////////////////////////////////////////////////////////////////////
               double errE=0,errEr=0,errU=0,errV=0,errZs=0,errh=0;
               const int nx=grid->nx;
               const int ny=grid->ny;
               
               for (int i=0;i<nx;i++)
               {
                    for (int j=0;j<ny;j++)
                    {
                         errE=max(errE,abs(E_m1[i*ny+j]-wave->E[i*ny+j]));
                         errEr=max(errEr,abs(Er_m1[i*ny+j]-wave->E_r[i*ny+j]));
                         errU=max(errU,abs(U_m1[i*ny+j]-hydro->U_c[i*ny+j]));                         
                         errV=max(errV,abs(V_m1[i*ny+j]-hydro->V_c[i*ny+j]));
                         errZs=max(errZs,abs(Zs_m1[i*ny+j]-hydro->Zs[i*ny+j]));
                         errh=max(errh,abs(h_m1[i*ny+j]-bathy->h[i*ny+j]));                         
                    }
               }                              
               printf("t: %.2f s;h0: %.2f m;dE: %.7f;dR: %.5f;dU: %.7f;dV: %.7f;dZs: %.7f;dh: %.10f \n",time,H_off_t,errE,errEr,errU,errV,errZs,errh);
               std::copy(wave->E,wave->E+size,E_m1);
               std::copy(wave->E_r,wave->E_r+size,Er_m1);
               std::copy(hydro->U_c,hydro->U_c+size,U_m1);
               std::copy(hydro->V_c,hydro->V_c+size,V_m1);
               std::copy(hydro->Zs,hydro->Zs+size,Zs_m1);
               std::copy(bathy->h,bathy->h+size,h_m1);
               ///////////////////////////////////////////////////////////////////////////////
               ///END Print results
               ///////////////////////////////////////////////////////////////////////////////

               ///Save time step
               if (saveInterval>0)
                    if (saveInterval*floor(time/saveInterval)>lastSavedTime)
                    {                         
                         string outputFile="000000000"+to_string((int)round(time));
                         outputFile=outputFile.substr(outputFile.length()-9);
                         outputFile=outputDir+PATH_SEPARATOR+"m70_"+outputFile+".nc";
                         saveTimeStep(outputFile,time,grid,bathy,wave,hydro,sediment);
                         lastSavedTime=saveInterval*floor(time/saveInterval);    
                                        
                    }              

               
          }
          //free auxiliary variables
          delete[] E_m1;delete[] Er_m1;delete[] U_m1;delete[] V_m1;delete[] Zs_m1;delete[] h_m1;          
          

     }

     void saveTimeStep(string outputFile,double time,Grid* grid,Bathymetry* bathy,  Wave* wave, Hydro* hydro,Sediment* sediment){
          
          MorfoIO::Netcdf nc=MorfoIO::Netcdf();
          nc.Open(outputFile,MorfoIO::Netcdf::OpenMode::Write);
          nc.WriteXY(grid->x,grid->y,grid->nx,grid->ny);
          nc.WriteVar("Z",bathy->Z); 
          nc.WriteVar("h",bathy->h); 
          nc.WriteVar("D",bathy->D); 
          nc.WriteVar("H",wave->H); 
          nc.WriteVar("E",wave->E); 
          // nc.WriteVar("Sxx",wave->Sxx); 
          // nc.WriteVar("Syy",wave->Syy); 
          // nc.WriteVar("Sxy",wave->Sxy); 
          // nc.WriteVar("Urms",wave->Urms); 
          nc.WriteVar("Cg",wave->Cg); 
          nc.WriteVar("lK",wave->K); 
          // nc.WriteVar("f",hydro->f);
          // nc.WriteVar("d",hydro->d);
          // nc.WriteVar("up",hydro->up);
          // nc.WriteVar("low",hydro->low);

          nc.WriteVar("dQx_dx",hydro->dQx_dx);
          nc.WriteVar("dQy_dy",hydro->dQy_dy);
          
          nc.WriteVar("Theta",wave->Theta); 
          nc.WriteVar("Dw",wave->Dw); 
          nc.WriteVar("E_r",wave->E_r); 
          nc.WriteVar("U",hydro->U_c); 
          nc.WriteVarY("U0",hydro->U); 
          nc.WriteVarY("Qx0",hydro->Qx); 
          nc.WriteVar("U_x",hydro->U+grid->ny); 
          nc.WriteVar("Qx_x",hydro->Qx+grid->ny); 
          nc.WriteVar("V",hydro->V_c); 
          nc.WriteVarY("V0",hydro->V); 
          nc.WriteVarY("Qy0",hydro->Qy); 
          nc.WriteVar("V_y",hydro->V+grid->ny); 
          nc.WriteVar("Qy_y",hydro->Qy+grid->ny); 
          nc.WriteVar("Zs",hydro->Zs); 
          
          nc.WriteVar("qx",sediment->qx); 
          nc.WriteVar("qy",sediment->qy); 
          //nc.WriteVar("dh_dt",sediment->dh_dt); 

          nc.WriteVar("Gamma_Swash",sediment->Gamma); 
          nc.WriteVar("Alpha",sediment->Alpha); 
          nc.WriteVar("gamma",sediment->gamma); 


          nc.WriteAtt("time",time);
          nc.Close();
    
     }


     double resume(string outputDir,Bathymetry* bathy,  Wave* wave, Hydro* hydro,Sediment* sediment)
     {
       
       string resumeFile="";

       /// look in output folder the last m70_*.nc file
       int lastFileTime=-1;
       DIR* dirp = opendir(outputDir.c_str());
       struct dirent * dp;
        while ((dp = readdir(dirp)) != NULL) {
             string fileName=dp->d_name;
             int com_len=MorfoIO::compare(fileName,"m70_");
             if (com_len)
             {
                  int fileTime=stoi(fileName.substr(com_len,9));
                  if (fileTime>lastFileTime)
                  {
                    resumeFile=fileName;
                    lastFileTime=fileTime;
                  }
                  
             }
          }
          closedir(dirp);

          if (lastFileTime>-1)  // There is resume file
          {
               resumeFile=outputDir+PATH_SEPARATOR+resumeFile;
               bathy->setBathyFromResumeFile(resumeFile);
               wave->setWaveEnergyFromResumeFile(resumeFile);
               if (Options::Wave::WaveCurrentInteraction)
                    wave->setUVFromResumeFile(resumeFile);
               hydro->setHydroFromResumeFile(resumeFile);
               sediment->setSedimentFluxesFromResumeFile(resumeFile);

               MorfoIO::Netcdf ncResume=MorfoIO::Netcdf();
               ncResume.Open(resumeFile,MorfoIO::Netcdf::OpenMode::Read);
               double time=0;
               ncResume.ReadAtt("time",&time);
               ncResume.Close();
               return time;

          }
          else                //No  resume          
               return -1;               
          
          
     }

}