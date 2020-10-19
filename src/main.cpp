
#include "Morfo70.h"


using namespace std;


int main( int argc, const char* argv[] )
{    

    if  (!Morfo70::readAndCheckParameters(argc,argv))
        return 0;      

    for (string iniFile:Morfo70::iniFiles)    
        MorfoIO::loadIniFile(iniFile);       

    Grid* grid= new Grid(Morfo70::batFile);
    Bathymetry* bathy=new Bathymetry(Morfo70::batFile,Morfo70::Zs_off,grid);
    Wave* wave=new Wave(grid,bathy);
    Hydro* hydro=new Hydro(grid,bathy,wave);
    Sediment* sediment=new Sediment(hydro,wave,bathy,grid);
    
    Morfo70::Run(grid,bathy,wave,hydro,sediment,Morfo70::outputDir);
   
    
    delete(hydro);
    delete(wave);
    delete(bathy);
    delete(grid);
       
      return 0;
     
}