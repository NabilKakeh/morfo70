#include "MorfoIniFile.h"

using namespace std;


namespace MorfoIO{

void loadIniFile(string iniFile){

    fstream infile(iniFile);    
    string line;
    
    while(std::getline(infile,line)){
        string subString=trim(line); 

        readOffshoreData(subString);


        int com_len=compare(subString,"options.wave.");
        if (com_len)    
                readWaveOption(subString.substr(com_len));
        
        com_len=compare(subString,"parameters.wave.");
        if (com_len)    
                readWaveParameter(subString.substr(com_len));

        com_len=compare(subString,"options.rollers.");
        if (com_len)    
                readRollerOption(subString.substr(com_len));

        com_len=compare(subString,"parameters.rollers.");
        if (com_len)    
                readRollerParameter(subString.substr(com_len));

        com_len=compare(subString,"options.hydro.");
        if (com_len)    
                readHydroOption(subString.substr(com_len));

        com_len=compare(subString,"parameters.hydro.");
        if (com_len)    
                readHydroParameter(subString.substr(com_len));
        com_len=compare(subString,"options.bottom.");
        if (com_len)    
                readBottomOption(subString.substr(com_len));                
        com_len=compare(subString,"parameters.bottom.");
        if (com_len)    
                readBottomParameter(subString.substr(com_len)); 
        
        com_len=compare(subString,"parameters.sediment.CWS.");
        if (com_len)    
                readSedimentCWSParameter(subString.substr(com_len));     

        com_len=compare(subString,"parameters.sediment.Swash.");
        if (com_len)    
                readSedimentSwashParameter(subString.substr(com_len));     

            
        com_len=compare(subString,"options.run.");
        if (com_len)    
                readRunOption(subString.substr(com_len));
        com_len=compare(subString,"parameters.run.");
        if (com_len)    
                readRunParameter(subString.substr(com_len));
    }

    
}

void readWaveOption(string option){
    int com_len=compare(option,"BreakingModel=");
    if (com_len)
    {
        string breakModel=option.substr(com_len);
        if (compare(breakModel,"None"))        
            Options::Wave::BreakingModel=Options::Wave::BreakingModelType::None;
        if (compare(breakModel,"ChurchThornton"))        
            Options::Wave::BreakingModel=Options::Wave::BreakingModelType::ChurchThornton;
        if (compare(breakModel,"ThorntonGuza"))        
            Options::Wave::BreakingModel=Options::Wave::BreakingModelType::ThorntonGuza;
        if (compare(breakModel,"RegularWaves"))        
            Options::Wave::BreakingModel=Options::Wave::BreakingModelType::RegularWaves;
    }
    com_len=compare(option,"BottomFriction=");
    if (com_len)    
        Options::Wave::BottomFriction=stoi(option.substr(com_len));
    
    com_len=compare(option,"FrictionCoefficientModel=");
    if (com_len)
    {
        string fricModel=option.substr(com_len);
        if (compare(fricModel,"Constant"))        
            Options::Wave::fwModel=Options::Wave::FricCoefficientModelType::Constant;
            //{}//Options::Wave::FrictionCoefficientModel=Options::Wave::FrictionCoefficentModelType::Constant;
        if (compare(fricModel,"Souslby97"))        
            Options::Wave::fwModel=Options::Wave::FricCoefficientModelType::Soulsby97;
        //Options::Rollers::enumTipo=Options::Rollers::TipoEnum::Opcion1;
       
    }    
    
    com_len=compare(option,"WaveCurrentInteraction=");
    if (com_len)    
        Options::Wave::WaveCurrentInteraction=stoi(option.substr(com_len));
    

}

void readWaveParameter(string parameter)
{
    int com_len=compare(parameter,"B3=");
    if (com_len)
        Parameters::Wave::B3=stod(parameter.substr(com_len));
    
    com_len=compare(parameter,"gamma_b=");
    if (com_len)
        Parameters::Wave::gamma_b=stod(parameter.substr(com_len));
    
    com_len=compare(parameter,"n_mon=");
    if (com_len)
        Parameters::Wave::n_mon=stod(parameter.substr(com_len));
    
    com_len=compare(parameter,"m_mon=");
    if (com_len)
        Parameters::Wave::m_mon=stod(parameter.substr(com_len));
    
    com_len=compare(parameter,"ThetaMax=");
    if (com_len)
        Parameters::Wave::ThetaMax=stod(parameter.substr(com_len))*pi/180; //lee en grados
    
    com_len=compare(parameter,"SolverMaxIterations=");
    if (com_len)
        Parameters::Wave::SolverMaxIterations=stoi(parameter.substr(com_len));
    
    com_len=compare(parameter,"SolverTolerance=");
    if (com_len)
        Parameters::Wave::SolverTolerance=stod(parameter.substr(com_len));
    
    com_len=compare(parameter,"fw0=");
    if (com_len)
        Parameters::Wave::fw0=stod(parameter.substr(com_len));
}

void readRollerOption(string option){   
    int com_len=compare(option,"CalculateRollers=");
    if (com_len)    
        Options::Rollers::CalculateRollers=stoi(option.substr(com_len));

}

void readRollerParameter(string parameter)
{
    int com_len=compare(parameter,"sinBeta_r=");
    if (com_len)
        Parameters::Rollers::sinBeta_r=stod(parameter.substr(com_len));
    com_len=compare(parameter,"alpha_r=");
    if (com_len)
        Parameters::Rollers::alpha_r=stod(parameter.substr(com_len));
    
}

void readHydroOption(string option)
{
    int com_len=compare(option,"DragModel=");
    if (com_len)
    {
        string dragModel=option.substr(com_len);
        if (compare(dragModel,"Constant"))        
            Options::Hydro::DragModel=Options::Hydro::DragModelType::ConstantDrag;
        if (compare(dragModel,"LogDepth"))        
            Options::Hydro::DragModel=Options::Hydro::DragModelType::LogDepth;
        if (compare(dragModel,"ManningStrickler"))        
            Options::Hydro::DragModel=Options::Hydro::DragModelType::ManningStrickler;        
        
    }
   com_len=compare(option,"ViscosityModel=");
    if (com_len)
    {
        string viscModel=option.substr(com_len);
        if (compare(viscModel,"Constant"))        
            Options::Hydro::ViscosityModel=Options::Hydro::ViscosityModelType::ConstantViscosity;
        if (compare(viscModel,"BreakingRollers"))        
            Options::Hydro::ViscosityModel=Options::Hydro::ViscosityModelType::BreakingRollers;
        if (compare(viscModel,"BreakingDepth"))        
            Options::Hydro::ViscosityModel=Options::Hydro::ViscosityModelType::BreakingDepth;
        if (compare(viscModel,"BreakingHrms"))        
            Options::Hydro::ViscosityModel=Options::Hydro::ViscosityModelType::BreakingHrms;        
        
    }
    com_len=compare(option,"OffshoreBC=");
    if (com_len)
    {
        string offBC=option.substr(com_len);
        if (compare(offBC,"Sponge"))        
            Options::Hydro::OffshoreBC=Options::Hydro::OffshoreBCType::Sponge;
        if (compare(offBC,"ConstantZs"))        
            Options::Hydro::OffshoreBC=Options::Hydro::OffshoreBCType::ConstantFreeSurface;
    }

    com_len=compare(option,"CalculateHydro=");
    if (com_len)    
        Options::Hydro::CalculateHydro=stoi(option.substr(com_len));

}

void readHydroParameter(string parameter){    
     int com_len=compare(parameter,"Dmin=");
     if (com_len)
        Parameters::Hydro::Dmin=stod(parameter.substr(com_len));
    com_len=compare(parameter,"alpha_f=");
    if (com_len)
        Parameters::Hydro::alpha_f=stod(parameter.substr(com_len));
    com_len=compare(parameter,"cD0=");
    if (com_len)
        Parameters::Hydro::cD0=stod(parameter.substr(com_len));
    com_len=compare(parameter,"nu0=");
    if (com_len)
        Parameters::Hydro::nu0=stod(parameter.substr(com_len));
    com_len=compare(parameter,"M=");
    if (com_len)
        Parameters::Hydro::M=stod(parameter.substr(com_len));    
    com_len=compare(parameter,"ka=");
    if (com_len)
        Parameters::Hydro::ka=stod(parameter.substr(com_len));
    com_len=compare(parameter,"dt=");    
    if (com_len)
        Parameters::Hydro::dt=stod(parameter.substr(com_len));

    com_len=compare(parameter,"runIterations=");
    if (com_len)
        Parameters::Hydro::runIterations=stoi(parameter.substr(com_len));

    com_len=compare(parameter,"runIterationsInBedEvolution=");
    if (com_len)
        Parameters::Hydro::runIterationsInBedEvolution=stoi(parameter.substr(com_len));

}

void readBottomParameter(string parameter){
    int com_len=compare(parameter,"z0=");
    if (com_len)
         Parameters::Bottom::z0=stod(parameter.substr(com_len));
    
    com_len=compare(parameter,"d50=");
    if (com_len)
         Parameters::Bottom::d50=stod(parameter.substr(com_len));    

    com_len=compare(parameter,"p=");
    if (com_len)
         Parameters::Bottom::p=stod(parameter.substr(com_len));    
}

void readBottomOption(string option){
    int com_len=compare(option,"CalculateBedEvolution=");
    if (com_len)    
        Options::Bottom::CalulateBedEvolution=stoi(option.substr(com_len));
}

void readSedimentCWSParameter(string parameter){
    int com_len=compare(parameter,"alpha0_max=");
    if (com_len)
         Parameters::Sediment::CWS::alpha0_max=stod(parameter.substr(com_len));

    com_len=compare(parameter,"alpha0_min=");
    if (com_len)
         Parameters::Sediment::CWS::alpha0_min=stod(parameter.substr(com_len));

    com_len=compare(parameter,"D_alpha0_max=");
    if (com_len)
         Parameters::Sediment::CWS::D_alpha0_max=stod(parameter.substr(com_len));

    com_len=compare(parameter,"D_alpha0_min=");
    if (com_len)
         Parameters::Sediment::CWS::D_alpha0_min=stod(parameter.substr(com_len));

    com_len=compare(parameter,"gamma0=");
    if (com_len)
         Parameters::Sediment::CWS::gamma0=stod(parameter.substr(com_len));
    
    
}

void readSedimentSwashParameter(string parameter){
   

    int com_len=compare(parameter,"Gamma0=");
    if (com_len)
         Parameters::Sediment::Swash::Gamma0=stod(parameter.substr(com_len));

    com_len=compare(parameter,"A_swash=");
    if (com_len)
         Parameters::Sediment::Swash::A_swash=stod(parameter.substr(com_len));
    
}

void readSedimentOption(string option){
    int com_len=compare(option,"TransportModel=");
    if (com_len)
    {
        string transportModel=option.substr(com_len);
        if (compare(transportModel,"CWS"))        
           Options::Sediment::TransportModel=Options::Sediment::TransportModelType::CWS;                
    }

    com_len=compare(option,"CalculateSedimentTransport=");
    if (com_len)    
        Options::Sediment::CalculateSedimentTransport=stoi(option.substr(com_len));

    com_len=compare(option,"CrosshoreTransport=");
    if (com_len)    
        Options::Sediment::CrosshoreTransport=stoi(option.substr(com_len));

 }

void readRunOption(string option){
    int com_len=compare(option,"doResume=");
    if (com_len)
         Options::Run::doResume=stoi(option.substr(com_len));  
}

void readRunParameter(string parameter){
    int com_len=compare(parameter,"warmTime=");
    if (com_len)
         Parameters::Run::warmTime=stod(parameter.substr(com_len));
    
    com_len=compare(parameter,"stopTime=");
    if (com_len)
         Parameters::Run::stopTime=stod(parameter.substr(com_len));  
    
    com_len=compare(parameter,"bedStartTime=");
    if (com_len)
         Parameters::Run::bedStartTime=stod(parameter.substr(com_len));  

    com_len=compare(parameter,"bedStopTime=");
    if (com_len)
         Parameters::Run::bedStopTime=stod(parameter.substr(com_len));  
    
    com_len=compare(parameter,"saveIntervalHydro=");
    if (com_len)
         Parameters::Run::saveIntervalHydro=stod(parameter.substr(com_len));    
    
    com_len=compare(parameter,"saveIntervalMorfo=");
    if (com_len)
         Parameters::Run::saveIntervalMorfo=stod(parameter.substr(com_len));    
    
    com_len=compare(parameter,"dt=");
    if (com_len)
         Parameters::Run::dt=stod(parameter.substr(com_len));    
    
      
}

void readOffshoreData(string line){
    int com_len=compare(line,"H_off=");
    if (com_len)
         Morfo70::H_off=stod(line.substr(com_len));

    com_len=compare(line,"T_off=");
    if (com_len)
         Morfo70::T_off=stod(line.substr(com_len));
    
    com_len=compare(line,"Theta_off="); //Lee en grados
    if (com_len)
         Morfo70::Theta_off=stod(line.substr(com_len))*pi/180;

    com_len=compare(line,"Zs_off="); 
    if (com_len)
         Morfo70::Zs_off=stod(line.substr(com_len));
}

}


 
    
