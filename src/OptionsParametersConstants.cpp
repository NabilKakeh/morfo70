/*********************************************************************************
 * 
 *                  OptionsParametersConstants.h (Morfo70)
 * 
 *  Definiton of global constants, Options and Parameters  
 * 
 * *******************************************************************************/

#include "OptionsParametersConstants.h"
namespace Options{
    namespace Bottom{
		bool CalulateBedEvolution=false;
	}
    namespace Wave{
        BreakingModelType BreakingModel=BreakingModelType::ChurchThornton;
        FricCoefficientModelType fwModel=FricCoefficientModelType::Constant;
        bool BottomFriction=false;
        bool WaveCurrentInteraction=false;
    }
    namespace Rollers{
        bool CalculateRollers=false;
        //TipoEnum enumTipo=TipoEnum::Opcion1;
    }
    namespace  Hydro{
        DragModelType DragModel=DragModelType::LogDepth;
        ViscosityModelType ViscosityModel=ViscosityModelType::BreakingDepth;
        OffshoreBCType OffshoreBC=OffshoreBCType::Sponge;
        bool CalculateHydro=true;
    }
    namespace Sediment{	
		TransportModelType TransportModel=TransportModelType::CWS;
		bool CalculateSedimentTransport=true;
		bool CrosshoreTransport=false;
	}
    namespace Run{
        bool doResume=true;
    }
}

namespace Parameters {
    namespace Bottom{
        double z0=0.001;
        double d50=215E-6;
        double p=0.4;
    }
    namespace Wave{
        double B3=2.0;
        double gamma_b=0.475;
        double m_mon=20;
        double n_mon=30;
        double ThetaMax=85;
        int SolverMaxIterations=10;
        double SolverTolerance=1E-10;
        double fw0=0.01;
    }
    namespace Rollers{
        double sinBeta_r=0.05;
        double alpha_r=1;
    }
    namespace Hydro{
        double Dmin=0.01;
        double alpha_f=1.16;
        double cD0=0.01; //Constant drag
		double nu0=0.0001; //Constant drag
		double M=1; //Viscosity parameter
		//double z0=0.006; //Rugosidad de fondo para calcular drag (m)        
		double ka=0.0219; //Rugosidad aparente (en teoria ka=30*z0) para Manning-Strickler
        double dt=0.1;
        int runIterations=1;
        int runIterationsInBedEvolution=1;
        
    }

    namespace Sediment{
		namespace CWS{
		double alpha0_max=0.002;
		double alpha0_min=0.0;
		double D_alpha0_max=4;
        double D_alpha0_min=10;
		double gamma0=0.5;
		
		}
        namespace Swash{
            double Gamma0=1;
            double A_swash=0.3;
        }
	}
    
    namespace Run{
        double warmTime=3600;
        double stopTime=7200;
        double bedStartTime=3600;
		double bedStopTime=7200;
        double saveIntervalHydro=600;
		double saveIntervalMorfo=3600;
        double dt=0.1;
        
    }
    

}