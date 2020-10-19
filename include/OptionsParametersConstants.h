/*********************************************************************************
 * 
 *                  OptionsParametersConstants.h (Morfo70)
 * 
 *  Definiton of global constants, declaration of Options and Parameters 
 *  Defintion in :
 *      - OptionsParametersConstants.cpp
 * 
 * *******************************************************************************/
#ifndef OPTIONSPARAMETERSCONSTANTS_H_
#define OPTIONSPARAMETERSCONSTANTS_H_


double const g=9.8100000000000000;
double const pi=3.141592653589793;
double const rho=1025.00000000000;
double const s=2.65;
double const sqrt_pi=1.772453850905516;
double const sqrt_2=1.414213562373095;
double const cteEnergy=0.125*rho*g;
double const inv_2pi=0.5/pi;
double const tanPhi_i=0.6249; //tan(32)
double const inv_tanPhi_i=1/tanPhi_i;


namespace Options{
	namespace Bottom{
		extern bool CalulateBedEvolution;	
	}
    namespace Wave{
        enum BreakingModelType{
            None=0,
            ChurchThornton=1,
            ThorntonGuza=2,
			RegularWaves=3
        };
		enum FricCoefficientModelType{
			Constant,
			Soulsby97
		};
		extern FricCoefficientModelType fwModel;		
		extern BreakingModelType BreakingModel;		
        extern bool BottomFriction;
		extern bool WaveCurrentInteraction;
    }
    namespace Rollers{		
        extern bool CalculateRollers;
    }
	namespace Hydro{
		enum DragModelType{
			ConstantDrag=0,
			LogDepth=1,
			ManningStrickler=2
		};
		enum ViscosityModelType{
			ConstantViscosity=0,
			BreakingRollers=1,
			BreakingDepth=2,
			BreakingHrms=3
		};
		enum OffshoreBCType{
			Sponge=0,
			ConstantFreeSurface=1
		};
		extern DragModelType DragModel;
		extern ViscosityModelType ViscosityModel;		
		extern OffshoreBCType OffshoreBC;
		extern bool CalculateHydro;
	}
	namespace Sediment{
		 enum TransportModelType{
            CWS            
        };
		extern TransportModelType TransportModel;
		extern bool CalculateSedimentTransport;
		extern bool CrosshoreTransport;
	}
	namespace Run{        
		extern bool doResume;
    }
}

namespace Parameters{
	namespace Bottom {
		extern double z0; //rugosidad de fondo para calcular drag
		extern double d50; //grain size
		extern double p;  //porosity (~ 0.4)
	}
	namespace Wave {
		extern double B3;
		extern double gamma_b;
		extern double n_mon; //parameter for breaking dissipation for regular/monochomatric waves
		extern double m_mon; //parameter for breaking dissipation for regular/monochomatric waves
		extern double ThetaMax;
		extern int SolverMaxIterations;
		extern double SolverTolerance;
		extern double fw0; 
	}
	namespace Rollers{
		extern double sinBeta_r;
		extern double alpha_r;
	}
	namespace Hydro {
		extern double Dmin;
		extern double alpha_f; //Federsen
		extern double cD0; //Constant drag
		extern double nu0; //Constant drag
		extern double M; //Viscosity parameter		
		extern double ka; //Rugosidad aparente (en teoris ka=30*z0) para Manning-Strickler
		extern double dt; //paso de tiempo
		extern int runIterations; //numero de iteraciones hydro
		extern int runIterationsInBedEvolution; //numero de iteraciones hydro cuando hay cambio de fondo
	}
	
	
	namespace Run{
		extern double warmTime;
		extern double stopTime;
		extern double bedStartTime;
		extern double bedStopTime;
		extern double saveIntervalHydro;
		extern double saveIntervalMorfo;		
		extern double dt;
	
	}
	namespace Sediment{
		namespace CWS{
		extern double alpha0_max;
		extern double alpha0_min;
		extern double D_alpha0_max;
		extern double D_alpha0_min;
		extern double gamma0;		
		}
		namespace Swash{
			extern double Gamma0;
			extern double A_swash;
		}
	}
	
}
#endif