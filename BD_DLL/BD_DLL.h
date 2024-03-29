#pragma once

#include <msclr/marshal_cppstd.h>
#include "../../Core/SimulationParameters.h"
#include "../../Core/discrete/Simulation3D.h"

using namespace System::Collections::Generic;

namespace BDDLL {

	public ref class ManagedSimulationParameters
	{
	public:
		ManagedSimulationParameters() {
			nativeParameters = new SimulationParametersFull();
		}

		~ManagedSimulationParameters() {
			delete nativeParameters;
		}

		!ManagedSimulationParameters() {
			delete nativeParameters;
		}

		SimulationParametersFull* GetNativePointer() {
			return nativeParameters;
		}

	private:
		SimulationParametersFull* nativeParameters;
	};


	public ref class BallisticSimulation
	{
	public:
		BallisticSimulation();

		int DiscreteSimulation(float theta, 
			uint16_t L, 
			uint16_t H, 
			uint32_t reps, 
			float phi, 
			float turns, 
			uint32_t seed,
			uint16_t diffusion_steps, 
			List<System::Byte>^ species,
			List<List<float>^>^ weights,
			List<float>^ spread, 
			array<int16_t>^ inputGrid, 
			uint32_t inputGridPoints, 
			array<int16_t>^% outGrid, 
			int phi_num, 
			float phi_deg,
			bool thetaSweep,
			float thetaEnd,
			uint32_t stepper_resolution,
			ManagedSimulationParameters^ params,
			System::String^ system,  
			Acceleration acc, 
			Collision collision_method, 
			System::String^ directory);
	};
}
