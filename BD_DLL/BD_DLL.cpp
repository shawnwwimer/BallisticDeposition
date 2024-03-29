#include "pch.h"

#include "BD_DLL.h"

#include "../Core/discrete/Simulation3D.h"
#include "../Core/cts/SimulationContinuous.h"


BDDLL::BallisticSimulation::BallisticSimulation() {}

int BDDLL::BallisticSimulation::DiscreteSimulation(float theta,
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
	System::String^ directory)
{
	// Convert to native types
	std::string nativeSystem = msclr::interop::marshal_as<std::string>(system);
	std::string nativeDirectory = msclr::interop::marshal_as<std::string>(directory);
	SimulationParametersFull* nativeParameters = params->GetNativePointer();
	
	// Fill native type containers
	std::vector<int8_t> nativeSpecies;
	for each (System::Byte byteValue in species) {
		nativeSpecies.push_back(static_cast<int8_t>(byteValue));
	}
	std::vector<float> nativeSpread;
	for each (float floatValue in spread) {
		nativeSpread.push_back(floatValue);
	}
	std::vector<std::vector<float>> nativeWeights;
	for each (List<float> ^ innerList in weights) {
		std::vector<float> nativeInnerVector;
		for each (float floatValue in innerList) {
			nativeInnerVector.push_back(floatValue);
		}
		nativeWeights.push_back(nativeInnerVector);
	}

	pin_ptr<int16_t> pinnedInputGrid = &inputGrid[0];
	int16_t* nativeOutGrid = new int16_t[inputGridPoints + reps];
	
	int points;
	if (phi_num == 0) {
		points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, diffusion_steps, &nativeSpecies, &nativeSpread, &nativeWeights, pinnedInputGrid, inputGridPoints, &nativeOutGrid, 0, 0, stepper_resolution, nativeParameters, nativeSystem, false, thetaSweep, thetaEnd, acc, collision_method, nativeDirectory);
	}
	else {
		points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, diffusion_steps, &nativeSpecies, &nativeSpread, &nativeWeights, pinnedInputGrid, inputGridPoints, &nativeOutGrid, phi_num, phi_deg, stepper_resolution, nativeParameters, nativeSystem, true, thetaSweep, thetaEnd, acc, collision_method, nativeDirectory);
	}

	outGrid = gcnew array<int16_t>(inputGridPoints + reps);
	for (size_t i = 0; i < inputGridPoints; ++i) {
		outGrid[i] = nativeOutGrid[i];
	}
	delete[] nativeOutGrid;

	return 0;
}