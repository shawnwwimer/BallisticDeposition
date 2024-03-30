#include "pch.h"

#include "BD_DLL.h"

#include "../Core/discrete/Simulation3D.h"
#include "../Core/cts/SimulationContinuous.h"


BDDLL::BallisticSimulation::BallisticSimulation() {}

int BDDLL::BallisticSimulation::DiscreteSimulation(
	float theta,
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

int BDDLL::BallisticSimulation::initializeSimulation(uint16_t L, uint16_t W, uint16_t H) {
	inGrid = (int16_t*)malloc(sizeof(int16_t) * L * L * 4);
	if (inGrid == nullptr) {
		return -1;
	}
	outGrid = (int16_t**)malloc(sizeof(int16_t*));
	if (outGrid == nullptr) {
		free(inGrid);
		return -1;
	}
	*outGrid = nullptr;
	params = new SimulationParametersFull();

	uint32_t count = 0;
	for (uint32_t i = 0; i < L; i++) {
		for (uint32_t j = 0; j < W; j++) {
			inGrid[count * 4] = i;
			inGrid[count * 4 + 1] = j;
			inGrid[count * 4 + 2] = 0;
			inGrid[count * 4 + 3] = 1;
			count += 1;
		}
	}
	inGridPoints = count;
	initialized = true;
	return 0;
}

void BDDLL::BallisticSimulation::cleanup() {
	free(*outGrid);
	*outGrid = nullptr;
	free(outGrid);
	outGrid = nullptr;
	free(inGrid);
	inGrid = nullptr;
	inGridPoints = 0;
	delete params;
	params = nullptr;
	initialized = false;
	std::cout << "Cleaned up objects." << std::endl;
}

int BDDLL::BallisticSimulation::simulateFilm(
	float theta,
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
	int phi_num,
	float phi_deg,
	bool thetaSweep,
	float thetaEnd,
	uint32_t stepper_resolution,
	System::String^ system,
	Acceleration acc,
	Collision collision_method,
	System::String^ directory
) {
	int ret = this->initializeSimulation(L, L, H);
	if (ret != 0) {
		this->cleanup();
		return -1;
	}
	this->simulating_now = true;

	// Convert to native types
	std::string nativeSystem = msclr::interop::marshal_as<std::string>(system);
	std::string nativeDirectory = msclr::interop::marshal_as<std::string>(directory);
	//SimulationParametersFull* nativeParameters = params->GetNativePointer();

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

	uint32_t points = 0;
	if (phi_num == 0) {
		points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, diffusion_steps, &nativeSpecies, &nativeSpread, &nativeWeights, this->inGrid, inGridPoints, this->outGrid, 0, 0, stepper_resolution, this->params, nativeSystem, false, thetaSweep, thetaEnd, acc, collision_method, nativeDirectory);
	}
	else {
		points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, diffusion_steps, &nativeSpecies, &nativeSpread, &nativeWeights, this->inGrid, inGridPoints, this->outGrid, phi_num, phi_deg, stepper_resolution, this->params, nativeSystem, true, thetaSweep, thetaEnd, acc, collision_method, nativeDirectory);
	}
	inGridPoints += points;

	this->cleanup();
	this->simulating_now = false;
	return points;
}