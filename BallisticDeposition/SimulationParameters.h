#pragma once
#include <cstdint>
#include <vector>
#include <string>

class SimulationParametersIndividual
{
public:
	// Physical parameters
	uint16_t length;
	uint16_t width;
	uint16_t height;
	float theta;
	float phi;
	float turns;
	std::vector<float> * spread;
	uint32_t stepper_resolution;
	
	// Material parameters
	std::vector<int8_t> * species;
	uint16_t diffusion_steps;
	uint32_t deposited;
	uint32_t repetitions;
	std::vector<std::vector<float>>* weights;
	std::string system;

	// Simulation parameters
	uint32_t seed;

	// Timing
	double time_taken;
	uint32_t time_finished;
	

};

class SimulationParametersFull
{
	// Individual parameters
	SimulationParametersIndividual * parameters[10];
	uint8_t number_of_layers = 0;
	
public:
	//SimulationParametersFull();
	~SimulationParametersFull();
	void clearLayers();
	void addLayer(SimulationParametersIndividual * layer);
	void serialize();
	std::string serialization = "";
	uint32_t deposited = 0;
	float turns = 0;
};
