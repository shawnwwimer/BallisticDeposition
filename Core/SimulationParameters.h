#pragma once
#include <cstdint>
#include <vector>
#include <string>

class SimulationParametersIndividual
{
public:
	// Physical parameters
	uint16_t length = 0;
	uint16_t width = 0;
	uint16_t height = 0;
	float theta = 0;
	float phi = 0;
	float turns = 0;
	std::vector<float> * spread = nullptr;
	int phi_num = 0;
	float phi_deg = 0;
	float theta_end = 0;
	uint32_t stepper_resolution = 0;
	
	// Material parameters
	std::vector<int8_t> * species = nullptr;
	uint16_t diffusion_steps = 0;
	uint32_t deposited = 0;
	uint32_t repetitions = 0;
	std::vector<std::vector<float>>* weights = nullptr;
	std::string system = "";

	// Simulation parameters
	uint32_t seed = 0;
	uint8_t acceleration = 0;
	uint8_t collision_method = 0;

	// Timing
	double time_taken = 0;
	uint32_t time_finished = 0;
	

};

class SimulationParametersFull
{
	// Individual parameters
	SimulationParametersIndividual * parameters[32];
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

class ContinuousSimulationParameters {
public:
	// Physical parameters
	uint16_t length;
	uint16_t width;
	uint16_t height;
	float theta;

	// Material parameters
	std::vector<int8_t>* species;
	std::vector<float>* radii;
	float diffusion_length;
	uint32_t repetitions;
	std::string system;

	// Simulation parameters
	uint32_t seed;
	uint8_t bin_size;
	float length_scale;
	float cube_size;
	int diffusion_method;

	// Timing
	double time_taken;
	uint32_t time_finished;
};

class ContinuousSimulationParametersFull
{
	// Individual parameters
	ContinuousSimulationParameters* parameters[10];
	uint8_t number_of_layers = 0;

public:
	//SimulationParametersFull();
	~ContinuousSimulationParametersFull();
	void clearLayers();
	void addLayer(ContinuousSimulationParameters* layer);
	void serialize();
	std::string serialization = "";
	uint32_t deposited = 0;
};