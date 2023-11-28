#pragma once
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cstdint>
#include <random>
#include <chrono>

#include "nlopt.hpp"

#include "SimulationParameters.h"
#include "SlantedCorridors.h"
#include "OverlappingCubicSpacePartition.h"
#include "Matrix3DLateralPBC.h"


float distance3D(uint16_t* p1, uint16_t* p2);
float disance3D_PBC(uint16_t* p1, uint16_t* p2, uint16_t L);
float distanceSquare3D(uint16_t* p1, uint16_t* p2);
float distanceSquare3D_PBC(uint16_t* p1, uint16_t* p2, uint16_t L);

enum class DiffusionMethod {
	PotentialHoppingLUT = 0,
	ForcePushingLUT = 1,
	HopAndSettleLUT = 2,
	NumericalMinimization = 3,
};

// The grid and params are always saved
// Defaults:
//   destinations and collisions: true
//   atomic_potential, volume_potential, and priority: false
struct FilesToSave {
	bool destinations = true;
	bool collisions = true;
	bool atomic_potential = false;
	bool volume_potential = false;
	bool priority = false;
};

int obliqueDepositionContinuous(float theta, float L, float H, uint32_t reps, float bin_size, uint32_t seed, float diffusion_length, float length_scale, std::vector<int8_t>* species, std::vector<float>* radii, std::vector<float>* spread, std::vector<std::vector<float>>* weights, std::vector<std::vector<float>> inputGrid, ContinuousSimulationParametersFull* params, std::string& system, DiffusionMethod diffusion_method, FilesToSave * save);