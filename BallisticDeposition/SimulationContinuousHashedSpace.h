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
#include "SpaceHashMap.h"
#include "OverlappingCubicSpacePartition.h"
#include "Matrix3DLateralPBC.h"
#include "SimulationContinuous.h"
#include "FileWriting.h"

int obliqueDepositionContinuousHashed(float theta, float L, float H, uint32_t reps, float bin_size, uint32_t seed, float diffusion_length, float length_scale, std::vector<int8_t>* species, std::vector<float>* radii, std::vector<float>* spread, std::vector<std::vector<float>>* weights, std::vector<std::vector<float>> inputGrid, ContinuousSimulationParametersFull* params, std::string& system, DiffusionMethod diffusion_method, FilesToSave* save);