#pragma once
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cstdint>
#include <random>
#include <chrono>

#include "SimulationParameters.h"
#include "SlantedCorridors.h"
#include "CubicSpacePartition.h"


float distance3D(uint16_t* p1, uint16_t* p2);
float disance3D_PBC(uint16_t* p1, uint16_t* p2, uint16_t L);
float distanceSquare3D(uint16_t* p1, uint16_t* p2);
float distanceSquare3D_PBC(uint16_t* p1, uint16_t* p2, uint16_t L);

int obliqueDepositionContinuous(float theta, float L, float H, uint32_t reps, uint8_t bin_size, uint32_t seed, float diffusion_length, float length_scale, std::vector<int8_t>* species, std::vector<float>* radii, std::vector<std::vector<float>>* weights, std::vector<std::vector<float>> inputGrid, SimulationParametersFull* params, std::string& system);