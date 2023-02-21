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
#include "Surface3D.h"



uint16_t Bresenham3D(uint16_t* src, uint16_t* dest, uint16_t* points);

uint16_t* traversePathCalculated(uint16_t* points, uint16_t points_len, uint16_t* dest, uint16_t* collision, int8_t* grid, uint16_t L, uint16_t H);

uint16_t traversePathRealTime(int32_t* src, int32_t* dest, uint16_t* collision, int8_t* grid, uint16_t L, uint16_t H);

uint32_t denseToSparse(int8_t* grid, int16_t** sparse, uint16_t L, uint16_t H);
uint32_t denseToSparse(uint32_t* grid, uint32_t** sparse, uint16_t L, uint16_t H);
uint16_t sparseToDense(int16_t* sparse, uint32_t num_points, int8_t* grid, uint16_t L, uint16_t H);

int obliqueDeposition(float theta, uint16_t L, uint16_t H, uint32_t reps, float phi, float turns, uint32_t seed, uint16_t diffusion_steps, std::vector<int8_t> * species, std::vector<float> * spread, std::vector<std::vector<float>> * weights, int16_t* inputGrid, uint32_t inputGridPoints, int16_t** outGrid, uint32_t stepper_resolution, SimulationParametersFull* params, std::string& system);



