#pragma once
#include<cstdint>
#include <random>
#include <math.h>
#include <iostream>



class Surface3D
{
	uint16_t L;
	uint16_t H;
	uint16_t W;

	int8_t* grid = nullptr;
	int8_t* adjacency = nullptr;
	double* energy_grid = nullptr;
	std::vector<int8_t>* species = nullptr;
	std::vector<std::vector<float>>* weights = nullptr;

	uint16_t* neighbors = nullptr;
	float* energies = nullptr;

	uint32_t seed;
	std::mt19937* gen;
	//std::uniform_int_distribution<>* dist;



public:
	Surface3D(uint16_t length, uint16_t height, uint16_t width, int8_t* grid, std::vector<int8_t>* species, std::vector<std::vector<float>>* weights, uint32_t seed);

	~Surface3D();

	uint8_t add(uint16_t* center, int8_t sp);
	uint8_t remove(uint16_t* center, int8_t sp);
	uint8_t getNearestNeighbors(uint16_t* center);
	uint16_t* getAdjacentVacancy(uint16_t* center, int8_t sp);
	double calculateLocalEnergy(uint16_t* center, int8_t sp);


};

