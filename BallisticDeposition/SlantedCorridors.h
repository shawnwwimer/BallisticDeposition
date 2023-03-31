#pragma once
#include <cstdint>
#include <array>
#include <set>
#include <vector>
#include <iterator>
#include <iostream>

//#include "math_utils.h"

struct particle_priority
{
	uint32_t idx;
	float priority;
	particle_priority(uint32_t idx, float priority) : idx{ idx }, priority{ priority }{}
	//bool operator <(particle_priority a) { return (priority < a.priority); }
	//bool operator >(particle_priority a) { return (priority > a.priority); }
	//bool operator ==(particle_priority a) { return (priority == a.priority); }
	friend bool operator <(particle_priority a, particle_priority b) { return a.priority > b.priority; }
	friend bool operator >(particle_priority a, particle_priority b) { return a.priority < b.priority; }
	friend bool operator ==(particle_priority a, particle_priority b) { return (a.priority == b.priority) && (a.idx == b.idx); }
};

struct collision_description
{
	std::array<float, 3> position;
	int64_t idx;
	collision_description(float x, float y, float z, int64_t idx) : idx{ idx }
	{
		position[0] = x;
		position[1] = y;
		position[2] = z;
	}
};

float modulof(float val, float mod);


class SlantedCorridors
{
private:
	float L;
	float Lfloat;
	float H;
	float theta;
	float tan_theta;
	float sin_theta;
	float cos_theta;

	float max_sub;
	float max_z;

	float pi = 3.141592653;

	uint8_t bin_size;
	uint8_t bins_on_side;
	uint8_t bins_num;
	std::vector<particle_priority*> particles;
	//bool comp(particle_priority a, particle_priority b) { return a.priority > b.priority; }; // need a comparison for the sets
	std::vector<std::set<particle_priority>> bins;

	// For use by functions
	std::array<std::set<particle_priority>*, 8> bins_found = { 0 };
	uint8_t bins_found_num = 0;
	collision_description collision = collision_description(0, 0, 0, -1);
	//float dropped_point[3] = { 0.f };

	// Private methods
	uint32_t find_bin_idx(std::array<float, 3> position)
	{
		float adjusted_x = position[2] * tan_theta + position[0];
		uint8_t spacex = (uint8_t)floor(modulof(adjusted_x, L) / bin_size);
		uint8_t spacey = (uint8_t)floor(modulof(position[1], L) / bin_size);
		return spacex + spacey * bins_on_side;
	}
	std::set<particle_priority>* find_bin(std::array<float, 3> position)
	{
		return &bins[find_bin_idx(position)];
	}
	std::set<particle_priority>* find_bins(std::array<float, 3> position, float radius);
	
	//std::vector<uint32_t> get_nearby_particles(float* position, float radius);
	float calc_priority(std::array<float, 3> position)
	{
		float x1 = position[0] + position[2] * tan_theta;
		float sub = (L - x1) * sin_theta + max_sub * floor(position[2]/max_z);
		float above = position[2] / cos_theta;
		return sub + above;
		//return position[2];
	}


public:
	SlantedCorridors(float L, float H, float theta, uint8_t bin_size) : L{ L }, H{ H }, theta{ theta } , bin_size{ bin_size }
	{
		Lfloat = (float)L;
		tan_theta = tan(theta * pi / 180.0);
		sin_theta = sin(theta * pi / 180.0);
		cos_theta = cos(theta * pi / 180.0);
		max_sub = L * sin_theta;
		max_z = L * tan((90 - theta) * pi / 180.0);

		bins_on_side = (uint8_t)(L / bin_size);
		bins_num = bins_on_side * bins_on_side;
		/*
		bins = (std::set<particle_priority>**)malloc(sizeof(std::set<particle_priority>**));
		if (bins == nullptr) {
			throw;
		}
		*bins = (std::set<particle_priority>*)malloc(sizeof(std::set<particle_priority>*) * bins_num);
		
		for (uint8_t i = 0; i < bins_num; i++) {
			(bins)[i] = new std::set<particle_priority>;
		}
		*/
		for (uint8_t i = 0; i < bins_num; i++) {
			std::set<particle_priority>* b = new std::set<particle_priority>;
			bins.push_back(*b);
		}
		//bins.push_back(new std::set<particle_priority>);

		std::cout << "After iteration" << std::endl;

	}
	~SlantedCorridors()
	{
		//for (int i = 0; i < bins_num; i++) {
		//	delete(bins[i]);
		//}
		//free(*bins);

		/*for (auto p : particles) {
			delete(p);
		}*/
	}

	std::set<particle_priority>* add_to_bins(std::array<float, 3> position, float radius, uint32_t idx)
	{
		float priority = calc_priority(position) + idx/100000.0;
		find_bins(position, radius);
		particle_priority* x = new particle_priority(idx, priority);
		particles.push_back(x);

		for (int i = 0; i < bins_found_num; i++) {
			bins_found[i]->insert(*x);
		}

		return nullptr;
	}

	float get_particle_priority(uint32_t idx) {
		return particles[idx]->priority;
	}

	collision_description* drop_particle(std::array<float, 3> position, float radius, std::vector<std::vector<float>> atoms);




};

