#include "Surface3D.h"
#include "math_utils.h"

Surface3D::Surface3D(uint16_t length, uint16_t height, uint16_t width, int8_t* grid, std::vector<int8_t>* species, std::vector<std::vector<float>>* weights, uint32_t seed)
{
	// Initialize parameters
	L = length;
	H = height;
	W = width;
	this->grid = grid;
	this->species = species;
	this->weights = weights;

	// Start random number generator
	if (seed == 0) {
		std::random_device rd;
		this->seed = rd();
	}
	else {
		this->seed = seed;
	}
	this->gen = new std::mt19937(this->seed);
	//this->dist = new std::uniform_int_distribution<>(0, L - 1);

	// Allocate memory
	neighbors = (uint16_t*)malloc(sizeof(uint16_t) * 27 * 3); // 26 neighbors or 27 points maximum in 3x3x3 cube, three coordinates
	adjacency = (int8_t*)malloc(sizeof(int8_t) * L * W * H); // Size of volume
	energy_grid = (double*)malloc(sizeof(double) * L * W * H * species->size());
	energies = (float*)malloc(sizeof(float) * 27); // Number of points in 3x3x3

	// Initialize, first to zero and then update based on the grid
	uint16_t point[3] = { 0, 0, 0 };
	for (uint32_t i = 0; i < L*W*H; i++) {
		adjacency[i] = 0;

		wrap_index(i, point, L, H);
		for (int32_t s = 0; s < species->size(); s++) {
			energy_grid[i * species->size() + s] = 0;
		}
		for (int32_t s = 0; s < species->size(); s++) {
			if (point[2] > 0 && point[2] < H - 1)
			{
				energy_grid[i * species->size() + s] += (*weights)[(*species)[s]][0] * (6 + 12 / 1.414 + 8 / 1.732);
			}
			else
			{
				energy_grid[i * species->size() + s] += (*weights)[(*species)[s]][0] * (5 + 8 / 1.414 + 4 / 1.732);
			}
		}
}
	uint32_t jkjk = 0;
	for (uint32_t i = 0; i < L * W * H; i++) {
		if (grid[i] > 0) {
			wrap_index(i, point, L, H);
			add_directly(point, grid[i] - 1);
			jkjk++;
		}
	}
	/*
	for (uint32_t i = 0; i < L * W * H; i++) {
		wrap_index(i, point, L, H);
		for (uint8_t j = 0; j < species->size(); j++) {
			energy_grid[i * species->size() + j] = calculateLocalEnergy(point, j);
		}
	}
	*/
}

Surface3D::~Surface3D()
{
	free(neighbors);
	free(adjacency);
	free(energies);
	free(energy_grid);
	delete gen;
	//delete dist;
}

uint8_t Surface3D::add_directly(uint16_t* center, int8_t sp)
{
	uint8_t n = 0;
	uint16_t point[3] = { 0, 0, 0 };
	for (int32_t i = -1; i < 2; i++) {
		uint16_t ii = modulo(center[0] + i, L);
		point[0] = ii;
		for (int32_t j = -1; j < 2; j++) {
			uint16_t jj = modulo(center[1] + j, L);
			point[1] = jj;
			for (int32_t k = -1; k < 2; k++) {
				int32_t kk = center[2] + k;

				// Don't bother with itself or out of z-range
				if ((i == 0 && j == 0 && k == 0) || kk < 0 || kk > H - 1) {
					continue;
				}
				point[2] = (uint16_t)kk;
				uint32_t idx = flat_index(point, L, L);

				// add to coordination number
				adjacency[idx] += 1;
				n += 1;

				// update energies
				for (int32_t s = 0; s < species->size(); s++) {
					if ((i == 0 && j == 0) || (i == 0 && k == 0) || (j == 0 && k == 0)) {
						// Shares side
						energy_grid[idx * species->size() + s] += (*weights)[(*species)[s]][sp];
						energy_grid[idx * species->size() + s] -= (*weights)[(*species)[s]][0];
					}
					else if (i == 0 || j == 0 || k == 0) {
						// Shares edge
						energy_grid[idx * species->size() + s] += (*weights)[(*species)[s]][sp] / 1.414;
						energy_grid[idx * species->size() + s] -= (*weights)[(*species)[s]][0] / 1.414;
					}
					else {
						//Shares corner
						energy_grid[idx * species->size() + s] += (*weights)[(*species)[s]][sp] / 1.732;
						energy_grid[idx * species->size() + s] -= (*weights)[(*species)[s]][0] / 1.732;
					}
				}
			}
		}
	}

	return n;
}

uint8_t Surface3D::add(uint16_t* center, int8_t sp)
{
	uint8_t n = 0;
	uint16_t point[3] = { 0, 0, 0 };
	for (int32_t i = -1; i < 2; i++) {
		uint16_t ii = modulo(center[0] + i, L);
		point[0] = ii;
		for (int32_t j = -1; j < 2; j++) {
			uint16_t jj = modulo(center[1] + j, L);
			point[1] = jj;
			for (int32_t k = -1; k < 2; k++) {
				int32_t kk = center[2] + k;

				// Don't bother with itself or out of z-range
				if ((i == 0 && j == 0 && k == 0) || kk < 0 || kk > H - 1) {
					continue;
				}
				point[2] = (uint16_t)kk;
				uint32_t idx = flat_index(point, L, L);

				// add to coordination number
				adjacency[idx] += 1;
				n += 1;

				// update energies
				for (int32_t s = 0; s < species->size(); s++) {
					if ((i == 0 && j == 0) || (i == 0 && k == 0) || (j == 0 && k == 0)) {
						// Shares side
						energy_grid[idx * species->size() + s] += (*weights)[(*species)[s]][(*species)[sp]];
						energy_grid[idx * species->size() + s] -= (*weights)[(*species)[s]][0];
					}
					else if (i == 0 || j == 0 || k == 0) {
						// Shares edge
						energy_grid[idx * species->size() + s] += (*weights)[(*species)[s]][(*species)[sp]] / 1.414;
						energy_grid[idx * species->size() + s] -= (*weights)[(*species)[s]][0] / 1.414;
					}
					else {
						//Shares corner
						energy_grid[idx * species->size() + s] += (*weights)[(*species)[s]][(*species)[sp]] / 1.732;
						energy_grid[idx * species->size() + s] -= (*weights)[(*species)[s]][0] / 1.732;
					}
				}
			}
		}
	}

	return n;
}

uint8_t Surface3D::remove(uint16_t* center, int8_t sp)
{
	uint8_t n = 0;
	uint16_t point[3] = { 0, 0, 0 };
	for (int32_t i = -1; i < 2; i++) {
		uint16_t ii = modulo(center[0] + i, L);
		point[0] = ii;
		for (int32_t j = -1; j < 2; j++) {
			uint16_t jj = modulo(center[1] + j, L);
			point[1] = jj;
			for (int32_t k = -1; k < 2; k++) {
				int32_t kk = center[2] + k;

				// Don't bother with itself or out of z-range
				if ((i == 0 && j == 0 && k == 0) || kk < 0 || kk > H - 1) {
					continue;
				}
				point[2] = (uint16_t)kk;
				uint32_t idx = flat_index(point, L, L);

				// add to coordination number
				adjacency[idx] -= 1;
				n += 1;

				// update energies
				for (int32_t s = 0; s < species->size(); s++) {
					energy_grid[idx * species->size() + s] -= (*weights)[s][sp-1];
				}
			}
		}
	}

	return n;
}

uint8_t Surface3D::getNearestNeighbors(uint16_t* center)
{
	uint32_t n = 0;
	uint16_t point[3] = { 0, 0, 0 };
	for (int32_t i = -1; i < 2; i++) {
		uint16_t ii = modulo(center[0] + i, L);
		point[0] = ii;
		for (int32_t j = -1; j < 2; j++) {
			uint16_t jj = modulo(center[1] + j, L);
			point[1] = jj;
			for (int32_t k = -1; k < 2; k++) {
				int32_t kk = center[2] + k;

				// Don't bother with itself or out of z-range
				if ((i == 0 && j == 0 && k == 0) || kk < 0 || kk > H - 1) {
					continue;
				}
				point[2] = (uint16_t)kk;

				// add to neighbor list
				if (grid[flat_index(point, L, L)] > 0) {
					neighbors[n * 3] = ii;
					neighbors[n * 3 + 1] = jj;
					neighbors[n * 3 + 2] = kk;
					n += 1;
				}
			}
		}
	}

	return n;
}

uint16_t* Surface3D::getAdjacentVacancy(uint16_t* center, int8_t sp)
{
	// Get every vacancy and calculate energies
	float max = 0;
	float min = 0;
	float energy = 0;
	uint32_t idx = 0;
	uint32_t n = 0;
	uint16_t point[3] = { 0, 0, 0 };
	for (int32_t i = -1; i < 2; i++) {
		uint16_t ii = modulo(center[0] + i, L);
		point[0] = ii;
		for (int32_t j = -1; j < 2; j++) {
			uint16_t jj = modulo(center[1] + j, L);
			point[1] = jj;
			for (int32_t k = -1; k < 2; k++) {
				int32_t kk = center[2] + k;

				// Don't bother with out of z-range
				if (kk < 0 || kk > H - 1) {
					continue;
				}
				point[2] = (uint16_t)kk;

				// add to neighbor list
				idx = flat_index(point, L, L);
				if (grid[idx] == 0 && (adjacency[idx] > 0 || kk == 0)) {
					// Add vacancy to list
					neighbors[n * 3] = ii;
					neighbors[n * 3 + 1] = jj;
					neighbors[n * 3 + 2] = kk;

					// If using weights, add to energy list
					if (weights != nullptr) {
						energy = energy_grid[idx * species->size() + sp];//calculateLocalEnergy(center, sp);//
						if (energy > max) {
							max = energy;
						}
						else if (energy < min) {
							min = energy;
						}
						energies[n] = energy;
					}

					n += 1;
				}
			}
		}
	}

	if (n == 0) {
		std::cout << "Particle landed at (" << center[0] << ", " << center[1] << ", " << center[2] << ") with " << n << " adjacent vacancies." << std::endl;
	}

	
	if (weights == nullptr) {
		// Choose from any index
		std::uniform_int_distribution<> dist(0, n - 1);
		return &neighbors[dist(*gen)*3];
	}
	else {
		float total_energy = 0;
		// Reprocess weights if negative values exist
		if (max < 0) {
			// Shift energies to range (0, ..., max - min)
			for (uint8_t i = 0; i < n; i++) {
				energies[i] -= min;
				total_energy += energies[i];
			}
		}
		else if (min < 0) {
			// Set negative energies to 0 weight
			for (uint8_t i = 0; i < n; i++) {
				if (energies[i] < 0) {
					energies[i] = 0;
				}
				total_energy += energies[i];
			}
		} else {
			for (uint8_t i = 0; i < n; i++) {
				total_energy += energies[i];
			}
		}

		if (total_energy == 0) {
			// Choose from any index
			std::uniform_int_distribution<> dist(0, n - 1);
			return &neighbors[dist(*gen) * 3];
		}

		// Get a random-ish vacancy
		std::uniform_real_distribution<> dist(0, total_energy);
		float val = dist(*gen);
		float run = 0;
		for (uint32_t i = 0; i < n; i++) {
			run += energies[i];
			if (run >= val) {
				return &neighbors[i * 3];
			}
		}
		std::cout << "Overrun: val " << val << "run " << run << " n " << n << " &neighbors[(n-1) * 3] - neighbors " << &neighbors[(n-1) * 3] - neighbors << std::endl;
		return &neighbors[(n-1) * 3];
	}
}

double Surface3D::calculateLocalEnergy(uint16_t* center, int8_t sp) {
	double energy = 0;
	uint16_t point[3] = { 0, 0, 0 };
	for (int32_t i = -1; i < 2; i++) {
		uint16_t ii = modulo(center[0] + i, L);
		point[0] = ii;
		for (int32_t j = -1; j < 2; j++) {
			uint16_t jj = modulo(center[1] + j, L);
			point[1] = jj;
			for (int32_t k = -1; k < 2; k++) {
				int32_t kk = center[2] + k;

				// Don't bother with itself or out of z-range
				if ((i == 0 && j == 0 && k == 0) || kk < 0 || kk > H - 1) {
					continue;
				}
				point[2] = (uint16_t)kk;

				// get neighbor so we can skip lookup in weights
				int8_t neighbor = grid[flat_index(point, L, L)];
				if (neighbor > 0) {
					if ((i == 0 && j == 0) || (i == 0 && k == 0) || (j == 0 && k == 0)) {
						// Shares side
						energy += (*weights)[(*species)[sp]][neighbor];
					}
					else if (i == 0 || j == 0 || k == 0) {
						// Shares edge
						energy += (*weights)[(*species)[sp]][neighbor] / 1.414;
					}
					else {
						//Shares corner
						energy += (*weights)[(*species)[sp]][neighbor] / 1.732;
					}
				}
			}
		}
	}
	return energy;
}