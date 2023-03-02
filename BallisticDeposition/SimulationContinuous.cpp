#include "SimulationContinuous.h"
#include "cnpy.h"

int obliqueDepositionContinuous(float theta, float L, float H, uint32_t reps, uint8_t bin_size, uint32_t seed, uint16_t diffusion_length, float length_scale, std::vector<int8_t>* species, std::vector<float>* radii, std::vector<std::vector<float>>* weights, std::vector<std::vector<float>> inputGrid, SimulationParametersFull* params, std::string& system)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<std::vector<float>> atoms;

	SlantedCorridors corridors = SlantedCorridors(L, H, theta, bin_size);

	// TODO: set up potentials for diffusion


	// Create random number generator
	if (seed == 0) {
		std::random_device rd;
		seed = rd();
	}
	std::mt19937 gen(seed);
	std::uniform_real_distribution<float> dist(0, L);

	std::mt19937 gens(seed);
	std::uniform_int_distribution<> sist(0, species->size());

	// Print how long it took
	std::chrono::duration<double, std::milli> timep = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Set-up took " << timep.count() / 1000 << " s." << std::endl;

	// Begin
	start = std::chrono::high_resolution_clock::now();
	int update = 128; // number of printed updates
	std::array<float, 3> dest;
	std::vector<float> new_atom;
	float current_fiber_id = 1;
	float fiber = 0;
	double timel = 0;
	double timed = 0;
	for (int n = 0; n < reps; n++) {
		auto startl = std::chrono::high_resolution_clock::now();

		// Generate 
		dest = { dist(gen), dist(gen), 0 };//{ ((float)n)  / reps + L / 2 * (n % 2), (float)n / reps + L / 2 * (n % 2), 0 };//

		// Choose species
		int sp = 0;

		// Drop particle
		collision_description* collision = corridors.drop_particle(dest, (*radii)[sp], atoms);

		std::chrono::duration<double, std::milli> dural = std::chrono::high_resolution_clock::now() - startl;
		timel += dural.count() / 1000;

		// TODO: diffusion
		auto startd = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> durad = std::chrono::high_resolution_clock::now() - startd;
		timed += durad.count() / 1000;

		// New fiber?
		if (collision->idx == -1) {
			fiber = current_fiber_id;
			current_fiber_id += 1;
		}
		else {
			fiber = atoms[collision->idx][5];
		}

		// Add the new atom
		new_atom = { collision->position[0], collision->position[1], collision->position[2], (float)(*species)[sp], (*radii)[sp], fiber };
		atoms.push_back(new_atom);
		corridors.add_to_bins(collision->position, (*radii)[sp], n);

		// Periodic updates
		if (n % (reps / update) == (reps / update) - 1) {
			std::chrono::duration<double, std::milli> durup = std::chrono::high_resolution_clock::now() - start;
			std::cout << n + 1 << "/" << reps << " complete; " << durup.count() / 1000 << " seconds taken; " << (n + 1) / durup.count() * 1000 << " reps per second." << std::endl;
		}
	}

	// Finish up timing
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> dural = end - start;
	double time = dural.count() / 1000;

	std::chrono::system_clock::time_point epoch = std::chrono::system_clock::now();
	uint32_t epoch_time = epoch.time_since_epoch().count() * std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;

	float* outGrid = (float*)malloc(sizeof(float) * reps * 6);
	if (outGrid != NULL) {
		for (int j = 0; j < atoms.size(); j++) {
			for (int i = 0; i < 6; i++) {
				outGrid[j * 6 + i] = atoms[j][i];
			}
		}
		cnpy::npy_save("grid.npy", outGrid, { reps, 6 });
		free(outGrid);
	}

	// Print the timing
	std::cout << reps << " reps completed in " << time << " seconds; \n" << reps / time << " reps per second." << std::endl;
	std::cout << "Line finding and traversal were " << timel << " seconds of that time." << std::endl;
	std::cout << "Diffusion was " << timed << " seconds of that time." << std::endl;

	return reps;
}