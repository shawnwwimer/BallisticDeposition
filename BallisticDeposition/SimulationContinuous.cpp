#include "SimulationContinuous.h"
#include "cnpy.h"

int obliqueDepositionContinuous(float theta, float L, float H, uint32_t reps, uint8_t bin_size, uint32_t seed, float diffusion_length, float length_scale, std::vector<int8_t>* species, std::vector<float>* radii, std::vector<std::vector<float>>* weights, std::vector<std::vector<float>> inputGrid, SimulationParametersFull* params, std::string& system)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<std::vector<float>> atoms;

	SlantedCorridors corridors = SlantedCorridors(L, H, theta, bin_size);

	// Set up potentials for diffusion
	float* potentials = (float*)malloc(sizeof(float) * L * length_scale * L * length_scale * H * length_scale * 4);
	for (int i = 0; i < L * L * H * length_scale * length_scale * length_scale * 4; i++) {
		potentials[i] = 0;
	}

	float A = 1;
	float R = 1;
	float Vm = length_scale * (2 * (*radii)[0]);
	float zero_pot = Vm / pow(2, 1.0 / 6.0);

	uint32_t diameter = 2 * length_scale;
	if (diameter % 2 == 0) {
		diameter += 1;
	}

	float* region = (float* )malloc(sizeof(float) * diameter * diameter * diameter * 4);
	for (int k = 0; k < diameter; k++) {
		int kk = k - (diameter - 1) / 2;
		for (int j = 0; j < diameter; j++) {
			int jj = j - (diameter - 1) / 2;
			for (int i = 0; i < diameter; i++) {
				int ii = i - (diameter - 1) / 2;

				if (ii == 0 && jj == 0 && kk == 0) {
					region[i, j, k, 0] = 0;
					region[i, j, k, 1] = 0;
					region[i, j, k, 2] = 0;
					region[i, j, k, 3] = 1e6;
				}
				else {
					float r = sqrt(ii * ii + jj * jj + kk * kk);
					float mag = A * pow((zero_pot / r), 12) - R * pow((zero_pot / r), 6);
					region[k * diameter * diameter + j * diameter + i + 3] = mag;
					if (r < Vm) {
						region[k * diameter * diameter + j * diameter + i] = ((float)ii)/r;
						region[k * diameter * diameter + j * diameter + i + 1] = ((float)jj) / r;
						region[k * diameter * diameter + j * diameter + i + 2] = ((float)kk) / r;
					}
					else {
						region[k * diameter * diameter + j * diameter + i] = -((float)ii) / r;
						region[k * diameter * diameter + j * diameter + i + 1] = -((float)jj) / r;
						region[k * diameter * diameter + j * diameter + i + 2] = -((float)kk) / r;
					}
				}
			}
		}
	}
	int steps = round(length_scale * diffusion_length);
	int rads = round(length_scale * 2 * (*radii)[0]);
	float scaled_position[3] = { 0 };

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

		// TODO: DIFFUSION
		auto startd = std::chrono::high_resolution_clock::now();
		if (diffusion_length > 0) {
			scaled_position[0] = modulof(collision->position[0] * length_scale, L);
			scaled_position[1] = modulof(collision->position[1] * length_scale, L);
			scaled_position[2] = collision->position[2] * length_scale;

			int remaining_distance = steps;
			while (remaining_distance > 0) {
				int it = (remaining_distance > rads) ? remaining_distance : it;
				remaining_distance -= it;

				int minx = scaled_position[0] - it;
				int miny = scaled_position[1] - it;
				int minz = scaled_position[2] - it;
			}

		}
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
		cnpy::npy_save("structures/cts/grid.npy", outGrid, { reps, 6 });
		free(outGrid);
	}

	// Print the timing
	std::cout << reps << " reps completed in " << time << " seconds; \n" << reps / time << " reps per second." << std::endl;
	std::cout << "Line finding and traversal were " << timel << " seconds of that time." << std::endl;
	std::cout << "Diffusion was " << timed << " seconds of that time." << std::endl;

	free(potentials);
	return reps;
}