#include "SimulationContinuous.h"
#include "cnpy.h"
#include "zip.h"

#define COMPRESSION_LEVEL (Z_DEFAULT_COMPRESSION)
#define VERBOSE (false)

int writeFileToZipCTS(const char* zipname, const char* filename)
{
	int error = ZIP_ERRNO;

	zipFile zf = zipOpen64(zipname, APPEND_STATUS_ADDINZIP);

	// Check for a valid zipfile
	if (zf == NULL) {
		return ZIP_BADZIPFILE;
	}

	// Attempt to open the file
	std::fstream file(filename, std::ios::binary | std::ios::in);
	if (file.is_open()) {
		// Get size of file
		file.seekg(0, std::ios::end);
		size_t size = file.tellg();
		file.seekg(0, std::ios::beg);

		// Read in file 
		// TODO: might want to chunk it since some files may be large
		std::vector<char> buffer(size);
		if (size == 0 || file.read(&buffer[0], size)) {
			// Initialize the parameters for the local header
			zip_fileinfo zi = { 0 };

			// Open it inside the zip for writing
			if (Z_OK == zipOpenNewFileInZip64(zf, filename, &zi, NULL, 0, NULL, 0, NULL, Z_DEFLATED, COMPRESSION_LEVEL, 1)) {
				// Write in zip
				if (zipWriteInFileInZip(zf, size == 0 ? "" : &buffer[0], size)) {
					error = ZIP_ERRNO;
				}

				// Close it inside zip
				if (zipCloseFileInZip(zf)) {
					error = ZIP_OK;
				}
				// Close file
				file.close();
			}
		}
	}

	// Try to close zip
	if (zipClose(zf, NULL)) {
		return ZIP_ERRNO;
	}
	return error;
}

int obliqueDepositionContinuous(float theta, float L, float H, uint32_t reps, uint8_t bin_size, uint32_t seed, float diffusion_length, float length_scale, std::vector<int8_t>* species, std::vector<float>* radii, std::vector<std::vector<float>>* weights, std::vector<std::vector<float>> inputGrid, ContinuousSimulationParametersFull* params, std::string& system)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<std::vector<float>> atoms;

	SlantedCorridors corridors = SlantedCorridors(L, H, theta, bin_size);

	CubicSpacePartition cubes = CubicSpacePartition(L, H, 2.0);
	float s = 2*(*radii)[0] / pow(2, 1.0 / 6.0);
	float s6 = pow(s, 6);
	float s12 = pow(s6, 2);

	Matrix3DLateralPBC potentials = Matrix3DLateralPBC(L, L, H, length_scale);

	float A = 1;
	float R = 1;
	float Vm = length_scale * (2 * (*radii)[0]);
	float zero_pot = Vm / pow(2, 1.0 / 6.0);

	Matrix3DLateralPBC region = Matrix3DLateralPBC(2.f, 2.f, 2.f, length_scale);

	int diameter = region.get_Hs();

	//float* region = (float* )malloc(sizeof(float) * diameter * diameter * diameter * 4);
	for (int k = 0; k < diameter; k++) {
		int kk = k - (diameter - 1) / 2;
		for (int j = 0; j < diameter; j++) {
			int jj = j - (diameter - 1) / 2;
			for (int i = 0; i < diameter; i++) {
				int ii = i - (diameter - 1) / 2;

				if (ii == 0 && jj == 0 && kk == 0) {
					/*region[i, j, k, 0] = 0;
					region[i, j, k, 1] = 0;
					region[i, j, k, 2] = 0;
					region[i, j, k, 3] = 1e6;*/
					region(i, j, k) = 1e6;
				}
				else {
					float r = sqrt(ii * ii + jj * jj + kk * kk);
					float mag = A * pow((zero_pot / r), 12) - R * pow((zero_pot / r), 6);
					//region[k * diameter * diameter + j * diameter + i + 3] = mag;
					if (r < Vm) {
						/*region[k * diameter * diameter + j * diameter + i] = ((float)ii)/r;
						region[k * diameter * diameter + j * diameter + i + 1] = ((float)jj) / r;
						region[k * diameter * diameter + j * diameter + i + 2] = ((float)kk) / r;*/
						/*(*region(i, j, k))[0] = mag * ((float)ii) / r;
						(*region(i, j, k))[1] = mag * ((float)jj) / r;
						(*region(i, j, k))[2] = mag * ((float)kk) / r;*/
						region(i, j, k) = mag;
					}
					else {
						/*region[k * diameter * diameter + j * diameter + i] = -((float)ii) / r;
						region[k * diameter * diameter + j * diameter + i + 1] = -((float)jj) / r;
						region[k * diameter * diameter + j * diameter + i + 2] = -((float)kk) / r;*/
						/*(*region(i, j, k))[0] = mag * ((float)ii) / r;
						(*region(i, j, k))[1] = mag * ((float)jj) / r;
						(*region(i, j, k))[2] = mag * ((float)kk) / r;*/
						region(i, j, k) = mag;
					}
				}
			}
		}
	}
	int steps = round(length_scale * diffusion_length);
	int rads = round(length_scale * 2 * (*radii)[0]);
	float scaled_position[3] = { 0 };
	std::array<float, 3> direction = { 0, 0, 0 };
	std::array<float, 3> current_minimum = { 0, 0, 0 };
	std::array<int, 3> center = { 0, 0, 0 };

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
	int update = (reps > 128) ? 128 : 1; // number of printed updates
	std::array<float, 3> dest;
	std::vector<float> new_atom;
	float current_fiber_id = 1;
	float fiber = 0;
	double timel = 0;
	double timed = 0;
	for (int n = 0; n < reps; n++) {
		auto startl = std::chrono::high_resolution_clock::now();

		// Generate 
		//dest = { dist(gen), dist(gen), 0 };//{ ((float)n)  / reps + L / 2 * (n % 2), (float)n / reps + L / 2 * (n % 2), 0 };//
		dest = { 16.0f+0*powf(n/12.8f, 2), 0, 0 };

		// Choose species
		int sp = 0;

		// Drop particle
		collision_description* collision = corridors.drop_particle(&dest, (*radii)[sp], &atoms);
		
		std::chrono::duration<double, std::milli> dural = std::chrono::high_resolution_clock::now() - startl;
		timel += dural.count() / 1000;

		if (VERBOSE) {
			std::cout << "Particles: " << std::endl;
			for (int i = 0; i < atoms.size(); i++) {
				auto p = atoms[i];
				std::cout << "\t(" << p[0] << ", " << p[2] << ")," << std::endl;// " with priority " << corridors.get_particle_priority(i) << std::endl;
			}
		}

		if (VERBOSE) {
			std::cout << "New particle collided at [" << collision->position[0] << ", " << collision->position[1] << ", " << collision->position[2] << "] ";
			if (collision->idx > -1) {
				std::cout << "with atom at [" << atoms[collision->idx][0] << ", " << atoms[collision->idx][1] << ", " << atoms[collision->idx][2] << "] " << std::endl;
			}
			else {
				std::cout << std::endl;
			}
		}

		// DIFFUSION
		auto startd = std::chrono::high_resolution_clock::now();
		//std::array<float, 3> force = { 0, 0, 0 };
		float force_mag = 0;
		float distance = diffusion_length;
		float step_size = 0.01;
		/*while (distance > 0) {
			direction = { 0, 0, 0 };
			std::vector<int>* neighbors = cubes.find_nearest_bin(collision->position);
			if (VERBOSE) {
				std::cout << "Found " << neighbors->size() << " neighbors";
			}
			for (auto p : *neighbors) {
				float dist2 = pow(atoms[p][0] - collision->position[0], 2) + pow(atoms[p][1] - collision->position[1], 2) + pow(atoms[p][2] - collision->position[2], 2);
				float r = sqrt(dist2);
				float r6 = pow(dist2, 3);
				float r12 = pow(r6, 2);
				force_mag =  (2/r) *  (2*s12 / r12 - s6 / r6);
				direction[0] += (collision->position[0] - atoms[p][0]) / r * force_mag;
				direction[1] += (collision->position[1] - atoms[p][1]) / r * force_mag;
				direction[2] += (collision->position[2] - atoms[p][2]) / r * force_mag;
			}
			if (collision->position[2] < 2) {
				float r = collision->position[2];
				float r6 = pow(r, 6);
				float r12 = pow(r6, 2);
				force_mag = (2 / r) * (2 * s12 / r12 - s6 / r6);
				direction[2] += force_mag;
			}
			force_mag = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);

			if (VERBOSE) {
				std::cout << "\tForce: " << force_mag << "\tDirection: [" << direction[0] << ", " << direction[1] << ", " << direction[2] << "]; ";
			}
			
			if (force_mag != 0) {
				if (force_mag < 1) {
					collision->position[0] += step_size * direction[0];
					collision->position[1] += step_size * direction[1];
					collision->position[2] += step_size * direction[2];
					if (VERBOSE) {
						std::cout << "Traveled " << step_size * force_mag << " nm ";
					}
				}
				else {
					collision->position[0] += step_size * direction[0] / force_mag;
					collision->position[1] += step_size * direction[1] / force_mag;
					collision->position[2] += step_size * direction[2] / force_mag;
					if (VERBOSE) {
						std::cout << "Traveled " << step_size << " nm ";
					}
				}
				if (collision->position[0] < 0) {
					collision->position[0] += L;
				}
				else if (collision->position[0] > L) {
					collision->position[0] -= L;
				}
				if (collision->position[1] < 0) {
					collision->position[1] += L;
				}
				else if (collision->position[1] > L) {
					collision->position[1] -= L;
				}
				if (VERBOSE) {
					std::cout << "to [" << collision->position[0] << ", " << collision->position[1] << ", " << collision->position[2] << "] " << std::endl;
				}
			}
			else if (VERBOSE) {
				std::cout << std::endl;
			}
			distance -= step_size;
		}*/

		while (distance > 0) {
			// first jump near minimum
			current_minimum = collision->position;
			float remaining_distance = (distance > (*radii)[0]) + 1/length_scale/2 ? (*radii)[0] : distance - 1/length_scale/2;
			potentials.find_local_minimum(current_minimum, remaining_distance);

			for (float settle = 1 / length_scale/2; settle > 0; settle -= step_size) {
				direction = { 0, 0, 0 };
				std::vector<int>* neighbors = cubes.find_nearest_bin(current_minimum);
				if (VERBOSE) {
					std::cout << "Found " << neighbors->size() << " neighbors";
				}
				for (auto p : *neighbors) {
					float dist2 = pow(atoms[p][0] - current_minimum[0], 2) + pow(atoms[p][1] - current_minimum[1], 2) + pow(atoms[p][2] - current_minimum[2], 2);
					float r = sqrt(dist2);
					float r6 = pow(dist2, 3);
					float r12 = pow(r6, 2);
					force_mag = (2 / r) * (2 * s12 / r12 - s6 / r6);
					direction[0] += (current_minimum[0] - atoms[p][0]) / r * force_mag;
					direction[1] += (current_minimum[1] - atoms[p][1]) / r * force_mag;
					direction[2] += (current_minimum[2] - atoms[p][2]) / r * force_mag;
				}
				if (collision->position[2] < 2) {
					float r = collision->position[2];
					float r6 = pow(r, 6);
					float r12 = pow(r6, 2);
					force_mag = (2 / r) * (2 * s12 / r12 - s6 / r6);
					direction[2] += force_mag;
				}
				force_mag = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);

				if (VERBOSE) {
					std::cout << "\tForce: " << force_mag << "\tDirection: [" << direction[0] << ", " << direction[1] << ", " << direction[2] << "]; ";
				}

				if (force_mag != 0) {
					if (force_mag < 1) {
						current_minimum[0] += step_size * direction[0];
						current_minimum[1] += step_size * direction[1];
						current_minimum[2] += step_size * direction[2];
						if (VERBOSE) {
							std::cout << "Traveled " << step_size * force_mag << " nm ";
						}
					}
					else {
						current_minimum[0] += step_size * direction[0] / force_mag;
						current_minimum[1] += step_size * direction[1] / force_mag;
						current_minimum[2] += step_size * direction[2] / force_mag;
						if (VERBOSE) {
							std::cout << "Traveled " << step_size << " nm ";
						}
					}
					if (current_minimum[0] < 0) {
						current_minimum[0] += L;
					}
					else if (current_minimum[0] > L) {
						current_minimum[0] -= L;
					}
					if (current_minimum[1] < 0) {
						current_minimum[1] += L;
					}
					else if (current_minimum[1] > L) {
						current_minimum[1] -= L;
					}
					if (VERBOSE) {
						std::cout << "to [" << current_minimum[0] << ", " << current_minimum[1] << ", " << current_minimum[2] << "] " << std::endl;
					}

					if (isnan(current_minimum[0])) {
						break;
					}
				}
				else if (VERBOSE) {
					std::cout << std::endl;
				}
			}

			// quit early if we're not moving
			if (current_minimum[0] == collision->position[0] && current_minimum[1] == collision->position[1] && current_minimum[2] == collision->position[2]) {
				break;
			}
			else {
				collision->position[0] == current_minimum[0];
				collision->position[1] == current_minimum[1];
				collision->position[2] == current_minimum[2];
			}
			distance -= remaining_distance + 1 / length_scale/2;
		}
		

		if (VERBOSE) {
			std::cout << "Particle stuck at [" << collision->position[0] << ", " << collision->position[1] << ", " << collision->position[2] << "] " << std::endl;
		}
		
		
		/*if (diffusion_length > 0) {
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

		}*/
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
		corridors.add_to_bins(&collision->position, (*radii)[sp], n);
		cubes.add_to_bins(n, collision->position);

		// Update potentials
		if (diffusion_length > 0) {
			potentials.scaled_point(collision->position, center);
			for (int k = -diameter / 2; k < diameter / 2; k++) {
				int kk = center[2] + k;
				if (kk < 0 || kk > potentials.get_Hs()) {
					continue;
				}
				for (int j = -diameter / 2; j < diameter / 2; j++) {
					int jj = center[1] + j;
					for (int i = -diameter / 2; i < diameter / 2; i++) {
						int ii = center[0] + i;

						// add to potential field
						potentials(ii, jj, kk) += region(i + diameter / 2, j + diameter / 2, k + diameter / 2);
					}
				}
			}
		}

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

	// Update object for parameters
	ContinuousSimulationParameters * layer_params = new ContinuousSimulationParameters();
	layer_params->length = L;
	layer_params->width = L;
	layer_params->height = H;
	layer_params->theta = theta;
	layer_params->species = species;
	layer_params->radii = radii;
	layer_params->diffusion_length = diffusion_length;
	layer_params->repetitions = reps;
	layer_params->system = system;
	layer_params->seed = seed;
	layer_params->bin_size = bin_size;
	layer_params->cube_size = 2.0f;
	layer_params->time_taken = time;
	layer_params->time_finished = epoch_time;
	params->addLayer(layer_params);
	params->serialize();

	// Create filename
	std::string len_str = std::to_string(L);
	len_str.erase(len_str.find_last_not_of('0') + 1, std::string::npos);
	len_str.erase(len_str.find_last_not_of('.') + 1, std::string::npos);
	std::string theta_str = std::to_string(theta);
	theta_str.erase(theta_str.find_last_not_of('0') + 1, std::string::npos);
	theta_str.erase(theta_str.find_last_not_of('.') + 1, std::string::npos);
	std::string diff_str = std::to_string(diffusion_length);
	diff_str.erase(diff_str.find_last_not_of('0') + 1, std::string::npos);
	diff_str.erase(diff_str.find_last_not_of('.') + 1, std::string::npos);
	std::string filename = "structures/cts/STF_" + system + "_L" + len_str + "_Th" + theta_str + "_D" + diff_str + "_N" + std::to_string(params->deposited) + "_" + std::to_string(epoch_time);

	// Save objects
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

	//potentials.save_file("potential.npy");
	//region.save_file("region.npy");
	corridors.save_file("priority.npy");

	std::ofstream json_file;
	json_file.open("params.json");
	json_file << params->serialization;
	json_file.close();

	// Combine the files into one zip/npz file
	zipFile zf = zipOpen64((filename + ".simc").c_str(), 0);
	zipClose(zf, NULL);
	int err = writeFileToZipCTS((filename + ".simc").c_str(), "grid.npy");
	//err = writeFileToZip((filename + ".sim").c_str(), "structures/cts/diff.npy");
	err = writeFileToZipCTS((filename + ".simc").c_str(), "params.json");
	//err = writeFileToZipCTS((filename + ".simc").c_str(), "potential.npy");
	err = writeFileToZipCTS((filename + ".simc").c_str(), "priority.npy");

	// Remove the individual files
	// TODO: combine the writeFileToZip function and the writing of npy files so clean up isn't necessary
	std::remove("grid.npy");
	//std::remove("structures/cts/diff.npy");
	std::remove("params.json");
	//std::remove("potential.npy");
	std::remove("priority.npy");

	// Print the timing
	std::cout << reps << " reps completed in " << time << " seconds; \n" << reps / time << " reps per second." << std::endl;
	std::cout << "Line finding and traversal were " << timel << " seconds of that time." << std::endl;
	std::cout << "Diffusion was " << timed << " seconds of that time." << std::endl;

	//free(potentials);
	return reps;
}