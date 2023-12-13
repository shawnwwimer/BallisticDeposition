#include "SimulationContinuousHashedSpace.h"
#define VERBOSE (false)

int obliqueDepositionContinuousHashed(float theta, float L, float H, uint32_t reps, float bin_size, uint32_t seed, float diffusion_length, float length_scale, std::vector<int8_t>* species, std::vector<float>* radii, std::vector<float>* spread, std::vector<std::vector<float>>* weights, std::vector<std::vector<float>> inputGrid, ContinuousSimulationParametersFull* params, std::string& system, DiffusionMethod diffusion_method, FilesToSave* save)
{
	switch (diffusion_method) {
	case DiffusionMethod::PotentialHoppingLUT:
		std::cout << "Deposition at " << theta << " and " << diffusion_length << " nm diffusion length by hopping." << std::endl;
		break;
	case DiffusionMethod::HopAndSettleLUT:
		std::cout << "Deposition at " << theta << " and " << diffusion_length << " nm diffusion length by hopping and settling." << std::endl;
		break;
	case DiffusionMethod::ForcePushingLUT:
		std::cout << "Deposition at " << theta << " and " << diffusion_length << " nm diffusion length by settling." << std::endl;
		break;
	case DiffusionMethod::NumericalMinimization:
		std::cout << "Deposition at " << theta << " and " << diffusion_length << " nm diffusion length by numerical minimization." << std::endl;
		break;
	default:
		std::cout << "Deposition at " << theta << " and " << diffusion_length << " nm diffusion length by unknown method." << std::endl;
		break;
	}

	auto start = std::chrono::high_resolution_clock::now();
	std::vector<std::vector<float>> atoms;
	std::vector<std::vector<float>> landing_positions;
	std::vector<std::array<float, 3>> dests;

	float cell_size = 0.5f;
	SpaceHashMap hash_map = SpaceHashMap(L, H, &atoms, cell_size);

	OverlappingCubicSpacePartition cubes = OverlappingCubicSpacePartition(L, H, &atoms, 2.0);
	float s = 2 * (*radii)[0] / pow(2, 1.0 / 6.0);
	float s6 = pow(s, 6);
	float s12 = pow(s6, 2);

	float A = 1;
	float R = 1;
	float Vm = length_scale * (2 * (*radii)[0]);
	float zero_pot = Vm / pow(2, 1.0 / 6.0);

	Matrix3DLateralPBC region = Matrix3DLateralPBC(2.f + 1 / length_scale, 2.f + 1 / length_scale, 2.f + 1 / length_scale, length_scale);
	Matrix3DLateralPBC potentials = Matrix3DLateralPBC(L, L, H, length_scale);
	int diameter = region.get_Hs();
	if (diffusion_length > 0 && (diffusion_method == DiffusionMethod::PotentialHoppingLUT || diffusion_method == DiffusionMethod::HopAndSettleLUT)) {
		region.initialize();
		potentials.initialize();
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

						if (r == 0) {
							region(i, j, k) = 1e6;
						}
					}
				}
			}
		}
		for (int k = 0; k < (diameter - 1) / 2; k++) {
			float mag = A * pow((zero_pot / (k + 1)), 12) - R * pow((zero_pot / (k + 1)), 6);
			for (int i = 0; i < potentials.get_Ls(); i++) {
				for (int j = 0; j < potentials.get_Ws(); j++) {
					potentials(i, j, k) = mag;
				}
			}
		}
	}

	//float* region = (float* )malloc(sizeof(float) * diameter * diameter * diameter * 4);

	int steps = round(length_scale * diffusion_length);
	int rads = round(length_scale * 2 * (*radii)[0]);
	float scaled_position[3] = { 0 };
	std::array<float, 3> direction = { 0, 0, 0 };
	std::array<float, 3> current_minimum = { 0, 0, 0 };
	std::array<int, 3> center = { 0, 0, 0 };

	// Populate from input grid
	float current_fiber_id = 0;
	int n = 0;
	if (inputGrid.size() > 0) {
		for (auto atom : inputGrid) {
			if (atom[5] + 1 > current_fiber_id) {
				current_fiber_id = atom[5] + 1;
			}
			std::array<float, 3> position = { atom[0], atom[1], atom[2] };
			atoms.push_back(atom);
			hash_map.add_to_bin(&position, n);
			cubes.add_to_bins(n, &position);
			n++;
		}
	}

	// Create random number generator
	if (seed == 0) {
		std::random_device rd;
		seed = rd();
	}
	std::mt19937 gen(seed);
	std::uniform_real_distribution<float> dist(0, L);

	std::mt19937 gens(seed);
	std::uniform_int_distribution<> sist(0, species->size());

	std::mt19937 gend(seed + 1);
	std::uniform_real_distribution<float> dist_d(0, 1);

	for (int n = 0; n < reps; n++) {
		dests.push_back({ dist(gen), dist(gen), 0 });
		//dests.push_back({ 0, 1+(float)n/10, 0 });
	}

	// Print how long it took
	std::chrono::duration<double, std::milli> timep = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Set-up took " << timep.count() / 1000 << " s." << std::endl;

	// Begin
	start = std::chrono::high_resolution_clock::now();
	int update = (reps > 128) ? 128 : 1; // number of printed updates
	std::array<float, 3> dest;
	std::array<float, 3> src;
	std::array<float, 3> landing_position;
	std::vector<float> new_atom;
	std::vector<int> bins;
	float fiber = 0;
	double timel = 0;
	double timed = 0;
	float x1, y1, z1;
	collision_description* collision = new collision_description(-1, -1, -1, -1);
	int collision_idx = -1;
	float earliest_collision = 1e6;
	for (int n = 0; n < reps; n++) {

		auto startl = std::chrono::high_resolution_clock::now();

		// Generate 
		dest = dests[n];
		//dest = { 32, 0.25f+0.05f*n, 0 };//
		//dest = { target[n], 1, 0 };

		// Choose species
		int sp = 0;

		// Drop particle
		float slant = tan(theta * M_PI / 180.f - M_PI / 2.f);
		src[0] = ((H - 1) / slant + dest[0]);
		src[1] = dest[1];
		src[2] = H - 1;

		float dx = abs(dest[0] - src[0]), dy = abs(dest[1] - src[1]), dz = abs(dest[2] - src[2]);
		float x = fmod(src[0], L), y = src[1], z = src[2];
		float realz = z;
		if (x < 0) {
			x += L;
		}

		// get bin by modified Bresenham
		float p = 2 * dz * cell_size - dx;
		while (z > -cell_size) {
			// increase x by cell_size and restrict within L
			x = x + cell_size;
			realz += slant * cell_size;
			if (x > L) {
				x -= L;
			}

			// adjust p and z if necessary
			if (p >= 0) {
				z -= cell_size;
				if (z == -cell_size) {
					break;
				}
				p -= 2 * dx * cell_size;
			}

			// Get bins on same level
			int idx = hash_map.return_bin_idx(x, y, realz);
			float ydiff = fmodf(y, cell_size);
			if (ydiff < (*radii)[0] * 2) {
				if (y < cell_size) {
					bins.push_back(hash_map.return_bin_idx(x, L - 0.01f, realz));
				}
				else {
					bins.push_back(hash_map.return_bin_idx(x, y - cell_size, realz));
				}
			}
			if (ydiff > cell_size - (*radii)[0] * 2) {
				if (y > L - cell_size) {
					bins.push_back(hash_map.return_bin_idx(x, 0.01f, realz));
				}
				else {
					bins.push_back(hash_map.return_bin_idx(x, y + cell_size, realz));
				}
			}
			bins.push_back(idx);
			float zdiff = fmodf(realz, cell_size);//(-p / 1 / dx) * cell_size;
			if (cell_size-zdiff < (*radii)[0] * 2 && realz+cell_size < H) {
				// near the upper cell
				bins.push_back(hash_map.return_bin_idx(x, y, realz + cell_size));
				if (ydiff < (*radii)[0] * 2) {
					if (y < cell_size) {
						bins.push_back(hash_map.return_bin_idx(x, L-0.01f, realz + cell_size));
					}
					else {
						bins.push_back(hash_map.return_bin_idx(x, y-cell_size, realz + cell_size));
					}
				}
				if (ydiff > cell_size - (*radii)[0] * 2) {
					if (y > L - cell_size) {
						bins.push_back(hash_map.return_bin_idx(x, 0.01f, realz + cell_size));
					}
					else {
						bins.push_back(hash_map.return_bin_idx(x, y + cell_size, realz + cell_size));
					}
				}
			}
			if (zdiff < (*radii)[0] * 2 && realz > cell_size) {
				bins.push_back(idx - hash_map.bins_in_layer);
				if (ydiff < (*radii)[0] * 2) {
					if (y < cell_size) {
						bins.push_back(hash_map.return_bin_idx(x, L - 0.01f, realz - cell_size));
					}
					else {
						bins.push_back(hash_map.return_bin_idx(x, y - cell_size, realz - cell_size));
					}
				}
				if (ydiff > cell_size - (*radii)[0] * 2) {
					if (y > L - cell_size) {
						bins.push_back(hash_map.return_bin_idx(x, 0.01f, realz - cell_size));
					}
					else {
						bins.push_back(hash_map.return_bin_idx(x, y + cell_size, realz - cell_size));
					}
				}
			}
			// increment error appropriately
			p += 2 * dz * cell_size;

			// now determine potential collisions
			// iterate over bins
			for (int i : bins) {
				// iterate over all atoms in the bins
				for (int j : hash_map.bins[i]) {
					// sum of radii
					float rs = (*radii)[0] + atoms[j][4];

					// quadratic equation
					float a = -(1 * 1 + 0 * 0 + slant * slant);
					float b = -2 * (1 * (x - atoms[j][0]) + 0 * (y - atoms[j][1]) + slant * (realz - atoms[j][2]));
					float c = -(pow(x - atoms[j][0], 2) + pow(y - atoms[j][1], 2) + pow(realz - atoms[j][2], 2) - pow(rs, 2));
					float d = b * b - 4 * a * c;

					if (d >= 0) {
						d = pow(d, 0.5);
						// get variable for parameterized position
						float t1 = (-b - d) / (2 * a);
						float t2 = (-b + d) / (2 * a);
						float t = t1 < t2 ? t1 : t2;

						// add collision
						collision_idx = j;
						if (t < earliest_collision) {
							earliest_collision = t;
							collision->position = { modulof(x + 1 * t, L), y + 0 * t, realz + slant * t };
							collision->idx = j;
						}
					}
				}
			}
			bins.clear();

			if (earliest_collision < 1e-6) {
				break;
			}
		}

		// check that collision wouldn't occur under the surface
		if (collision_idx == -1 || collision->position[2] < (*radii)[0]) {
			collision->position = { modulof(dest[0] + (*radii)[0] / slant, L), dest[1], (*radii)[0] };
			collision->idx = -1;
		}
		collision_idx = -1;
		earliest_collision = 1e6;
		//collision_description* collision = corridors.drop_particle(&dest, (*radii)[sp], &atoms);

		std::chrono::duration<double, std::milli> dural = std::chrono::high_resolution_clock::now() - startl;
		timel += dural.count() / 1000;

		if (VERBOSE) {
			std::cout << "New particle collided at [" << collision->position[0] << ", " << collision->position[1] << ", " << collision->position[2] << "] ";
			if (collision->idx > -1) {
				std::cout << "with atom at [" << atoms[collision->idx][0] << ", " << atoms[collision->idx][1] << ", " << atoms[collision->idx][2] << "] " << std::endl;
			}
			else {
				std::cout << std::endl;
			}
		}

		landing_position = collision->position;

		//if (n == 50632) {
		//	printf("jesus");
		//}

		// DIFFUSION
		auto startd = std::chrono::high_resolution_clock::now();
		//std::array<float, 3> force = { 0, 0, 0 };
		float force_mag = 0;
		float distance = diffusion_length;
		float step_size = 0.01;

		// this branch is for diffusion that uses a LUT for the potential of a test particle
		if (diffusion_method == DiffusionMethod::PotentialHoppingLUT || diffusion_method == DiffusionMethod::HopAndSettleLUT) {
			while (distance > 0) {
				// start from current position
				current_minimum = collision->position;
				float remaining_distance = distance > 3 * (*radii)[0] ? 3 * (*radii)[0] : distance;
				if (current_minimum[2] < 0) {
					current_minimum[2] = (*radii)[0];
				}
				// first jump near minimum
				potentials.find_local_minimum(current_minimum, remaining_distance, dist_d(gend));

				if (diffusion_method == DiffusionMethod::HopAndSettleLUT) {
					for (float settle = 1 / length_scale / 2; settle > 0; settle -= step_size) {
						direction = { 0, 0, 0 };
						std::vector<int>* neighbors = cubes.find_nearest_bin(&current_minimum);
						if (VERBOSE) {
							std::cout << "Found " << neighbors->size() << " neighbors";
						}
						for (auto p : *neighbors) {
							float dx = current_minimum[0] - atoms[p][0];
							float dy = current_minimum[1] - atoms[p][1];
							if (dx > L / 2) {
								dx -= L;
							}
							else if (dx < -L / 2) {
								dx += L;
							}
							if (dy > L / 2) {
								dy -= L;
							}
							else if (dy < -L / 2) {
								dy += L;
							}
							float dist2 = pow(dx, 2) + pow(dy, 2) + pow(atoms[p][2] - current_minimum[2], 2);
							float r = sqrt(dist2);
							float r6 = pow(dist2, 3);
							float r12 = pow(r6, 2);
							force_mag = -(2 / r) * (s6 / r6 - 2 * s12 / r12);
							direction[0] += dx / r * force_mag;
							direction[1] += dy / r * force_mag;
							direction[2] += (current_minimum[2] - atoms[p][2]) / r * force_mag;
						}
						if (current_minimum[2] < (*radii)[0] && direction[2] < 0) {
							direction = { 0, 0, 1 };
						}
						else if (current_minimum[2] < 2) {
							float r = current_minimum[2] + (*radii)[0];
							float r6 = pow(r, 6);
							float r12 = pow(r6, 2);
							force_mag = (2 / r) * (s6 / r6 - 2 * s12 / r12);
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
						}
						else if (VERBOSE) {
							std::cout << std::endl;
						}
					}
				}

				// quit early if we're not moving
				if (current_minimum[0] == collision->position[0] && current_minimum[1] == collision->position[1] && current_minimum[2] == collision->position[2]) {
					break;
				}
				else if (isnan(current_minimum[0]) || current_minimum[2] < (*radii)[0]) {
					break;
				}
				else {
					collision->position[0] = current_minimum[0];
					collision->position[1] = current_minimum[1];
					collision->position[2] = current_minimum[2];
				}
				distance -= remaining_distance;
			}
		}
		// this branch is for diffusion that directly implements the force without a LUT
		else if (diffusion_method == DiffusionMethod::ForcePushingLUT) {
			while (distance > 0) {
				current_minimum = collision->position;
				direction = { 0, 0, 0 };
				std::vector<int>* neighbors = cubes.find_nearest_bin(&current_minimum);
				if (VERBOSE) {
					std::cout << "Found " << neighbors->size() << " neighbors";
				}
				for (auto p : *neighbors) {
					float dist2 = pow(atoms[p][0] - collision->position[0], 2) + pow(atoms[p][1] - collision->position[1], 2) + pow(atoms[p][2] - collision->position[2], 2);
					float r = sqrt(dist2);
					float r6 = pow(dist2, 3);
					float r12 = pow(r6, 2);
					force_mag = (2 / r) * (2 * s12 / r12 - s6 / r6);
					direction[0] += (collision->position[0] - atoms[p][0]) / r * force_mag;
					direction[1] += (collision->position[1] - atoms[p][1]) / r * force_mag;
					direction[2] += (collision->position[2] - atoms[p][2]) / r * force_mag;
				}
				if (current_minimum[2] < 2) {
					float r = current_minimum[2];
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

					if (isnan(current_minimum[0]) || current_minimum[2] < (*radii)[0]) {
						break;
					}
				}
				else if (VERBOSE) {
					std::cout << std::endl;
				}

				// quit early if we're not moving
				if (current_minimum[0] == collision->position[0] && current_minimum[1] == collision->position[1] && current_minimum[2] == collision->position[2]) {
					break;
				}
				else {
					collision->position[0] = current_minimum[0];
					collision->position[1] = current_minimum[1];
					collision->position[2] = current_minimum[2];
				}
				distance -= step_size;
			}
		}
		// this branch is for diffusion that numerically minimizes the potential
		else if (diffusion_method == DiffusionMethod::NumericalMinimization) {
			while (distance > 0) {
				// start from current position
				current_minimum = collision->position;
				float remaining_distance = distance > 3 * (*radii)[0] ? 3 * (*radii)[0] : distance;
				if (current_minimum[2] < 0) {
					current_minimum[2] = (*radii)[0];
				}
				std::vector<double>* minimum = cubes.find_local_minimum(&current_minimum, distance);
				collision->position = { modulof((float)(*minimum)[0], L), modulof((float)(*minimum)[1], L), (float)(*minimum)[2] };
				distance -= remaining_distance;
			}
		}

		if (VERBOSE) {
			std::cout << "Particle stuck at [" << collision->position[0] << ", " << collision->position[1] << ", " << collision->position[2] << "] " << std::endl;
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
		//corridors.add_to_bins(&collision->position, (*radii)[sp], n + inputGrid.size());
		hash_map.add_to_bin(&collision->position, n + inputGrid.size());
		landing_positions.push_back({ landing_position[0], landing_position[1], landing_position[2] });

		if (diffusion_length > 0 && (diffusion_method == DiffusionMethod::HopAndSettleLUT || diffusion_method == DiffusionMethod::ForcePushingLUT || diffusion_method == DiffusionMethod::NumericalMinimization)) {
			cubes.add_to_bins(n + inputGrid.size(), &collision->position);
		}

		// Update potentials
		if (diffusion_length > 0 && (diffusion_method == DiffusionMethod::PotentialHoppingLUT || diffusion_method == DiffusionMethod::HopAndSettleLUT)) {
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
	ContinuousSimulationParameters* layer_params = new ContinuousSimulationParameters();
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
	layer_params->length_scale = length_scale;
	layer_params->cube_size = 2.0f;
	layer_params->diffusion_method = (int)diffusion_method;
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
	zipFile zf = zipOpen64((filename + ".simc").c_str(), 0);
	zipClose(zf, NULL);
	int err;

	// grid
	float* outGrid = (float*)malloc(sizeof(float) * atoms.size() * 6);
	if (outGrid != NULL) {
		for (int j = 0; j < atoms.size(); j++) {
			for (int i = 0; i < 6; i++) {
				outGrid[j * 6 + i] = atoms[j][i];
			}
		}
		cnpy::npy_save("grid.npy", outGrid, { reps, 6 });
		free(outGrid);
	}
	err = writeFileToZipCTS((filename + ".simc").c_str(), "grid.npy");
	if (err == ZIP_ERRNO) {
		std::cout << "Couldn't add to zip file correctly." << std::endl;
	}
	std::remove("grid.npy");

	// priority
	//if (save != nullptr && save->priority) {
	//	corridors.save_file("priority.npy");
	//	err = writeFileToZipCTS((filename + ".simc").c_str(), "priority.npy");
	//	if (err == ZIP_ERRNO) {
	//		std::cout << "Couldn't add to zip file correctly." << std::endl;
	//	}
	//	std::remove("priority.npy");
	//}

	// collision positions
	if (save == nullptr || save->collisions) {
		float* colout = (float*)malloc(sizeof(float) * reps * 3);
		for (int i = 0; i < reps; i++) {
			for (int j = 0; j < 3; j++) {
				colout[i * 3 + j] = landing_positions[i][j];
			}
		}
		cnpy::npy_save("collisions.npy", colout, { reps, 3 });
		free(colout);
		err = writeFileToZipCTS((filename + ".simc").c_str(), "collisions.npy");
		if (err == ZIP_ERRNO) {
			std::cout << "Couldn't add to zip file correctly." << std::endl;
		}
		std::remove("collisions.npy");
	}

	// desinations
	if (save != nullptr && save->destinations) {
		float* destout = (float*)malloc(sizeof(float) * reps * 3);
		for (int i = 0; i < reps; i++) {
			for (int j = 0; j < 3; j++) {
				destout[i * 3 + j] = dests[i][j];
			}
		}
		cnpy::npy_save("destinations.npy", destout, { reps, 3 });
		free(destout);
		err = writeFileToZipCTS((filename + ".simc").c_str(), "destinations.npy");
		if (err == ZIP_ERRNO) {
			std::cout << "Couldn't add to zip file correctly." << std::endl;
		}
		std::remove("destinations.npy");
	}

	// Parameters
	std::ofstream json_file;
	json_file.open("params.json");
	json_file << params->serialization;
	json_file.close();
	err = writeFileToZipCTS((filename + ".simc").c_str(), "params.json");
	if (err == ZIP_ERRNO) {
		std::cout << "Couldn't add to zip file correctly." << std::endl;
	}
	std::remove("params.json");
	params->clearLayers();

	// volume potential
	if (save != nullptr && save->volume_potential) {
		if (potentials.save_file("potential.npy")) {
			err = writeFileToZipCTS((filename + ".simc").c_str(), "potential.npy");
			if (err == ZIP_ERRNO) {
				std::cout << "Couldn't add to zip file correctly." << std::endl;
			}
			std::remove("potential.npy");
		}
	}

	// atomic potential
	if (save != nullptr && save->atomic_potential) {
		if (region.save_file("region.npy")) {
			err = writeFileToZipCTS((filename + ".simc").c_str(), "region.npy");
			if (err == ZIP_ERRNO) {
				std::cout << "Couldn't add to zip file correctly." << std::endl;
			}
			std::remove("region.npy");
		}
	}

	// Print the timing
	std::cout << reps << " reps completed in " << time << " seconds; \n" << reps / time << " reps per second." << std::endl;
	std::cout << "Line finding and traversal were " << timel << " seconds of that time." << std::endl;
	std::cout << "Diffusion was " << timed << " seconds of that time." << std::endl;

	return reps;
}