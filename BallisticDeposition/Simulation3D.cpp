#include "Simulation3D.h"
#include "cnpy.h"
//#include "math_utils.h"

uint16_t Bresenham3D(uint16_t* src, uint16_t* dest, uint16_t* points) {
	// Distance
	uint16_t d[3] = { (uint16_t)abs((int)dest[0] - (int)src[0]), (uint16_t)abs((int)dest[1] - (int)src[1]), (uint16_t)abs((int)dest[2] - (int)src[2]) };

	// Initialize 
	uint16_t* row = &points[0];
	row[0] = src[0];
	row[1] = src[1];
	row[2] = src[2];
	uint16_t x1 = src[0];
	uint16_t x2 = dest[0];
	uint16_t y1 = src[1];
	uint16_t y2 = dest[1];
	uint16_t z1 = src[2];
	uint16_t z2 = dest[2];

	// Determine directions
	int xs, ys, zs;
	if (dest[0] > src[0]) {
		xs = 1;
	}
	else {
		xs = -1;
	}
	if (dest[1] > src[1]) {
		ys = 1;
	}
	else {
		ys = -1;
	}
	if (dest[2] > src[2]) {
		zs = 1;
	}
	else {
		zs = -1;
	}

	// Row counter
	int i = 1;

	int32_t p1, p2;
	if (d[0] >= d[1] && d[0] > d[2]) {
		// Driving axis is x-axis
		p1 = 2 * d[1] - d[0];
		p2 = 2 * d[2] - d[0];
		while (x1 != x2) {
			x1 += xs;
			// Update space
			if (p1 >= 0) {
				y1 += ys;
				p1 -= 2 * d[0];
			}
			if (p2 >= 0) {
				z1 += zs;
				p2 -= 2 * d[0];
			}
			p1 += 2 * d[1];
			p2 += 2 * d[2];

			// Add to points
			row = &points[i * 3];
			row[0] = x1;
			row[1] = y1;
			row[2] = z1;
			i += 1;
		}
	}
	else if (d[1] >= d[0] && d[1] > d[2]) {
		// Driving axis is y-axis
		p1 = 2 * d[0] - d[1];
		p2 = 2 * d[2] - d[1];
		while (y1 != y2) {
			y1 += ys;
			// Update space
			if (p1 >= 0) {
				x1 += xs;
				p1 -= 2 * d[1];
			}
			if (p2 >= 0) {
				z1 += zs;
				p2 -= 2 * d[1];
			}
			p1 += 2 * d[0];
			p2 += 2 * d[2];

			// Add to points
			row = &points[i * 3];
			row[0] = x1;
			row[1] = y1;
			row[2] = z1;
			i += 1;
		}
	}
	else {
		// Driving axis is y-axis
		p1 = 2 * d[0] - d[2];
		p2 = 2 * d[1] - d[2];
		while (z1 != z2) {
			z1 += zs;
			// Update space
			if (p1 >= 0) {
				x1 += xs;
				p1 -= 2 * d[2];
			}
			if (p2 >= 0) {
				y1 += ys;
				p2 -= 2 * d[2];
			}
			p1 += 2 * d[0];
			p2 += 2 * d[1];

			// Add to points
			row = &points[i * 3];
			row[0] = x1;
			row[1] = y1;
			row[2] = z1;
			i += 1;
		}
	}

	// Return the number of points involved
	return i;
}

uint16_t* traversePathCalculated(uint16_t* points, uint16_t points_len, uint16_t* dest, uint16_t* collision, int8_t* grid, uint16_t L, uint16_t H) {
	for (int i = 0; i < points_len; i++) {
		uint16_t* row = &points[i * 3];
		if (grid[row[2]*H*L + ((row[1] + dest[1]) % L) * L + ((row[0] + dest[0]) % L)] > 0) {
			// Intersection before reaching end destination
			if ((i == 0) || (points[(i - 1)*3 + 2] >= H)) {
				// Can't assign outside of grid
				return nullptr;
			}
			row = &points[(i - 1) * 3];
			collision[0] = row[0] % L;
			collision[1] = row[1] % L;
			collision[2] = row[2];
			return collision;
		}
		else if (row[0] == 0 && row[1] == 0 && row[2] == 0) {
			// Lands at destination
			collision[0] = dest[0];
			collision[1] = dest[1];
			collision[2] = dest[2];
			return collision;
		}
	}

	return nullptr;
}

uint16_t traversePathRealTime(int32_t* src, int32_t* dest, uint16_t* collision, int8_t* grid, uint16_t L, uint16_t H) {
	// Distance
	int32_t d[3] = { abs(dest[0] - src[0]), abs(dest[1] - src[1]), abs(dest[2] - src[2]) };

	// Initialize 
	int x1 = src[0];
	int x2 = dest[0];
	int y1 = src[1];
	int y2 = dest[1];
	int z1 = src[2];
	int z2 = dest[2];

	// Determine directions
	int xs, ys, zs;
	if (dest[0] > src[0]) {
		xs = 1;
	}
	else {
		xs = -1;
	}
	if (dest[1] > src[1]) {
		ys = 1;
	}
	else {
		ys = -1;
	}
	if (dest[2] > src[2]) {
		zs = 1;
	}
	else {
		zs = -1;
	}

	// Row counter
	int i = 1;

	int32_t p1, p2;
	if (d[0] >= d[1] && d[0] > d[2]) {
		// Driving axis is x-axis
		p1 = 2 * d[1] - d[0];
		p2 = 2 * d[2] - d[0];
		int xlast = src[0];
		int ylast = src[1];
		int zlast = src[2];
		int xidx = (x1 % L + L) % L;
		int yidx = (y1 % L + L) % L;
		while (x1 != x2) {
			if (grid[z1 * L * L + yidx * L + xidx] > 0) {
				collision[0] = (xlast % L + L) % L;
				collision[1] = (ylast % L + L) % L;
				collision[2] = zlast;
				return i;
			}
			xlast = x1;
			ylast = y1;
			zlast = z1;
			x1 += xs;
			// Update space
			if (p1 >= 0) {
				y1 += ys;
				p1 -= 2 * d[0];
			}
			if (p2 >= 0) {
				z1 += zs;
				p2 -= 2 * d[0];
			}
			xidx = (x1 % L + L) % L;
			yidx = (y1 % L + L) % L;
			p1 += 2 * d[1];
			p2 += 2 * d[2];
			i += 1;
		}
		if (grid[z1 * L * L + yidx * L + xidx] > 0) {
			// collide just before surface
			collision[0] = (xlast % L + L) % L;
			collision[1] = (ylast % L + L) % L;
			collision[2] = zlast;
		}
		else {
			// Collide at surface
			collision[0] = xidx;
			collision[1] = yidx;
			collision[2] = z1;
		}
		
		i += 1;
	}
	else if (d[1] >= d[0] && d[1] > d[2]) {
		// Driving axis is y-axis
		p1 = 2 * d[0] - d[1];
		p2 = 2 * d[2] - d[1];
		int xlast = src[0];
		int ylast = src[1];
		int zlast = src[2];
		int xidx = (x1 % L + L) % L;
		int yidx = (y1 % L + L) % L;
		while (y1 != y2) {
			if (grid[z1 * L * L + yidx * L + xidx] > 0) {
				collision[0] = (xlast % L + L) % L;
				collision[1] = (ylast % L + L) % L;
				collision[2] = zlast;
				return i;
			}
			xlast = x1;
			ylast = y1;
			zlast = z1;
			y1 += ys;
			// Update space
			if (p1 >= 0) {
				x1 += xs;
				p1 -= 2 * d[1];
			}
			if (p2 >= 0) {
				z1 += zs;
				p2 -= 2 * d[1];
			}
			xidx = (x1 % L + L) % L;
			yidx = (y1 % L + L) % L;
			p1 += 2 * d[0];
			p2 += 2 * d[2];
			i += 1;
		}
		if (grid[z1 * L * L + yidx * L + xidx] > 0) {
			// collide just before surface
			collision[0] = (xlast % L + L) % L;
			collision[1] = (ylast % L + L) % L;
			collision[2] = zlast;
		}
		else {
			// Collide at surface
			collision[0] = xidx;
			collision[1] = yidx;
			collision[2] = z1;
		}
		i += 1;
	}
	else {
		// Driving axis is z-axis
		p1 = 2 * d[0] - d[2];
		p2 = 2 * d[1] - d[2];
		int xlast = src[0];
		int ylast = src[1];
		int zlast = src[2];
		int xidx = (x1 % L + L) % L;
		int yidx = (y1 % L + L) % L;
		while (z1 != z2) {
			if (grid[z1 * L * L + yidx * L + xidx] > 0) {
				collision[0] = (xlast % L + L) % L;
				collision[1] = (ylast % L + L) % L;
				collision[2] = zlast;
				return i;
			}
			xlast = x1;
			ylast = y1;
			zlast = z1;
			z1 += zs;
			// Update space
			if (p1 >= 0) {
				x1 += xs;
				p1 -= 2 * d[2];
			}
			if (p2 >= 0) {
				y1 += ys;
				p2 -= 2 * d[2];
			}
			xidx = (x1 % L + L) % L;
			yidx = (y1 % L + L) % L;
			p1 += 2 * d[0];
			p2 += 2 * d[1];
			i += 1;
		}
		if (grid[z1 * L * L + yidx * L + xidx % L] > 0) {
			// collide just before surface
			collision[0] = (xlast % L + L) % L;
			collision[1] = (ylast % L + L) % L;
			collision[2] = zlast;
		}
		else {
			// Collide at surface
			collision[0] = xidx;
			collision[1] = yidx;
			collision[2] = z1;
		}
		i += 1;
	}

	// Return the number of points involved
	return i;
}

uint32_t denseToSparse(int8_t* grid, int16_t** sparse, uint16_t L, uint16_t H) {
	// Count the number of points (and therefore rows) needed
	uint32_t points = 0;
	for (uint32_t k = 0; k < H; k++) {
		for (uint32_t j = 0; j < L; j++) {
			for (uint32_t i = 0; i < L; i++) {
				uint32_t idx = k * L * L + j * L + i;
				if (grid[idx] > 0) {
					points += 1;
				}
			}
		}
	}

	// Create sparse grid
	*sparse = (int16_t*)malloc(sizeof(int16_t) * points * 4);
	if (*sparse == nullptr) {
		return 0;
	}
	uint32_t counter = 0;
	for (uint32_t k = 0; k < H; k++) {
		for (uint32_t j = 0; j < L; j++) {
			for (uint32_t i = 0; i < L; i++) {
				uint32_t idx = k * L * L + j * L + i;
				if (grid[idx] > 0) {
					(*sparse)[counter * 4] = i;
					(*sparse)[counter * 4 + 1] = j;
					(*sparse)[counter * 4 + 2] = k;
					(*sparse)[counter * 4 + 3] = grid[idx];
					counter += 1;
				}
			}
		}
	}

	return points;
}

uint32_t denseToSparse(uint32_t* grid, uint32_t** sparse, uint16_t L, uint16_t H) {
	// Count the number of points (and therefore rows) needed
	uint32_t points = 0;
	for (uint32_t k = 0; k < H; k++) {
		for (uint32_t j = 0; j < L; j++) {
			for (uint32_t i = 0; i < L; i++) {
				uint32_t idx = k * L * L + j * L + i;
				if (grid[idx] > 0) {
					points += 1;
				}
			}
		}
	}

	// Create sparse grid
	*sparse = (uint32_t*)malloc(sizeof(uint32_t) * points * 4);
	if (*sparse == nullptr) {
		return 0;
	}
	uint32_t counter = 0;
	for (uint32_t k = 0; k < H; k++) {
		for (uint32_t j = 0; j < L; j++) {
			for (uint32_t i = 0; i < L; i++) {
				uint32_t idx = k * L * L + j * L + i;
				if (grid[idx] > 0) {
					(*sparse)[counter * 4] = i;
					(*sparse)[counter * 4 + 1] = j;
					(*sparse)[counter * 4 + 2] = k;
					(*sparse)[counter * 4 + 3] = grid[idx];
					counter += 1;
				}
			}
		}
	}

	return points;
}

uint16_t sparseToDense(int16_t* sparse, uint32_t num_points, int8_t* grid, uint16_t L, uint16_t H) {
	uint16_t maxh = 0;
	uint32_t counter = 0;
	int16_t row[4] = { 0, 0, 0, 0 };
	for (uint32_t i = 0; i < num_points; i++) {
		row[0] = sparse[i * 4];
		row[1] = sparse[i * 4 + 1];
		row[2] = sparse[i * 4 + 2];
		row[3] = sparse[i * 4 + 3];
		grid[row[2] * L * L + row[1] * L + row[0]] = row[3];
		if (row[2] > maxh) {
			maxh = row[2];
		}
	}
	return maxh;
}



int obliqueDeposition(float theta, uint16_t L, uint16_t H, uint32_t reps, float phi, float turns, uint32_t seed, uint16_t diffusion_steps, std::vector<int8_t> * species, std::vector<float> * spread, std::vector<std::vector<float>> * weights, int16_t* inputGrid, uint32_t inputGridPoints, int16_t** outGrid, uint32_t stepper_resolution, SimulationParametersFull* params, std::string &system) {
	auto start = std::chrono::high_resolution_clock::now();

	// Magic numbers
	double pi = 3.141592653;
	uint16_t update = 128;

	// Height ratcheting
	uint16_t maxh = 0;

	// Allocate memory for the grids we'll use
	int8_t* grid = (int8_t*)malloc(sizeof(int8_t) * L * L * H);
	uint32_t* ordered = (uint32_t*)malloc(sizeof(uint32_t) * L * L * H);
	for (uint32_t i = 0; i < (uint32_t)L * (uint32_t)L * (uint32_t)H; i++) {
		grid[i] = 0;
		ordered[i] = 0;
	}

	// Initialize with inputGrid if necessary
	if (inputGrid != nullptr) {
		maxh = sparseToDense(inputGrid, inputGridPoints, grid, L, H);
	}

	// Free *outputGrid memory if allocated
	// Done after initialization with inputGrid in case outGrid was used for it
	if (*outGrid != nullptr)
	{
		free(*outGrid);
		*outGrid = nullptr;
	}

	// Set up rotation values
	double dphi = turns * 360 / reps;
	if (dphi != 0.0) {
		if (stepper_resolution == 0) {
			printf("Rotating %e degrees between consecutive particles.\n", dphi);
		}
		else {
			printf("Rotating %e degrees between consecutive particles; clipped to increments of %e.\n", dphi, 360.0 / stepper_resolution);
		}
	}
	double theta_rad = theta * pi / 180.0;
	double phi_rad = phi * pi / 180.0;
	double dphi_rad = dphi * pi / 180.0;
	double phi_inc_rad = phi_rad;

	// Create surface
	Surface3D * surface = new Surface3D(L, H, L, grid, species, weights, seed);

	// Create random number generator
	if (seed == 0) {
		std::random_device rd;
		seed = rd();
	}
	std::mt19937 gen(seed);
	std::uniform_int_distribution<> dist(0, L - 1);

	std::mt19937 gens(seed);
	std::uniform_int_distribution<> sist(0, species->size());

	// Create path
	int32_t src[3] = { (int16_t)(round(maxh / tan(theta_rad - pi / 2.0) * cos(phi_rad))), (int16_t)(round(maxh / tan(theta_rad - pi / 2.0) * sin(phi_rad))), maxh };
	int32_t dest[3] = { 0,0,0 };
	uint16_t collision[3] = { 0, 0, 0 };
	uint16_t diffusing[3] = { 0, 0, 0 };
	uint16_t* vacancy;
	//uint16_t* points = (uint16_t*)malloc(sizeof(uint16_t) * (src[0] + src[1] + src[2] + 10) * 3);
	//uint16_t path_length = Bresenham3D(src, dest, points);

	// Print how long it took
	std::chrono::duration<double, std::milli> timep = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Set-up took " << timep.count() / 1000 << " s." << std::endl;

	start = std::chrono::high_resolution_clock::now();
	uint32_t deposited = 0;
	uint32_t not_deposited = 0;
	double timel = 0;
	double timed = 0;
	// Start deposition
	for (uint32_t n = 0; n < reps; n++) {
		auto startl = std::chrono::high_resolution_clock::now();

		// Generate random destination and calculated source point
		dest[0] = dist(gen);
		dest[1] = dist(gen);
		dest[2] = 0;

		// Calculate current phi and calculate source
		phi_inc_rad += dphi_rad;
		uint16_t v_offset = ceil(-L * 1.15 * tan(theta_rad - pi / 2.0));
		src[0] = (int16_t)(round((maxh + v_offset) / tan(theta_rad - pi / 2.0) * cos(phi_inc_rad))) + dest[0];
		src[1] = (int16_t)(round((maxh + v_offset) / tan(theta_rad - pi / 2.0) * sin(phi_inc_rad))) + dest[1];
		src[2] = (maxh + v_offset);

		// Send particle along
		traversePathRealTime(src, dest, collision, grid, L, H);

		std::chrono::duration<double, std::milli> dural = std::chrono::high_resolution_clock::now() - startl;
		timel += dural.count()/1000;

		// Diffuse particle
		auto startd = std::chrono::high_resolution_clock::now();
		diffusing[0] = collision[0];
		diffusing[1] = collision[1];
		diffusing[2] = collision[2];
		for (uint16_t d = 1; d < diffusion_steps + 1; d++) {
			vacancy = surface->getAdjacentVacancy(diffusing, (*species)[0]);
			diffusing[0] = vacancy[0];
			diffusing[1] = vacancy[1];
			diffusing[2] = vacancy[2];
		}
		std::chrono::duration<double, std::milli> durad = std::chrono::high_resolution_clock::now() - startd;
		timed += durad.count()/1000;

		// Add particle
		int8_t sp = (*species)[0];
		if (grid[diffusing[2] * L * L + diffusing[1] * L + diffusing[0]] == 0) {
			grid[diffusing[2] * L * L + diffusing[1] * L + diffusing[0]] = sp;
			ordered[diffusing[2] * L * L + diffusing[1] * L + diffusing[0]] = n + 1;
			surface->add(diffusing, sp);
			deposited += 1;
			if (diffusing[2] > maxh) {
				maxh = diffusing[2];
			}
		}
		else {
			not_deposited += 1;
			traversePathRealTime(src, dest, collision, grid, L, H);
		}

		// Periodic updates
		if (n % (reps/update) == (reps/update) - 1) {
			std::chrono::duration<double, std::milli> durup = std::chrono::high_resolution_clock::now() - start;
			std::cout << n + 1 << "/" << reps << " complete; " << durup.count() / 1000 << " seconds taken; " << (n + 1) / durup.count() * 1000 << " reps per second." << std::endl;
		}
	}

	

	// Finish up timing
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> dural = end - start;
	double time = dural.count()/1000;

	std::chrono::system_clock::time_point epoch = std::chrono::system_clock::now();
	uint32_t epoch_time = epoch.time_since_epoch().count() * std::chrono::system_clock::period::num / std::chrono::system_clock::period::den;


	// Update object for parameters
	SimulationParametersIndividual * layer_params = new SimulationParametersIndividual();
	layer_params->length = L;
	layer_params->width = L;
	layer_params->height = H;
	layer_params->theta = theta;
	layer_params->phi = phi;
	layer_params->turns = turns;
	layer_params->spread = spread;
	layer_params->stepper_resolution = stepper_resolution;
	layer_params->species = species;
	layer_params->diffusion_steps = diffusion_steps;
	layer_params->deposited = deposited;
	layer_params->repetitions = reps;
	layer_params->weights = weights;
	layer_params->system = system;
	layer_params->seed = seed;
	layer_params->time_taken = time;
	layer_params->time_finished = epoch_time;
	params->addLayer(layer_params);
	params->serialize();

	// Create filename
	std::string theta_str = std::to_string(theta);
	theta_str.erase(theta_str.find_last_not_of('0') + 1, std::string::npos);
	theta_str.erase(theta_str.find_last_not_of('.') + 1, std::string::npos);
	std::string turns_str = "";
		if (abs(params->turns) >= 1e-5) {
			turns_str += "_x" + std::to_string(params->turns);
			turns_str.erase(turns_str.find_last_not_of('0') + 1, std::string::npos);
			turns_str.erase(turns_str.find_last_not_of('.') + 1, std::string::npos);
		}
	std::string filename = "structures/STF_" + system + "_L" + std::to_string(L) + turns_str + "_Th" + theta_str + "_D" + std::to_string(diffusion_steps) + "_N" + std::to_string(params->deposited) + "_" + std::to_string(epoch_time);

	// Save objects
	uint32_t point_total = denseToSparse(grid, outGrid, L, H);
	cnpy::npz_save(filename + ".npz", "arr_0", *outGrid, { point_total, 4 });
	//free(sparse);
	uint32_t* sparse2 = nullptr;
	uint32_t point_total_ordered = denseToSparse(ordered, &sparse2, L, H);
	cnpy::npz_save(filename + "_ordered.npz", "arr_0", sparse2, { point_total_ordered, 4 });
	free(sparse2);

	std::ofstream json_file;
	json_file.open(filename + ".json");
	json_file << params->serialization;
	json_file.close();

	// Print the timing
	std::cout << reps << " reps completed in " << time << " seconds; \n" << reps / time << " reps per second." << std::endl;
	std::cout << "Line finding and traversal were " << timel << " seconds of that time." << std::endl;
	std::cout << "Diffusion was " << timed << " seconds of that time." << std::endl;

	// Free memory
	delete surface;
	free(grid);
	free(ordered);
	//free(points);
	return point_total;
}