#include "CubicSpacePartition.h"

double local_potential(unsigned n, const double* x, double* grad, void* my_func_data) {
	std::vector<std::vector<float>*>* positions = (std::vector<std::vector<float>*>*) my_func_data;

	float val = 0;
	float r2, r6 = 0;
	float s6, s12 = 0;
	float s6r6, s12r12 = 0;
	float dx, dy = 0;

	if (grad) {
		for (int i = 0; i < 3; i++) {
			grad[i] = 0;
		}
	}
	
	
	for (int i = 0; i < positions->size(); i++) {
		
		auto p = (*positions)[i];
		dx = (*p)[0] - x[0];
		dy = (*p)[1] - x[1];
		if (dx > 96 / 2) {
			dx -= 96;
		}
		else if (dx < -96 / 2) {
			dx += 96;
		}
		if (dy > 96 / 2) {
			dy -= 96;
		}
		else if (dy < -96 / 2) {
			dy += 96;
		}
		r2 = pow(dx, 2) + pow(dy, 2) + pow((*p)[2] - x[2], 2);
		r6 = pow(r2, 3);
		s6 = pow(2 * (*p)[4], 6) / 2;
		s12 = pow(s6, 2);
		s6r6 = s6 / r6;
		s12r12 = pow(s6r6, 2);
		val += s12r12 - s6r6;
		if (grad) {
			grad[0] += (6 * s6 * r6 - 12 * s12) / r6 / r6 / r2 * dx;
			grad[1] += (6 * s6 * r6 - 12 * s12) / r6 / r6 / r2 * dy;
			grad[2] += (6 * s6 * r6 - 12 * s12) / r6 / r6 / r2 * (x[2] - (*p)[2]);
		}
	}
	return val;
}

std::vector<double>* CubicSpacePartition::find_local_minimum(std::array<float, 3>* position, float distance) {
	// Get the bin
	int binidx = find_nearest_bin_idx(position);
	std::vector<int>* bin = &bins[binidx];

	// Build list of neighbor positions
	std::vector<std::vector<float>*> neighbors;
	for (auto& p : *bin) {
		neighbors.push_back(&(*atoms)[p]);
		//std::array<float, 3> npos = { (*atoms)[p][0], (*atoms)[p][1], (*atoms)[p][2] };
		//int nidx = find_nearest_bin_idx(npos);
		//if (nidx != binidx) {
		//	nidx = nidx;
		//}
	}

	// Get the neighbor list
	//std::vector<std::vector<float>*> neighbors = neighbor_lists[find_nearest_bin_idx(position)];

	// Create function
	opt.set_min_objective(local_potential, &neighbors);

	// Set bounds
	double x_min = (*position)[0] - distance > 0 ? (*position)[0] - distance : 0;
	double y_min = (*position)[1] - distance > 0 ? (*position)[1] - distance : 0;
	double z_min = (*position)[2] - distance > 0 ? (*position)[2] - distance : 0;
	double x_max = (*position)[0] + distance < L ? (*position)[0] + distance : L;
	double y_max = (*position)[1] + distance < L ? (*position)[1] + distance : L;
	double z_max = (*position)[2] + distance < H ? (*position)[2] + distance : H;
	opt.set_lower_bounds({ x_min, y_min, z_min });
	opt.set_upper_bounds({ x_max, y_max, z_max });

	// Minimize
	double residual;
	minimization_result = { (*position)[0], (*position)[1], (*position)[2] };
	nlopt::result result = opt.optimize(minimization_result, residual);

	return &minimization_result;
}