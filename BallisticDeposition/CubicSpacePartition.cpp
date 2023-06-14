#include "CubicSpacePartition.h"

double local_potential(unsigned n, const double* x, double* grad, void* my_func_data) {
	std::vector<std::vector<float>>* positions = (std::vector<std::vector<float>>*) my_func_data;

	float val = 0;
	float r2, r6 = 0;
	float s6, s12 = 0;
	float s6r6, s12r12 = 0;
	for (int i = 0; i < positions->size(); i++) {
		auto p = (*positions)[i];
		r2 = pow(p[0] - x[0], 2) + pow(p[1] - x[1], 2) + pow(p[2] - x[2], 2);
		r6 = pow(r2, 3);
		s6 = pow(2 * p[4], 6) / 2;
		s12 = pow(s6, 2);
		s6r6 = s6 / r2;
		s12r12 = pow(s6r6, 2);
		val += s12r12 - s6r6;
		if (grad) {
			grad[0] += (6 * s6r6 / r2 - 12 * s12r12 / r2) * (x[0] - p[0]);
			grad[1] += (6 * s6r6 / r2 - 12 * s12r12 / r2) * (x[1] - p[1]);
			grad[2] += (6 * s6r6 / r2 - 12 * s12r12 / r2) * (x[2] - p[2]);
		}
	}
	return val;
}

std::vector<double>* CubicSpacePartition::find_local_minimum(std::array<float, 3> position, float distance) {
	// Get the bin
	std::vector<int>* bin = find_nearest_bin(position);

	// Build list of neighbor positions
	std::vector<std::vector<float>> neighbors;
	for (auto& p : *bin) {
		neighbors.push_back((*atoms)[p]);
	}

	// Create function
	opt.set_min_objective(local_potential, &neighbors);

	// Set bounds
	double x_min = position[0] - distance;
	double y_min = position[1] - distance;
	double z_min = position[2] - distance > 0 ? position[2] - distance : 0;
	double x_max = position[0] + distance;
	double y_max = position[1] + distance;
	double z_max = position[2] + distance < H ? position[2] + distance : H;
	opt.set_lower_bounds({ x_min, y_min, z_min });
	opt.set_upper_bounds({ x_max, y_max, z_max });

	// Minimize
	double residual;
	minimization_result = { position[0], position[1], position[2] };
	nlopt::result result = opt.optimize(minimization_result, residual);

	return &minimization_result;
}