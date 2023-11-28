#pragma once
#include <vector>
#include <array>

#include <nlopt.hpp>


class OverlappingCubicSpacePartition
{
private:
	float L;
	float H;
	float cube_size;
	int bins_on_side, bins_in_set, bins_num;
	std::vector<std::vector<int>> bins;
	std::vector<std::vector<std::vector<float>*>> neighbor_lists;
	std::vector<std::vector<float>>* atoms;
	nlopt::opt opt = nlopt::opt(nlopt::LN_COBYLA, 3);
	std::vector<double> minimization_result = { 0, 0, 0 };
	
public:
	OverlappingCubicSpacePartition(float L, float H, std::vector<std::vector<float>>* atoms, float cube_size = 2.0f) : L{ L }, H{ H }, cube_size{ cube_size }, atoms{ atoms }{
		bins_on_side = L / cube_size;
		bins_in_set = bins_on_side * bins_on_side * H / cube_size;
		bins_num = bins_in_set * 8;
		for (int i = 0; i < bins_num; i++) {
			std::vector<int>* b = new std::vector<int>;
			bins.push_back(*b);
			delete(b);
			std::vector<std::vector<float>*>* n = new std::vector<std::vector<float>*>;
			neighbor_lists.push_back(*n);
			delete(n);
		}

		// Set optimizer parameters
		opt.set_xtol_abs(1e-3);
		opt.set_xtol_rel(1e-3);
		opt.set_maxeval(400);
		if (opt.get_algorithm() == nlopt::LN_AUGLAG) {
			nlopt::opt local = nlopt::opt(nlopt::LN_NELDERMEAD, 3);
			local.set_xtol_abs(1e-3);
			local.set_xtol_rel(1e-3);
			local.set_maxeval(400);
			opt.set_local_optimizer(local);
		}
		else if (opt.get_algorithm() == nlopt::LD_AUGLAG) {
			nlopt::opt local = nlopt::opt(nlopt::LD_MMA, 3);
			local.set_xtol_abs(1e-4);
			local.set_xtol_rel(1e-4);
			local.set_maxeval(400);
			opt.set_local_optimizer(local);
		}
	}

	void add_to_bins(int idx, std::array<float, 3> * position) {
		float nx = (*position)[0] / cube_size;
		float ny = (*position)[1] / cube_size;
		float nz = (*position)[2] / cube_size;
		int xidx = floor(nx);
		int yidx = floor(ny);
		int zidx = floor(nz);
		bool offz = nz > 1 / cube_size;

		int xmod = xidx + (nx - xidx > 1 / cube_size ? 0 : -1);
		int ymod = yidx + (ny - yidx > 1 / cube_size ? 0 : -1);
		int zmod = zidx + (nz - zidx > 1 / cube_size ? 0 : -1);

		if (xmod == bins_on_side) {
			xmod = 0;
		}
		else if (xmod < 0) {
			xmod = bins_on_side - 1;
		}

		if (ymod == bins_on_side) {
			ymod = 0;
		}
		else if (ymod < 0) {
			ymod = bins_on_side - 1;
		}

		if (zmod < 0) {
			zmod = 0;
		}

		bins[zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx].push_back(idx); // axis-aligned
		bins[(zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xmod) + bins_in_set].push_back(idx); // x-offset
		bins[(zidx * bins_on_side * bins_on_side + ymod * bins_on_side + xidx) + bins_in_set * 2].push_back(idx); // y-offset
		bins[(zidx * bins_on_side * bins_on_side + ymod * bins_on_side + xmod) + bins_in_set * 3].push_back(idx); // xy-offset
		if (offz) {
			bins[(zmod * bins_on_side * bins_on_side + yidx * bins_on_side + xidx) + bins_in_set * 4].push_back(idx); // z-offset
			bins[(zmod * bins_on_side * bins_on_side + yidx * bins_on_side + xmod) + bins_in_set * 5].push_back(idx); // xz-offset
			bins[(zmod * bins_on_side * bins_on_side + ymod * bins_on_side + xidx) + bins_in_set * 6].push_back(idx); // yz-offset
			bins[(zmod * bins_on_side * bins_on_side + ymod * bins_on_side + xmod) + bins_in_set * 7].push_back(idx); // xyz-offset
		}

		//neighbor_lists[zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx].push_back(&((*atoms)[idx])); // axis-aligned
		//neighbor_lists[(zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xmod) + bins_in_set].push_back(&((*atoms)[idx])); // x-offset
		//neighbor_lists[(zidx * bins_on_side * bins_on_side + ymod * bins_on_side + xidx) + bins_in_set * 2].push_back(&((*atoms)[idx])); // y-offset
		//neighbor_lists[(zidx * bins_on_side * bins_on_side + ymod * bins_on_side + xmod) + bins_in_set * 3].push_back(&((*atoms)[idx])); // xy-offset
		//if (zmod > 0) {
		//	neighbor_lists[(zmod * bins_on_side * bins_on_side + yidx * bins_on_side + xidx) + bins_in_set * 4].push_back(&((*atoms)[idx])); // z-offset
		//	neighbor_lists[(zmod * bins_on_side * bins_on_side + yidx * bins_on_side + xmod) + bins_in_set * 5].push_back(&((*atoms)[idx])); // xz-offset
		//	neighbor_lists[(zmod * bins_on_side * bins_on_side + ymod * bins_on_side + xidx) + bins_in_set * 6].push_back(&((*atoms)[idx])); // yz-offset
		//	neighbor_lists[(zmod * bins_on_side * bins_on_side + ymod * bins_on_side + xmod) + bins_in_set * 7].push_back(&((*atoms)[idx])); // xyz-offset
		//}
	}

	int find_nearest_bin_idx(std::array<float, 3> * position) {
		float nx = round((*position)[0] / cube_size * 2)/2;
		float ny = round((*position)[1] / cube_size * 2)/2;
		float nz = round((*position)[2] / cube_size * 2)/2;
		
		// if nx is odd we're offset in x
		// if ny is odd we're offset in y
		// if nz is odd we're offset in z
		int factor = 0;
		if (int(nx*2) % 2 == 1) {
			if (int(ny*2) % 2 == 1) {
				if (fmod(nz, 1) == 0) {
					factor = 4; // xyz-offset
				}
				else {
					factor = 0; // xy-offset
				}
			}
			else {
				if (fmod(nz, 1) == 0) {
					factor = 6; // xz-offset
				}
				else {
					factor = 2; // x-offset
				}
			}
		}
		else {
			if (int(ny*2) % 2 == 1) {
				if (fmod(nz, 1) == 0) {
					factor = 5; // yz-offset
				}
				else {
					factor = 1; // y-offset
				}
			}
			else {
				if (fmod(nz, 1) == 0) {
					factor = 7; // z-offset
				}
				else {
					factor = 3;  // axis-aligned
				}
			}
		}

		if (nz == 0) {
			factor -= 4;
		}

		int xidx = int(ceil(nx - 1.f)) % bins_on_side;
		int yidx = int(ceil(ny - 1.f)) % bins_on_side;
		int zidx = int(ceil(nz - 1.f));
		if (xidx < 0) {
			xidx = bins_on_side - 1;
		}
		if (yidx < 0) {
			yidx = bins_on_side - 1;
		}
		if (zidx < 0) {
			zidx = 0;
		}

		return zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx + bins_in_set * factor;
	}

	std::vector<int>* find_nearest_bin(std::array<float, 3>* position) {
		return &bins[find_nearest_bin_idx(position)];
	}

	std::vector<double>* find_local_minimum(std::array<float, 3>* position, float distance);
};