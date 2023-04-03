#pragma once
#include <vector>
#include <array>

class CubicSpacePartition
{
private:
	float L;
	float H;
	float cube_size;
	int bins_on_side;
	int bins_num;
	std::vector<std::vector<int>> bins;

	
	
public:
	CubicSpacePartition(float L, float H, float cube_size = 2.0f) : L{ L }, H{ H }, cube_size{ cube_size }{
		bins_on_side = L / cube_size * 2;
		bins_num = bins_on_side * bins_on_side * H / cube_size * 2;
		for (int i = 0; i < bins_num; i++) {
			std::vector<int>* b = new std::vector<int>;
			bins.push_back(*b);
		}
	}

	void add_to_bins(int idx, std::array<float, 3> position) {
		float nx = position[0] / cube_size * 2;
		float ny = position[1] / cube_size * 2;
		float nz = position[2] / cube_size * 2;
		int xidx = round(nx);
		int yidx = round(ny);
		int zidx = round(nz);

		int xmod = xidx < nx ? 1 : -1;
		int ymod = yidx < ny ? 1 : -1;
		int zmod = zidx < nz ? 1 : -1;

		if (xmod + xidx == bins_on_side) {
			xmod = 0;
		}
		else if (xmod + xidx < 0) {
			xmod = bins_on_side - 1;
		}
		else {
			xmod += xidx;
		}
		if (ymod + yidx == bins_on_side) {
			ymod = 0;
		}
		else if (ymod + yidx < 0) {
			ymod = bins_on_side - 1;
		}
		else {
			ymod += yidx;
		}

		bins[zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx].push_back(idx); // add nearest
		bins[zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xmod].push_back(idx); // add nearest
		bins[zidx * bins_on_side * bins_on_side + ymod * bins_on_side + xidx].push_back(idx); // add nearest
		bins[zidx * bins_on_side * bins_on_side + ymod * bins_on_side + xmod].push_back(idx); // add nearest
		if (zmod > 0) {
			bins[zmod * bins_on_side * bins_on_side + yidx * bins_on_side + xidx].push_back(idx); // add nearest
			bins[zmod * bins_on_side * bins_on_side + yidx * bins_on_side + xmod].push_back(idx); // add nearest
			bins[zmod * bins_on_side * bins_on_side + ymod * bins_on_side + xidx].push_back(idx); // add nearest
			bins[zmod * bins_on_side * bins_on_side + ymod * bins_on_side + xmod].push_back(idx); // add nearest
		}
	}

	std::vector<int>* find_nearest_bin(std::array<float, 3> position) {
		float nx = position[0] / cube_size * 2;
		float ny = position[1] / cube_size * 2;
		float nz = position[2] / cube_size * 2;
		int xidx = round(nx);
		int yidx = round(ny);
		int zidx = round(nz);
		return &bins[zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx];
	}
};