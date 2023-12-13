#pragma once
#include <vector>
#include <array>
#define _USE_MATH_DEFINES
#include <cmath>


class SpaceHashMap
{
private:
	float L;
	float H;
	float cube_size;
	
	//std::vector<std::vector<float>>* atoms;

public:
	SpaceHashMap(float L, float H, std::vector<std::vector<float>>* atoms, float cube_size = 0.5f) : L{ L }, H{ H }, cube_size{ cube_size }/*, atoms{atoms}*/{
		bins_on_side = L / cube_size;
		bins_in_layer = bins_on_side * bins_on_side;
		bins_num = bins_in_layer * H / cube_size;
		yswitch = bins_on_side * (bins_on_side - 1);
		for (int i = 0; i < bins_num; i++) {
			std::vector<int>* b = new std::vector<int>;
			bins.push_back(*b);
			delete(b);
		}
	}

	void add_to_bin(std::array<float, 3>* position, int idx) {
		int xidx = (int) ((*position)[0] / cube_size);
		int yidx = (int) ((*position)[1] / cube_size);
		int zidx = (int) ((*position)[2] / cube_size);

		bins[zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx].push_back(idx);
	}

	std::vector<int>* return_bin(std::array<float, 3>* position) {
		int xidx = (int) ((*position)[0] / cube_size);
		int yidx = (int) ((*position)[1] / cube_size);
		int zidx = (int) ((*position)[2] / cube_size);
		return &bins[zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx];
	}

	int return_bin_idx(std::array<float, 3>* position) {
		int xidx = (int) ((*position)[0] / cube_size);
		int yidx = (int) ((*position)[1] / cube_size);
		int zidx = (int) ((*position)[2] / cube_size);
		return zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx;
	}

	int return_bin_idx(float x, float y, float z) {
		int xidx = (int) (x / cube_size);
		int yidx = (int) (y / cube_size);
		int zidx = (int) (z / cube_size);
		return zidx * bins_on_side * bins_on_side + yidx * bins_on_side + xidx;
	}

	int bins_on_side, bins_in_layer, bins_num;
	int yswitch;
	std::vector<std::vector<int>> bins;
};