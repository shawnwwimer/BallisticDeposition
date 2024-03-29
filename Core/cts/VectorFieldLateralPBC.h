#pragma once

#include <stdlib.h>
#include <array>
#include "../math_utils.h"

class VectorFieldLateralPBC
{
private:
	std::array<float, 3>** arr;
	float L, W, H;
	float scale;
	int Ls, Ws, Hs, WLs;

public:
	/// <summary>
	/// Create vector field that obeys PBC in the lateral dimensions.
	/// </summary>
	/// <param name="L">x-period</param>
	/// <param name="W">y-period</param>
	/// <param name="H">z-extent</param>
	/// <param name="scale">Some number such that L/scale, W/scale, and H/scale are all integers larger than the dimension. Default: 0.1f. </param>
	VectorFieldLateralPBC(float L, float W, float H, float scale = 10.f) : L{ L }, W{ W }, H{ H }, scale{ scale }{
		Ls = L * scale;
		Ws = W * scale;
		Hs = H * scale;
		WLs = Ls * Ws;
		arr = (std::array<float, 3>**)malloc(sizeof(std::array<float, 3>*) * Ls * Ws * Hs);
		for (int i = 0; i < Ls * Ws * Hs; i++) {
			arr[i] = new std::array<float, 3>();
			(*arr[i])[0] = 0;
			(*arr[i])[1] = 0;
			(*arr[i])[2] = 0;
		}
	}

	~VectorFieldLateralPBC() {
		for (int i = 0; i < Ls * Ws * Hs; i++) {
			delete arr[i];
		}

		free(arr);
	}

	/// <summary>
	/// Get vector at specific point in matrix.
	/// </summary>
	/// <param name="x">x-coordinate</param>
	/// <param name="y">y-coordinate</param>
	/// <param name="z">z-coordinate</param>
	/// <returns>nullptr if z out of range, std::array&lt;float, 3&gt;* if not.</returns>
	std::array<float, 3>* operator()(int x, int y, int z) {
		if (z > H || z < 0) {
			return nullptr;
		}

		if (x > L || x < 0) {
			x = modulof(x, L);
		}

		if (y > W || y < 0) {
			y = modulof(y, W);
		}

		return arr[x + y * Ls + z * WLs];
	}

	int get_Ls() {
		return Ls;
	}

	int get_Ws() {
		return Ws;
	}

	int get_Hs() {
		return Hs;
	}


	/*int find_local_minimum(std::array<float, 3>& xyz, float radius) {
		if (xyz[2] > H || xyz[2] < 0) {
			return 1;
		}

		std::array<int, 3> minimum = { xyz[0] * scale, xyz[1] * scale, xyz[2] * scale };

		int x = minimum[0];
		int y = minimum[0];
		int z = minimum[0];

		float min = arr(x, y, z);

		int rad_disc = radius * scale;

		for (int k = -rad_disc; k < rad_disc; k++) {
			int kk =
				for (int j = -rad_disc; j < rad_disc; j++) {
					for (int i = -rad_disc; i < rad_disc; i++) {

					}
				}

		}

		if (x > L || x < 0) {
			x = modulof(x, L);
		}

		if (y > W || y < 0) {
			y = modulof(y, W);
		}
	}*/
};
