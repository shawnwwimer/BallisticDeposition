#pragma once
#include <stdlib.h>
#include <array>
#include <vector>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <random>
#include "cnpy.h"
#include "../math_utils.h"


class Matrix3DLateralPBC
{
private:
	float * arr = nullptr;
	float L, W, H;
	float scale;
	size_t Ls, Ws, Hs, WLs;
	bool init = false;
	std::vector<std::array<float, 3>> points;
	std::array<std::array<float, 3>, 128> minima;

public:
	/// <summary>
	/// Create matrix that obeys PBC in the lateral dimensions.
	/// </summary>
	/// <param name="L">x-period</param>
	/// <param name="W">y-period</param>
	/// <param name="H">z-extent</param>
	/// <param name="scale">Some number such that L/scale, W/scale, and H/scale are all integers larger than the dimension. Default: 0.2f. </param>
	Matrix3DLateralPBC(float L, float W, float H, float scale = 5.f) : L{ L }, W{ W }, H{ H }, scale{ scale }{
		Ls = L * (double)scale + 0.5;
		Ws = W * (double)scale + 0.5;
		Hs = H * (double)scale + 0.5;
		WLs = Ls * Ws;
	}
	
	void initialize() {
		if (!init) {
			arr = (float*)malloc(sizeof(float) * Ls * Ws * Hs);
			if (arr != nullptr) {
				for (int i = 0; i < Ls * Ws * Hs; i++) {
					arr[i] = 0;
				}
				init = true;
			}
		}
	}

	~Matrix3DLateralPBC() {
		free(arr);
	}

	/// <summary>
	/// Get value at specific point in matrix.
	/// </summary>
	/// <param name="x">x-coordinate</param>
	/// <param name="y">y-coordinate</param>
	/// <param name="z">z-coordinate</param>
	/// <returns>nullptr if z out of range, std::array&lt;float, 3&gt;* if not.</returns>
	inline float& operator()(int x, int y, int z) {
		if (z > Hs || z < 0) {
			return arr[0, 0, 0];
		}

		if (x > Ls || x < 0) {
			x = modulo(x, Ls);
		}

		if (y > Ws || y < 0) {
			y = modulo(y, Ws);
		}

		return arr[x + y * Ls + z * WLs];
	}

	/// <summary>
	/// Get value at specific point in matrix.
	/// </summary>
	/// <param name="xyz">Coordinate array</param>
	/// <returns>nullptr if z out of range, std::array&lt;float, 3&gt;* if not.</returns>
	inline float& operator()(std::array<int, 3> xyz) {
		int x = xyz[0];
		int y = xyz[1];
		int z = xyz[2];
		if (z > Hs || z < 0) {
			return arr[0, 0, 0];
		}

		if (x > Ls || x < 0) {
			x = modulo(x, Ls);
		}

		if (y > Ws || y < 0) {
			y = modulo(y, Ws);
		}

		return arr[x + y * Ls + z * WLs];
	}

	inline void real_point(std::array<int, 3>& xyz, std::array<float, 3>& output) {
		output = { float(xyz[0]) / scale, float(xyz[1]) / scale, float(xyz[2]) / scale };
	}

	inline void scaled_point(std::array<float, 3>& xyz, std::array<int, 3>& output) {
		output = { int(xyz[0] * scale), int(xyz[1] * scale), int(xyz[2] * scale) };
	}

	inline int get_Ls() {
		return Ls;
	}

	inline int get_Ws() {
		return Ws;
	}

	inline int get_Hs() {
		return Hs;
	}

	/// <summary>
	/// Find the local minimum.
	/// </summary>
	/// <param name="xyz">Three-element float array of location in nm to be updated.</param>
	/// <param name="radius">Radius in nm over which to search for minimum.</param>
	void find_local_minimum(std::array<float, 3>& xyz, float radius, float choice=1.f) {
		if (xyz[2] > H || xyz[2] < 0 || isnan(xyz[2])) {
			throw std::out_of_range("Z coordinate out of range.");
		}

		int x = round(xyz[0] * scale);
		int y = round(xyz[1] * scale);
		int z = round(xyz[2] * scale);

		/*if (x > L || x < 0) {
			x = modulof(x, L);
		}

		if (y > W || y < 0) {
			y = modulof(y, W);
		}*/

		// initialize minimum
		std::array<int, 3> minimum = { x, y, z };
		float min_val = arr[x + y * Ls + z * WLs];
		float min_distance2 = 0;
		minima[0] = xyz;
		int num_minima = 1;

		// only search in the nearby cube of discretized radius
		int rad_disc = radius * scale;
		for (int k = 0; k < rad_disc; k++) {
			// create and ignore out-of-range z-index, then update flat
			int kk = z + k;
			if (kk < 0 || kk > Hs - 1) {
				continue;
			}
			int idxz = kk * WLs;
			for (int j = rad_disc; j > -rad_disc-1; j--) {
				// create and range y-index, then update flat
				int jj = y + j;
				if (jj > Ws || jj < 0) {
					jj = modulof(jj, Ws);
				}
				int idxy = idxz + jj * Ls;
				for (int i = rad_disc; i > -rad_disc-1; i--) {
					// create and range x-index, then update flat
					int ii = x + i;
					if (ii > Ls || ii < 0) {
						ii = modulof(ii, Ls);
					}
					int idxx = idxy + ii;

					// if it's the minimum update
 					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}
				}
			}

			// create and ignore out-of-range z-index, then update flat
			kk = z - k;
			if (kk < 0 || kk > Hs - 1) {
				continue;
			}
			idxz = kk * WLs;
			for (int j = -rad_disc; j < rad_disc + 1; j++) {
				// create and range y-index, then update flat
				int jj = y + j;
				if (jj > Ws || jj < 0) {
					jj = modulof(jj, Ws);
				}
				int idxy = idxz + jj * Ls;
				for (int i = rad_disc; i > -rad_disc - 1; i--) {
					// create and range x-index, then update flat
					int ii = x + i;
					if (ii > Ls || ii < 0) {
						ii = modulof(ii, Ls);
					}
					int idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}
				}
			}
		}

		uint32_t idx = (uint32_t)(floorf(choice * (num_minima)));

		/*for (int k = 0; k < rad_disc; k++) {
			// create and ignore out-of-range +z-index, then update flat
			int kk = z + k;
			if (kk > Hs - 1) {
				continue;
			}
			int idxz = kk * WLs;
			for (int j = 0; j < rad_disc; j++) {
				// create and range +y-index, then update flat
				int jj = y + j;
				if (jj > Ws - 1) {
					jj = modulof(jj, Ws);
				}
				int idxy = idxz + jj * Ls;
				for (int i = 0; i < rad_disc; i++) {
					// create and range +x-index, then update flat
					int ii = x + i;
					if (ii > Ls - 1) {
						ii = modulof(ii, Ls);
					}
					int idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}

					// create and range -x-index, then update flat
					ii = x - i;
					if (ii < 0) {
						ii = modulof(ii, Ls);
					}
					idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}
				}

				// create and range -y-index, then update flat
				jj = y - j;
				if (jj < 0) {
					jj = modulof(jj, Ws);
				}
				idxy = idxz + jj * Ls;
				for (int i = 0; i < rad_disc; i++) {
					// create and range +x-index, then update flat
					int ii = x + i;
					if (ii > Ls || ii < 0) {
						ii = modulof(ii, Ls);
					}
					int idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}

					// create and range -x-index, then update flat
					ii = x - i;
					if (ii > Ls || ii < 0) {
						ii = modulof(ii, Ls);
					}
					idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}
				}
			}

			// create and ignore out-of-range -z-index, then update flat
			kk = z - k;
			if (kk < 0) {
				continue;
			}
			idxz = kk * WLs;
			for (int j = 0; j < rad_disc; j++) {
				// create and range +y-index, then update flat
				int jj = y + j;
				if (jj > Ws - 1) {
					jj = modulof(jj, Ws);
				}
				int idxy = idxz + jj * Ls;
				for (int i = 0; i < rad_disc; i++) {
					// create and range +x-index, then update flat
					int ii = x + i;
					if (ii > Ls - 1) {
						ii = modulof(ii, Ls);
					}
					int idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}

					// create and range -x-index, then update flat
					ii = x - i;
					if (ii < 0) {
						ii = modulof(ii, Ls);
					}
					idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}
				}

				// create and range -y-index, then update flat
				jj = y - j;
				if (jj < 0) {
					jj = modulof(jj, Ws);
				}
				idxy = idxz + jj * Ls;
				for (int i = 0; i < rad_disc; i++) {
					// create and range +x-index, then update flat
					int ii = x + i;
					if (ii > Ls - 1) {
						ii = modulof(ii, Ls);
					}
					int idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}

					// create and range -x-index, then update flat
					ii = x - i;
					if (ii < 0) {
						ii = modulof(ii, Ls);
					}
					idxx = idxy + ii;

					// if it's the minimum update
					if (arr[idxx] < min_val) {
						min_distance2 = k * k + j * j + i * i;
						min_val = arr[idxx];
						minimum = { ii, jj, kk };
						//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						num_minima = 1;
					}
					else if (arr[idxx] == min_val) {
						float distance2 = k * k + j * j + i * i;
						if (distance2 < min_distance2) {
							min_distance2 = distance2;
							min_val = arr[idxx];
							minimum = { ii, jj, kk };
							//xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							minima[0] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima = 1;
						}
						else if (distance2 == min_distance2) {
							minima[num_minima] = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
							num_minima += 1;
						}
					}
				}
			}
		}*/

		/*for (int size = 1; size < rad_disc + 1; size++) {
			for (int i = -x - size + 1; i < x + size; i++) {
				int ii = i;
				if (ii > Ls || ii < 0) {
					ii = modulof(ii, Ls);
				}
				int idxx = ii;

				for (int j = -y - size + 1; j < y + size; j++) {
					int jj = j;
					if (jj > Ws || jj < 0) {
						jj = modulof(jj, Ws);
					}
					int idxy = idxx + jj * Ls;

					for (int k = -z - size + 1; k < z + size; k++) {
						if (k < 0 || k > Hs - 1) {
							continue;
						}

						int idx = idxy + k * WLs;
						// if it's the minimum update
						if (arr[idx] < min_val) {
							min_val = arr[idxx];
							minimum = { ii, jj, k };
							xyz = { minimum[0] / scale, minimum[1] / scale, minimum[2] / scale };
						}
					}
				}
			}
		}*/


		points.push_back(xyz);
		xyz = minima[idx];
	}

	bool save_file(const char* fname) {
		if (init) {
			/*float* out = (float*)malloc(sizeof(float) * Ls * Ws * Hs);
			for (int i = 0; i < Ls; i++) {
				for (int j = 0; j < Ws; j++) {
					for (int k = 0; k < Hs; k++) {
						out[k * Ls * Ws + j * Ls + i] = points[;
					}
				}
			}*/
			cnpy::npy_save(fname, arr, { Ls, Ws, Hs });
			//free(out);
			return true;
		}
		return false;

	}
};

