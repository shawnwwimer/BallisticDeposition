#include "math_utils.h"



uint32_t flat_index(uint16_t* point, uint16_t L, uint16_t W) {
	return (point[2] * L * L + point[1] * L + point[0]);
}

uint32_t flat_index_PBC(uint16_t* point, uint16_t L, uint16_t W) {
	return (point[2] * L * L + modulo(point[1], L) * L + modulo(point[0], L));
}

void wrap_index(uint32_t idx, uint16_t* idcs, uint16_t L, uint16_t H) {
	idcs[0] = idx % L;
	idcs[1] = (idx / L) % L;
	idcs[2] = (idx / L / L) % L;
}