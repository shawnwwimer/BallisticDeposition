#pragma once

#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_
#include<cstdint>

template<typename T>
inline T modulo(T val, uint16_t mod) {
	return ((val % mod + mod) % mod);
}

uint32_t flat_index(uint16_t* point, uint16_t L, uint16_t W);

uint32_t flat_index_PBC(uint16_t* point, uint16_t L, uint16_t W);

void wrap_index(uint32_t idx, uint16_t* idcs, uint16_t L, uint16_t H);



#endif