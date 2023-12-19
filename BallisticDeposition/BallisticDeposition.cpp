// BallisticDeposition.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "Simulation3D.h"
#include "SimulationContinuous.h"
#include "SimulationContinuousHashedSpace.h"

// should be modified to take radius into account
float estimate_binsize(float flux_angle, float radius) {
    // First find the minimum BS
    float BS = 2 * radius / cos(flux_angle * M_PI / 180);

    // Now get the next-largest power of 2
    float f = floor(log2(BS) + 1);

    // Return that power of 2
    return pow(2, f);
}

int main()
{
    bool cts_simulation = true;
    if (cts_simulation) {
        float theta = 85;
        float L = 64;
        float H = 32;

        uint32_t reps = 128*128*8*2;
        uint8_t bin_size = 2;
        uint32_t seed = 0;//1277363101;
        float diffusion_length = 0;
        std::vector<int8_t> species = { 1 };
        std::vector<float> radii = { 0.111 }; // Si: 0.111; Ag: 0.144; Ti: 0.147
        std::vector<float> spread = { 5, 5 };
        std::vector<std::vector<float>> weights = { {{1, .1}, {.1, 1}} };
        std::vector<std::vector<float>> inputGrid;
        ContinuousSimulationParametersFull params;
        std::string system = "Si";

        FilesToSave save_params;
        save_params.priority = false;
        //save_params.collisions = false;
        //inputGrid = { {51.87532, 1, 3.860409, 1, 0.147, 1} };
        int treps = reps;
        std::vector<float> thetas = { 85 };
        std::vector<float> diffs = { 0.2 };
        for (float d : diffs) {
            for (float t : thetas) {
                obliqueDepositionContinuousHashed(t, L, H, treps, estimate_binsize(t, radii[0]), seed, d, 11, &species, &radii, &spread, &weights, inputGrid, &params, system, DiffusionMethod::HopAndSettleLUT, &save_params);

                params.clearLayers();
            }
        }
    }
    else {
        float theta = 85;
        uint16_t L = 768;
        uint16_t H = 800;
        uint32_t reps = 8192 * 32 * 8 * 4 * 2;
        float phi = 0;
        float turns = 2;
        uint32_t seed = 0;
        uint16_t diffusion_steps = 5;
        std::vector<int8_t> sSi = { 1 };
        std::vector<int8_t> sAg = { 2 };
        std::vector<int8_t> sZrO2 = { 3, 4, 4 };
        std::vector<float> spread = { 1e-6, 1e-6 };
        std::vector<float> spread1 = { 2, 2 };
        std::vector<float> spread2 = { 4, 4 };
        std::vector<std::vector<float>> Si = { { {1, 0}, { 0, 1 }} };
        std::vector<std::vector<float>> SiAg = { {1, 0, 0}, {0, 1, 0}, {0.2, -.3, 3} };
        std::vector<std::vector<float>> ZrO2 = { {1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0.2, -0.3, 3, -.3, -.3}, {0, 0, 0, 1, 5}, {0, 0, 0, 5, -0.3} }; // 0/Si/Ag/Zr/O

        uint32_t stepper_resolution = 0;
        SimulationParametersFull params;
        int16_t** outGrid = (int16_t**)malloc(sizeof(int16_t*));
        *outGrid = nullptr;
        std::string ySi = "Si";
        std::string ySiAg = "Si_Ag";
        std::string yZrO2 = "ZrO2";
        int16_t* inputGrid = (int16_t*)malloc(sizeof(int16_t) * L * L * 4);
        uint32_t count = 0;
        for (uint32_t i = 0; i < L; i++) {
            for (uint32_t j = 0; j < L; j++) {
                inputGrid[count * 4] = i;
                inputGrid[count * 4 + 1] = j;
                inputGrid[count * 4 + 2] = 0;
                inputGrid[count * 4 + 3] = 1;
                count += 1;
            }
        }

        uint32_t inputGridPoints = count;
        uint32_t points = 0;
        std::vector<float> thetas = { 80, 82, 84, 85, 86, 88, 89 };
        std::vector<int> Ds = { 0, 1, 2, 5, 10, 15, 20, 25, 30, 50, 75, 100 };
        float spreads[6] = { 1.f, 2.f, 3.f, 4.f, 5.f };

        
        for (int n = 0; n < 5; n++) {
            for (int d = 0; d < Ds.size(); d++) {
                for (int t = 0; t < thetas.size(); t++) {
                    for (int s = 0; s < 1; s++) {
                        spread[0] = spreads[s];
                        spread[1] = spreads[s];
                        std::cout << "Simulating theta " << thetas[t] << " degrees and " << Ds[d] << " diffusion steps." << std::endl;
                        points = obliqueDeposition(thetas[t], L, H, reps, 0, 1, seed, Ds[d], &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
                        params.clearLayers();
                        //points = obliqueDeposition(thetas[t], L, H, reps, 0, , seed, Ds[i], &sZrO2, &spread, &ZrO2, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, yZrO2, false, false, Acceleration::NONE);
                        //params.clearLayers();
                    }
                }
            }
        }

        
    }

    return 0;
}