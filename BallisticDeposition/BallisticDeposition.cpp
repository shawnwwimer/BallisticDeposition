// BallisticDeposition.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "Simulation3D.h"
#include "SimulationContinuous.h"

int main()
{
    bool cts_simulation = true;
    if (cts_simulation) {
        float theta = 85;
        float L = 128;
        float H = 128;
        uint32_t reps = 8192*16*8;
        uint8_t bin_size = 4;
        uint32_t seed = 0;
        float diffusion_length = 0.5;
        std::vector<int8_t> species = { 1 };
        std::vector<float> radii = { 0.15 };
        std::vector<std::vector<float>> weights = { {{1, .1}, {.1, 1}} };
        std::vector<std::vector<float>> inputGrid = { {0} };
        SimulationParametersFull params;
        std::string system = "Si";

        obliqueDepositionContinuous(theta, L, H, reps, bin_size, seed, diffusion_length, 10, &species, &radii, &weights, inputGrid, &params, system);
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
        float spreads[6] = { 1e-6, 1.f, 2.f, 3.f, 4.f, 5.f };

    }

    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
