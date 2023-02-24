// BallisticDeposition.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "Simulation3D.h"
#include "SimulationContinuous.h"

int main()
{
    bool discrete = false;

    if (discrete) {
        float theta = 85;
        uint16_t L = 768;
        uint16_t H = 384;
        uint32_t reps = 8192 * 32 * 8 * 4;// * 2;
        float phi = 0;
        float turns = 0;
        uint32_t seed = 3286398647;
        uint16_t diffusion_steps = 5;
        std::vector<int8_t> species = { 1 };
        std::vector<float> spread = { 0.000001, 0.000001 };
        std::vector<std::vector<float>> weights = { { {1, .1}, { .1, 1 }} };
        uint32_t stepper_resolution = 0;
        SimulationParametersFull params;
        int16_t** outGrid = (int16_t**)malloc(sizeof(int16_t*));
        *outGrid = nullptr;
        std::string system = "Si";
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
        for (int n = 0; n < 1; n++) {
            std::cout << "Starting simulation " << n + 1 << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, diffusion_steps, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, stepper_resolution, &params, system, true, false);
            params.clearLayers();
            phi += 15;
        }
    }
    else {
        float theta = 85;
        float L = 32;
        float H = 32;
        uint32_t reps = 8192;
        uint8_t bin_size = 8;
        uint32_t seed = 3286398647;
        float diffusion_length = 0;
        std::vector<int8_t> species = { 1 };
        std::vector<float> radii = { 0.15 };
        std::vector<std::vector<float>> weights = { { {1, .1}, { .1, 1 }} };
        std::vector<std::vector<float>> inputGrid = { {0} };
        SimulationParametersFull params;
        std::string system = "Si";

        obliqueDepositionContinuous(theta, L, H, reps, bin_size, seed, diffusion_length, &species, &radii, &weights, inputGrid, &params, system);
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
