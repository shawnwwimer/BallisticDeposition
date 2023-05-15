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
        float L = 88;
        float H = 24;
        uint32_t reps = 65536*4*4;
        uint8_t bin_size = 2;
        uint32_t seed = 0;
        float diffusion_length = 1;
        std::vector<int8_t> species = { 1 };
        std::vector<float> radii = { 0.147 }; // Si: 0.111; Ag: 0.144; Ti: 0.147
        std::vector<std::vector<float>> weights = { {{1, .1}, {.1, 1}} };
        std::vector<std::vector<float>> inputGrid = { {0} };
        ContinuousSimulationParametersFull params;
        std::string system = "Si";

        std::vector<float> thetas = { 85 };
        for (float t : thetas) {
            std::cout << "Deposition at " << t << std::endl;
            obliqueDepositionContinuous(t, L, H, reps, bin_size, seed, diffusion_length, 25, &species, &radii, &weights, inputGrid, &params, system);
            params.clearLayers();
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
        float spreads[6] = { 1e-6, 1.f, 2.f, 3.f, 4.f, 5.f };

        /*points = obliqueDeposition(85.0476, 768, 800, reps * 40.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        points = obliqueDeposition(83.0166, 768, 800, reps * 130.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        points = obliqueDeposition(87.2449, 768, 800, reps * 13.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        params.clearLayers();
        reps *= 2;
        points = obliqueDeposition(85.0476, 768, 800, reps * 40.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        points = obliqueDeposition(83.0166, 768, 800, reps * 130.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        points = obliqueDeposition(87.2449, 768, 800, reps * 13.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        params.clearLayers();
        reps *= 2;
        points = obliqueDeposition(85.0476, 768, 800, reps * 40.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        points = obliqueDeposition(83.0166, 768, 800, reps * 130.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        points = obliqueDeposition(87.2449, 768, 800, reps * 13.0 / 183, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        params.clearLayers();*/

        //points = obliqueDeposition(25, 768, 800, reps, 0, 0, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        //params.clearLayers();
        //points = obliqueDeposition(89, 768, 800, reps / 2, 0, 3, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, true, 85, Acceleration::NONE);
        //params.clearLayers();
        //points = obliqueDeposition(70, 768, 800, reps, 0, 0, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        //params.clearLayers();
        //points = obliqueDeposition(75, 768, 800, reps, 0, 0, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, 0, Acceleration::NONE);
        //params.clearLayers();

        /*points = obliqueDeposition(70, 768, 800, reps, 0, 0, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, true, 89, Acceleration::NONE);
        points = obliqueDeposition(86, 768, 800, reps / 4, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        points = obliqueDeposition(87, 768, 800, reps / 4, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        points = obliqueDeposition(88, 768, 800, reps / 4, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        points = obliqueDeposition(89, 768, 800, reps / 4, 0, 0, seed, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();*/

        //points = obliqueDeposition(88, 768, 800, reps / 2, 0, 0, 10001, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::DEC);
        //params.clearLayers();
        //points = obliqueDeposition(88, 768, 800, reps, 0, 0, 10002, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::DEC);
        //params.clearLayers();
        //points = obliqueDeposition(88, 768, 800, reps, 0, 0, 10003, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::DEC);
        //params.clearLayers();
        //points = obliqueDeposition(88, 768, 800, reps / 4, 0, 1000, 10000, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        //points = obliqueDeposition(88, 768, 800, reps / 2, 0, 0, 10001, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::DEC);
        //params.clearLayers();
        //points = obliqueDeposition(88, 768, 800, reps / 4, 0, 1000, 10000, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        //points = obliqueDeposition(88, 768, 800, reps, 0, 0, 10002, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::DEC);
        //params.clearLayers();
        //points = obliqueDeposition(88, 768, 800, 1048576, 0, 0, 10000, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::DEC);
        //points = obliqueDeposition(88, 768, 800, reps, 0, 0, 10003, 10, &sSi, &spread, &Si, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::DEC);
        //params.clearLayers();




        //for (int d = 0; d < Ds.size(); d++) {
        //    for (int t = 0; t < thetas.size(); t++) {
        //        if (thetas[t] < 89 && !(Ds[d] == 1 || Ds[d] == 2 || Ds[d] == 5)) {
        //            continue;
        //        }
        //        for (int s = 0; s < 1; s++) {
        //            spread[0] = spreads[s];
        //            spread[1] = spreads[s];
        //            std::cout << "Simulating theta " << thetas[t] << " degrees and " << Ds[d] << " diffusion steps." << std::endl;
        //            points = obliqueDeposition(thetas[t], L, H, reps, 0, 0, seed, Ds[d], &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        //            params.clearLayers();
        //            //points = obliqueDeposition(thetas[t], L, H, reps, 0, 0, seed, Ds[i], &sZrO2, &spread, &ZrO2, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, yZrO2, false, false, Acceleration::NONE);
        //            //params.clearLayers();
        //        }
        //    }
        //}
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

        //
        //for (int n = 0; n < 5; n++) {
        //    for (int d = 0; d < Ds.size(); d++) {
        //        for (int t = 0; t < thetas.size(); t++) {
        //            for (int s = 0; s < 1; s++) {
        //                spread[0] = spreads[s];
        //                spread[1] = spreads[s];
        //                std::cout << "Simulating theta " << thetas[t] << " degrees and " << Ds[d] << " diffusion steps." << std::endl;
        //                points = obliqueDeposition(thetas[t], L, H, reps, 0, 0, seed, Ds[d], &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        //                params.clearLayers();
        //                //points = obliqueDeposition(thetas[t], L, H, reps, 0, 0, seed, Ds[i], &sZrO2, &spread, &ZrO2, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, yZrO2, false, false, Acceleration::NONE);
        //                //params.clearLayers();
        //            }
        //        }
        //    }
        //}



        /*points = obliqueDeposition(theta, L, H, reps, 0, 2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 0, -2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 90, 2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 90, -2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 180, 2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 180, -2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 270, 2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 270, -2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 45, 2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 45, -2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 135, 2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 135, -2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 225, 2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 225, -2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 315, 2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();
        points = obliqueDeposition(theta, L, H, reps, 315, -2, seed, 10, &sSi, &spread, &Si, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        params.clearLayers();*/
        /*points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 10, &sSi, &spread, &SiAg, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps/8, 0, 0, seed, 30, &sAg, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 10, &sSi, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps/8, 0, 0, seed, 30, &sAg, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        params.clearLayers();*/
        /*points = obliqueDeposition(theta, L, H, reps, 0, 0.85, seed, 10, &sSi, &spread, &SiAg, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps/8, 0.85*360, 0.15, seed, 30, &sAg, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps, 0, 0.85, seed, 10, &sSi, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps/8, 0.85 * 360, 0.15, seed, 30, &sAg, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        params.clearLayers();*/
        /*points = obliqueDeposition(theta, L, H, reps, 90, 0, seed, 10, &sSi, &spread, &SiAg, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps/8, 90, 0, seed, 30, &sAg, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps, 90, 0, seed, 10, &sSi, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps/8, 90, 0, seed, 30, &sAg, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        params.clearLayers();*/
        /*points = obliqueDeposition(theta, L, H, reps, 0, -0.85, seed, 10, &sSi, &spread, &SiAg, inputGrid, inputGridPoints, outGrid, 0, 0, stepper_resolution, &params, ySi, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps/8, -0.85 * 360, -0.15, seed, 30, &sAg, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps, 0, -0.85, seed, 10, &sSi, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);
        points = obliqueDeposition(theta, L, H, reps/8, -0.85 * 360, -0.15, seed, 30, &sAg, &spread, &SiAg, *outGrid, points, outGrid, 0, 0, stepper_resolution, &params, ySiAg, false, false, Acceleration::NONE);*/
        //points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 5, &species, &spread, &ZrO2, points, *outGrid, outGrid, 0, 0, stepper_resolution, &params, system, false, false, Acceleration::NONE);
        //for (int p = 1; p < 2; p++) {
        //    std::cout << "Starting simulation " << p << "." << std::endl;
        //    uint32_t points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 5, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::ACC);
        //    params.clearLayers();
        //    points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 5, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::DEC);
        //    params.clearLayers();
        //} // Si ACC DEC D5
        //
        //for (int p = 1; p < 2; p++) {
        //    std::cout << "Starting simulation " << p << "." << std::endl;
        //    uint32_t points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 10, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::ACC);
        //    params.clearLayers();
        //    points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 10, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::DEC);
        //    params.clearLayers();
        //} // Si ACC DEC D10

        //for (int p = 1; p < 2; p++) {
        //    std::cout << "Starting simulation " << p << "." << std::endl;
        //    uint32_t points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 5, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::ACC);
        //    points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 5, &species, &spread, &weights, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::DEC);
        //    params.clearLayers();
        //} // Si BICONE D5

        //for (int p = 1; p < 2; p++) {
        //    std::cout << "Starting simulation " << p << "." << std::endl;
        //    uint32_t points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 10, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::ACC);
        //    points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 10, &species, &spread, &weights, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::DEC);
        //    params.clearLayers();
        //} // Si BICONE D10

        //for (int p = 1; p < 2; p++) {
        //    std::cout << "Starting simulation " << p << "." << std::endl;
        //    uint32_t points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 5, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::DEC);
        //    points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 5, &species, &spread, &weights, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::ACC);
        //    params.clearLayers();
        //} // Si BICONE D5

        //for (int p = 1; p < 2; p++) {
        //    std::cout << "Starting simulation " << p << "." << std::endl;
        //    uint32_t points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 10, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::DEC);
        //    points = obliqueDeposition(theta, L, H, reps, phi, turns, seed, 10, &species, &spread, &weights, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false, Acceleration::ACC);
        //    params.clearLayers();
        //} // Si BICONE D10
        /*
        for (int p = 1; p < 17; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, phi, 1, seed, 5, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            params.clearLayers();
        } // Si sweep D5 helix

        for (int p = 1; p < 17; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, phi, 1, seed, 10, &species, &spread, &weights, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            params.clearLayers();
        } // Si sweep D10 helix



        std::vector<int8_t> speciesAg = { 2 };

        system = "Si_Ag";

        for (int p = 1; p < 1; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 5, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0, 0, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false);

            params.clearLayers();
        } // Si_Ag no sweep D5

        for (int p = 1; p < 1; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 10, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0, 0, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false);

            params.clearLayers();
        } // Si_Ag no sweep D10

        for (int p = 2; p < 3; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 5, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps/8, 0, 0, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false);

            params.clearLayers();
        } // Si_Ag sweep Si D5

        for (int p = 2; p < 3; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 5, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0, 0, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);

            params.clearLayers();
        } // Si_Ag sweep both D5

        for (int p = 2; p < 3; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 10, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0, 0, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false);

            params.clearLayers();
        } // Si_Ag sweep Si D10

        for (int p = 2; p < 3; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, 0, seed, 10, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0, 0, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);

            params.clearLayers();
        } // Si_Ag sweep both D10

        for (int p = 1; p < 17; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, .85, seed, 5, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0.85*360, .15, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false);

            params.clearLayers();
        } // Si_Ag sweep Si D5 helix

        for (int p = 1; p < 17; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, .85, seed, 5, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0.85 * 360, .15, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);

            params.clearLayers();
        } // Si_Ag sweep both D5 helix

        for (int p = 1; p < 17; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, .85, seed, 10, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0.85 * 360, .15, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, false, false);

            params.clearLayers();
        } // Si_Ag sweep Si D10 helix

        for (int p = 1; p < 17; p++) {
            std::cout << "Starting simulation " << p << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, .85, seed, 10, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 0.85 * 360, .15, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, p * 4, 45, stepper_resolution, &params, system, true, false);

            params.clearLayers();
        } // Si_Ag sweep both D10 helix


        for (int n = 0; n < 16; n++) {
            std::cout << "Starting simulation " << n + 1 << "." << std::endl;
            uint32_t points = obliqueDeposition(theta, L, H, reps, 0, turns*.85, seed, diffusion_steps + n, &species, &spread, &weights_SiAg, inputGrid, inputGridPoints, outGrid, stepper_resolution, &params, system, true, false);
            points = obliqueDeposition(theta, L, H, reps / 8, 360*.85, turns*.15, seed, 30, &speciesAg, &spread, &weights_SiAg, *outGrid, points, outGrid, stepper_resolution, &params, system, true, false);

            params.clearLayers();
        }
        */
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
