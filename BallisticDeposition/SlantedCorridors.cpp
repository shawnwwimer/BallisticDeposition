#include "SlantedCorridors.h"
#include <set>
#include <algorithm>
#include <iterator>



std::set<particle_priority>* SlantedCorridors::find_bins(std::array<float, 3>* position, float radius)
{
    std::set<particle_priority>* bin = find_bin(position);
    bins_found[0] = bin;
    bins_found_num = 1;

    bool found = false;

    // Order chosen to minimize individual parameter changes
    // --
    std::array<float, 3> corner = { modulof((*position)[0] - radius, L), modulof((*position)[1] - radius, L), (*position)[2] };
    bin = find_bin(&corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // +-
    corner[0] = modulof((*position)[0] + radius, L);
    bin = find_bin(&corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // ++
    corner[1] = modulof((*position)[1] + radius, L);
    bin = find_bin(&corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // -+
    corner[0] = modulof((*position)[0] - radius, L);
    bin = find_bin(&corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    /*
    // ---
    std::array<float, 3>corner = {modulof(position[0] - radius, L), modulof(position[1] - radius, L), position[2] - radius};
    bin = find_bin(corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;


    // +--
    corner[0] = modulof(position[0] + radius, L);
    bin = find_bin(corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // +-+
    corner[2] = position[2] + radius;
    bin = find_bin(corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // +++
    corner[1] = modulof(position[1] + radius, L);
    bin = find_bin(corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // ++-
    corner[2] = position[2] - radius;
    bin = find_bin(corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // -+-
    corner[0] = modulof(position[0] - radius, L);
    bin = find_bin(corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // -++
    corner[2] = position[2] + radius;
    bin = find_bin(corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;

    // --+
    corner[1] = modulof(position[1] - radius, L);
    bin = find_bin(corner);
    for (int i = 0; i < bins_found_num; i++) {
        if (bins_found[i] == bin) {
            found = true;
        }
    }
    if (!found) {
        bins_found[bins_found_num] = bin;
        bins_found_num += 1;
    }
    found = false;
    */

    return nullptr;
}

std::set<particle_priority>* SlantedCorridors::find_neighborhood(std::array<float, 3>* position) {
    uint32_t idx = find_bin_idx(position);
    bins_found[0] = &bins[idx];
    uint32_t row = idx / bins_on_side;
    uint32_t column = idx % bins_on_side;
    uint32_t left = modulo(column % bins_on_side - 1, bins_on_side) + bins_on_side * row;
    uint32_t right = (column % bins_on_side + 1) % bins_on_side + bins_on_side * row;
    bins_found[1] = &bins[left];
    bins_found[2] = &bins[right];
    bins_found[3] = &bins[modulo(idx + bins_on_side, bins_num)];
    bins_found[4] = &bins[modulo(idx - bins_on_side, bins_num)];
    bins_found[5] = &bins[modulo(left + bins_on_side, bins_num)];
    bins_found[6] = &bins[modulo(left - bins_on_side, bins_num)];
    bins_found[7] = &bins[modulo(right + bins_on_side, bins_num)];
    bins_found[8] = &bins[modulo(right - bins_on_side, bins_num)];
    bins_found_num = 9;

    return nullptr;
}

// Atoms is (nx6) float array
collision_description* SlantedCorridors::drop_particle(std::array<float, 3>* position, float radius, std::vector<std::vector<float>>* atoms)
{
    // Combine bins into temporary big bin
    find_neighborhood(position);
    //std::set<particle_priority>** fbins = bins_found;
    //std::set<particle_priority> b;
    //for (int i = 0; i < bins_found_num; i++) {
    //    std::set_union(std::begin(b), std::end(b), std::begin(*bins_found[i]), std::end(*bins_found[i]), std::inserter(b, std::begin(b)));
    //}

    // Need to know relative position
    float adjusted_x = (*position)[2] * tan_theta + (*position)[0];
    int factor = (int)(adjusted_x / L);
    adjusted_x = modulof(adjusted_x, Lfloat);

    // If adjusted_x too near L we need to set 'straddle' high
    bool straddle = false;
    if (adjusted_x < radius || adjusted_x > L - radius) {
        straddle = true;
    }

    // Get sizes and set up indices for bins found
    std::vector<std::set<particle_priority>::iterator> iterators;
    int total = 0;
    for (int i = 0; i < bins_found_num; i++) {
        iterators.push_back(bins_found[i]->begin());
        total += bins_found[i]->size();
    }

    for (int i = 0; i < total; i++) {
        // find bin with top particle
        float max = 0;
        int bin = -1;
        for (int j = 0; j < bins_found_num; j++) {
            // can't read past end of bin
            if (iterators[j] != bins_found[j]->end()) {
                float temp_pri = particles[(*iterators[j]).idx]->priority;
                // if the particle is near the PB we need to adjust some of the priorities as we go
                if (straddle) {
                    float adj_x_p = particles[(*iterators[j]).idx]->adjx;
                    if (straddle && adj_x_p > L - bin_size) {
                        temp_pri += L;
                    }
                }

                if (temp_pri > max) {
                    bin = j;
                    max = particles[(*iterators[j]).idx]->priority;
                }
            }
        }

        if (bin == -1) {
            break;
        }

        // get top particle from bin and increment the iterator
        particle_priority p = *iterators[bin];
        iterators[bin]++;

        int idx = p.idx;//atoms.size() - p.idx - 1;
        //std::vector<float> particle = atoms[idx];

        float rs = radius + (*atoms)[idx][4];


        float x_on_particle_z = ((*position)[2] - (*atoms)[idx][2]) * tan_theta + (*position)[0];
        int factor = (int)(x_on_particle_z / L);
        x_on_particle_z = fmod(x_on_particle_z, Lfloat);

        float offset[2] = { 0 };

        if ((*atoms)[idx][0] - x_on_particle_z < -Lfloat / 2.0) {
            offset[0] = -Lfloat;
        }
        else if ((*atoms)[idx][0] - x_on_particle_z > Lfloat / 2.0) {
            offset[0] = Lfloat;
        }
        if ((*atoms)[idx][1] - (*position)[1] < -Lfloat / 2.0) {
            offset[1] = -Lfloat;
        }
        else if ((*atoms)[idx][1] - (*position)[1] > Lfloat / 2.0) {
            offset[1] = Lfloat;
        }
        //float d = sqrt(pow(sin_theta * ((*position)[2] - (*atoms)[idx][2]) - cos_theta * ((*position)[0] - (*atoms)[idx][0] + offset[0]), 2) + pow((*position)[1] - (*atoms)[idx][1] + offset[1], 2));

        // Find distance squared value
        float radii2 = rs * rs;
        float rotx = ((*position)[0] - factor * Lfloat - (*atoms)[idx][0] + offset[0]) * cos_theta + ((*position)[2] - (*atoms)[idx][2]) * sin_theta;
        float d2 = rotx * rotx + ((*position)[1] - (*atoms)[idx][1] + offset[1]) * ((*position)[1] - (*atoms)[idx][1] + offset[1]);

        if (d2 < rs*rs) {
            // Calculate position
            float rotz = sqrt(radii2 - d2);
            float z = rotx * sin_theta + rotz * cos_theta + (*atoms)[idx][2];
            float x = rotx * cos_theta - rotz * sin_theta + (*atoms)[idx][0];
            //float z = (sqrt(pow(rs, 2) - pow((*position)[1] - (*atoms)[idx][1] + offset[1], 2)) + (*atoms)[idx][2] * sin_theta + cos_theta * ((*position)[0] + (*position)[2] * tan_theta)) / (sin_theta + cos_theta * tan_theta);
            //float x = (*position)[0] + ((*position)[2] - z) * tan_theta;

            if (z < radius) {
                continue; // Would only collide below the substrate
            }
            collision.position = { modulof(x, L), (*position)[1], z };
            collision.idx = p.idx;//atoms.size() - p.idx - 1;
            //if ((x - (*atoms)[idx][0]) * (x - (*atoms)[idx][0]) + ((*position)[1] - (*atoms)[idx][1] + offset[1]) * ((*position)[1] - (*atoms)[idx][1] + offset[1]) + (z - (*atoms)[idx][2]) * (z - (*atoms)[idx][2])-.001 > radii2) {
            //    std::cout << "Too far!" << std::endl;
            //}
            return &collision;
        }
    }
    float x = (*position)[0] - radius * tan_theta;

    collision.position = { modulof(x, L), (*position)[1], radius };
    collision.idx = -1;
    return &collision;
}
