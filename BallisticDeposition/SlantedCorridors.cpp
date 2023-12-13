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
    int32_t idx = find_bin_idx(position);
    bins_found[0] = &bins[idx];
    int32_t row = idx / bins_on_side;
    int32_t column = idx % bins_on_side;
    int32_t left = modulo(column % bins_on_side - 1, bins_on_side) + bins_on_side * row;
    int32_t right = (column % bins_on_side + 1) % bins_on_side + bins_on_side * row;
    int32_t top = modulo(idx + bins_on_side, bins_num);
    int32_t bottom = modulo(idx - bins_on_side, bins_num);
    bins_found[1] = &bins[left];
    bins_found[2] = &bins[right];
    bins_found[3] = &bins[top];
    bins_found[4] = &bins[bottom];
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
    //float adjusted_x = (*position)[2] * tan_theta + (*position)[0];
    //adjusted_x = modulof(adjusted_x, Lfloat);

    // If adjusted_x too near L we need to set 'straddle'
    int straddle = 0;
    if ((*position)[0] < bin_size) {
        straddle = -1;
    }
    else if ((*position)[0] >= L - bin_size) {
        straddle = 1;
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
                if (straddle != 0) {
                    float adj_x_p = particles[(*iterators[j]).idx]->adjx;
                    if (straddle == 1 && adj_x_p < bin_size) {
                        temp_pri -= L * sin_theta;
                    }
                    else if (straddle == -1 && adj_x_p >= L - bin_size) {
                        temp_pri += L * sin_theta;
                    }
                }

                if (temp_pri > max) {
                    bin = j;
                    max = temp_pri;
                }
            }
        }

        if (bin == -1) {
            break;
        }

        // get top particle from bin and increment the iterator
        particle_priority p = *iterators[bin];
        iterators[bin]++;
        int idx = p.idx;

        // Distance to determine collision
        float rs = radius + (*atoms)[idx][4];

        // Get nearest y-distance between particles
        float projy = (*position)[1] - (*atoms)[idx][1];
        if (projy > Lfloat / 2) {
            projy -= L;
        }
        else if (projy < -Lfloat / 2) {
            projy += L;
        }
        
        // If the distance is too large we can skip this particle
        if (projy > rs || -projy < -rs) {
            continue;
        }

        float x_on_particle_z = ((*position)[2] - (*atoms)[idx][2]) * tan_theta + (*position)[0];
        x_on_particle_z = modulof(x_on_particle_z, Lfloat);

        // Find distance squared value
        float radii2 = rs * rs;

        // Get nearest position of incoming particle respecting PBC
        float projx = x_on_particle_z - (*atoms)[idx][0];
        if (projx > Lfloat / 2) {
            projx -= L;
        }
        else if (projx < -Lfloat / 2) {
            projx += L;
        }
        projx += (*atoms)[idx][0];
        projy += (*atoms)[idx][1];
        
        // Now get the projection of the stationary particle onto the line of travel
        // projy doesn't change because the line has constant y-coordinate
        float dot = ((*atoms)[idx][0] - projx) * sin_theta;// +((*atoms)[idx][2] - projz) * -cos_theta;
        projx += dot * sin_theta;
        float projz = (*atoms)[idx][2] + dot * -cos_theta;

        // Find distance squared back to line
        float d2 = pow(projx - (*atoms)[idx][0], 2) + pow(projy - (*atoms)[idx][1], 2) + pow(projz - (*atoms)[idx][2], 2);

        if (d2 < rs*rs) {
            // Calculate distance displaced along line and then position
            float r = sqrt(radii2 - d2);
            float x = projx - r * sin_theta;
            float z = projz + r * cos_theta;

            if (z < radius) {
                continue; // Would first collide below the substrate
            }
            collision.position = { modulof(x, L), (*position)[1], z };
            collision.idx = p.idx;
            return &collision;
        }
    }
    float x = (*position)[0] - radius * tan_theta;

    collision.position = { modulof(x, L), (*position)[1], radius };
    collision.idx = -1;
    return &collision;
}
