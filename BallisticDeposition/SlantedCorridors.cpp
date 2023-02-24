#include "SlantedCorridors.h"
#include <set>
#include <algorithm>
#include <iterator>

float modulof(float val, float mod)
{
    return fmodf((fmodf(val, mod) + mod), mod);
}

std::set<particle_priority>* SlantedCorridors::find_bins(std::array<float, 3> position, float radius)
{
    std::set<particle_priority>* bin = find_bin(position);
    bins_found[0] = bin;
    bins_found_num = 1;
    
    bool found = false;

    // Order chosen to minimize individual parameter changes
    // ---
    std::array<float, 3>corner = { modulof(position[0] - radius, L), modulof(position[1] - radius, L), position[2] - radius };
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

    return nullptr;
}

// Atoms is (nx6) float array
collision_description* SlantedCorridors::drop_particle(std::array<float, 3> position, float radius, std::vector<std::array<float, 6>> atoms)
{
    // Combine bins into temporary big bin
    find_bins(position, radius);
    //std::set<particle_priority>** fbins = bins_found;
    std::set<particle_priority> b;
    for (int i = 0; i < bins_found_num; i++) {
        std::set_union(std::begin(b), std::end(b), std::begin(*bins_found[i]), std::end(*bins_found[i]), std::inserter(b, std::begin(b)));
    }
    
    // Need to know relative position
    float adjusted_x = position[2] * tan_theta + position[0];
    int factor = (int)(adjusted_x / L);
    adjusted_x = modulof(adjusted_x, Lfloat);
    

    for (auto p : b) {
        std::array<float, 6> particle = atoms[p.idx];

        float rs = radius + particle[4];

        float x_on_particle_z = (position[2] - particle[2]) * tan_theta + position[0];
        int factor = (int)(x_on_particle_z / L);
        x_on_particle_z = fmod(x_on_particle_z, Lfloat);

        float offset[3] = { 0 };

        if (particle[0] - x_on_particle_z < -Lfloat / 2.0) {
            offset[0] = -Lfloat;
        }
        else if (particle[0] - x_on_particle_z > Lfloat / 2.0) {
            offset[0] = Lfloat;
        }
        if (particle[1] - position[1] < -Lfloat / 2.0) {
            offset[1] = -Lfloat;
        }
        else if (particle[1] - position[1] > Lfloat / 2.0) {
            offset[1] = Lfloat;
        }

        // Find distance squared value
        float radii2 = rs * rs;
        float rotx = (position[0] - factor * Lfloat - particle[0] + offset[0]) * cos_theta + (position[2] - particle[2]) * sin_theta;
        float d2 = rotx * rotx + (position[1] - particle[1] + offset[1]) * (position[1] - particle[1] + offset[1]);

        if (radii2 > d2) {
            // Calculate position
            float rotz = sqrt(radii2 - d2);
            float z = rotx * sin_theta + rotz * cos_theta + particle[2];
            float x = rotx * cos_theta - rotz * sin_theta + particle[0];

            if (z < radius) {
                continue; // Would only collide below the substrate
            }
            collision.position = { x, position[1], z };
            collision.idx = p.idx;
            return &collision;
        }
    }
    float x = adjusted_x - radius * tan_theta;
    
    collision.position = { x, position[1], radius };
    collision.idx = -1;
    return &collision;
}
