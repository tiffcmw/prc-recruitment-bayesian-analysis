#ifndef PRC1SYSTEM_H
#define PRC1SYSTEM_H

#include <vector>

class PRC1System {
public:
    // positions where the top and bottom heads of each microtubule are attached.
    // positions are relative to the top and bottom microtubules respectively;
    // the left side of each microtubule is the 0 position
    std::vector<int> top_positions;
    std::vector<int> bottom_positions;
    
    // whether the top or bottom heads of the microtubule is attached
    std::vector<bool> top_is_attached;
    std::vector<bool> bottom_is_attached;

    // location in nm where each site actually is relative to the left of the respective microtubule
    std::vector<double> sites;

    // keeps track whether each of the top or bottom sites are taken
    std::vector<bool> top_sites_are_taken;
    std::vector<bool> bottom_sites_are_taken;

    // horizontal displacement of the top microtubule's left end from the bottom microtubule's left end (if positive, top is further right than bottom)
    double microtubule_offset;

    int num_sites;
    int num_agents;
    int num_heads_attached_top;
    int num_heads_attached_bottom;

    // relevant constants
    double site_spacing;
    double spring_constant;
    double rest_length;
    double k_B_T;
    double microtubule_seperation;

    PRC1System(double microtubule_length, double site_spacing, double offset, double spring_constant, double rest_length, double k_B_T, double microtubule_seperation);
    
};

#endif