#include <vector>
#include "prc1System.h"


PRC1System::PRC1System(double microtubule_length, double site_spacing, double offset,  double spring_constant, double rest_length, double k_B_T, double microtubule_seperation) {
    microtubule_offset = offset;
    num_agents = 0;
    num_heads_attached_top = 0;
    num_heads_attached_bottom = 0;

    top_sites_are_taken = {false};
    bottom_sites_are_taken = {false};
    sites = {};
    double cur_site = 0;
    while (cur_site < microtubule_length) {
        sites.push_back(cur_site);
        top_sites_are_taken.push_back(false);
        bottom_sites_are_taken.push_back(false);
        cur_site += site_spacing;
    }
    num_sites = sites.size();

    PRC1System::site_spacing = site_spacing;
    PRC1System::spring_constant = spring_constant;
    PRC1System::rest_length = rest_length;
    PRC1System::k_B_T = k_B_T;
    PRC1System::microtubule_seperation = microtubule_seperation;
}