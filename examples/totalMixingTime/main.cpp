// system includes
#include <iostream>

// marathon includes
#include "marathon/mixing_time.h"
#include "marathon/binary_matrix/fixed_margin/switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/edge_switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/curveball.h"
#include "marathon/binary_matrix/fixed_margin/sys_curveball.h"

// auxiliary functions
#include "helper.h"

int main(int argc, char **argv) {

    // parse command line arguments
    if (!parse_arguments(argc, argv)) {
        print_help_message();
        return -1;
    }

    // create Markov chain
    std::unique_ptr<marathon::MarkovChain> mc;
    switch (chain) {
        case classical_switch:
            mc = std::make_unique<marathon::binary_matrix::fixed_margin::SwitchChain>(inst);
            break;
        case edge_switch:
            mc = std::make_unique<marathon::binary_matrix::fixed_margin::EdgeSwitchChain>(inst);
            break;
        case curveball:
            mc = std::make_unique<marathon::binary_matrix::fixed_margin::Curveball>(inst);
            break;
        case sys_curveball:
            mc = std::make_unique<marathon::binary_matrix::fixed_margin::SysCurveball>(inst);
            break;
        default:
            std::cerr << "Undefined behavior: CHAIN unknown: " << argv[1] << std::endl;
            return -1;
    }

    // construct state graph
    marathon::StateGraph sg(*mc);

    // determine number of states
    const size_t N = sg.getNumStates();
    // calculate total mixing time
    marathon::MixingTimeCalculator<double> mtc(sg);
    int t;
    if(N < 10000) {
        // more efficient but requires about 32 N^2 byte of main memory
        t = mtc.totalMixingTimeDense(eps);
    }
    else {
        // larger running time but less memory intense
        t = mtc.totalMixingTime(eps);
    }

    std::cout << "ds        " << inst << "\n";
    std::cout << "mixalgo   " << std::string(argv[1]) << "\n";
    std::cout << "eps       " << eps << "\n";
    std::cout << "numstates " << N << "\n";
    std::cout << "tmixtime  " << t << "\n";

    return 0;
}
