// system includes
#include <iostream>

// marathon includes
#include "marathon/mixing_time.h"
#include "marathon/binary_matrix/fixed_margin/switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/edge_switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/curveball.h"
#include "marathon/binary_matrix/fixed_margin/sys_curveball.h"
#include "marathon/binary_matrix/fixed_margin/quasi_curveball.h"
#include "marathon/binary_matrix/fixed_margin/quasi_curveball2.h"
#include "marathon/binary_matrix/fixed_margin/quasi_curveball3.h"
#include "marathon/binary_matrix/fixed_margin/quasi_curveball4.h"
#include "marathon/binary_matrix/fixed_margin/quasi_curveball5.h"
#include "marathon/binary_matrix/fixed_margin/quasi_curveball6.h"

// auxiliary functions
#include "helper.h"

void runTMixTime(std::string &inst) {
    // create Markov chain
    std::unique_ptr<marathon::MarkovChain> mc;
    std::string chainName("");
    switch (chain) {
    case classical_switch:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::SwitchChain>(inst);
    chainName = "classical_switch";
    break;
    case edge_switch:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::EdgeSwitchChain>(inst);
    chainName = "edge_switch";
    break;
    case curveball:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::Curveball>(inst);
    chainName = "curveball";
    break;
    case sys_curveball:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::SysCurveball>(inst);
    chainName = "sys_curveball";
    break;
    case quasi_curveball:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::QuasiCurveball>(inst);
    chainName = "quasi_curveball";
    break;
    case quasi_curveball2:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::QuasiCurveball2>(inst);
    chainName = "quasi_curveball2";
    break;
    case quasi_curveball3:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::QuasiCurveball3>(inst);
    chainName = "quasi_curveball3";
    break;
    case quasi_curveball4:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::QuasiCurveball4>(inst);
    chainName = "quasi_curveball4";
    break;
    case quasi_curveball5:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::QuasiCurveball5>(inst);
    chainName = "quasi_curveball5";
    break;
    case quasi_curveball6:
    mc = std::make_unique<marathon::binary_matrix::fixed_margin::QuasiCurveball6>(inst);
    chainName = "quasi_curveball6";
    break;
    default:
    std::cerr << "Undefined behavior: CHAIN unknown" << std::endl;
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


    std::stringstream res("");
    res << "\"" << inst << "\"" << ";";
    res << chainName << ";";
    res << eps << ";";
    res << N << ";";
    res << t << "\n";
    std::cout << res.str();
}

int main(int argc, char **argv) {

    // parse command line arguments
    if (!parse_arguments(argc, argv)) {
        print_help_message();
        return -1;
    }

    std::cout << "ds;mixchain;eps;numstates;tmixtime" << "\n";

    for( auto &inst : insts ) {
        try {
            runTMixTime(inst);
        } catch (std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    return 0;
}
