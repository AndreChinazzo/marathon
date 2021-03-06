//
// Created by rechner on 3/1/18.
//

#ifndef MARATHON_HELPER_H
#define MARATHON_HELPER_H

#include <iostream>
#include <fstream>
#include <sys/stat.h>


// From https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool fileExists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

// Markov chains
enum chain_t {
    classical_switch,
    edge_switch,
    curveball,
    sys_curveball,
    quasi_curveball,
    quasi_curveball2,
    quasi_curveball3,
    quasi_curveball4,
    quasi_curveball5,
    quasi_curveball6,
} chain;

std::vector<std::string> insts;       // vector of strings encodeding problem instance
double eps = -1;        // maximal distance to uniform distribution (error term)

/**
 * Parse command line arguments and extract parameters.
 * @param argc Number of arguments.
 * @param argv Argument list.
 * @return True, if arguments are specified correctly. False, otherwise.
 */
bool parse_arguments(int argc, char **argv) {
    if (argc != 4)
        return false;

    // parse problem type
    if (strcmp(argv[1], "classical-switch") == 0)
        chain = classical_switch;
    else if (strcmp(argv[1], "edge-switch") == 0)
        chain = edge_switch;
    else if (strcmp(argv[1], "curveball") == 0)
        chain = curveball;
    else if (strcmp(argv[1], "sys-curveball") == 0)
        chain = sys_curveball;
    else if (strcmp(argv[1], "quasi-curveball") == 0)
        chain = quasi_curveball;
    else if (strcmp(argv[1], "quasi-curveball2") == 0)
        chain = quasi_curveball2;
    else if (strcmp(argv[1], "quasi-curveball3") == 0)
        chain = quasi_curveball3;
    else if (strcmp(argv[1], "quasi-curveball4") == 0)
        chain = quasi_curveball4;
    else if (strcmp(argv[1], "quasi-curveball5") == 0)
        chain = quasi_curveball5;
    else if (strcmp(argv[1], "quasi-curveball6") == 0)
        chain = quasi_curveball6;
    else {
        std::cerr << "Unknown CHAIN specifier: " << argv[1] << std::endl;
        return false;
    }

    eps = atof(argv[2]);

    if ( fileExists(std::string(argv[3])) ) {
        std::ifstream instsFile( argv[3] );
        std::string inst;
        while (std::getline(instsFile, inst)) {
            insts.push_back( inst );
        }
    }
    else // Single instance
    {
        insts.push_back(std::string(argv[3]));
    }
    return true;
}

void print_help_message() {
    std::cout << "Usage: totalMixingTime CHAIN EPSILON INSTANCE" << std::endl;
    std::cout << "Calculate the total mixing time of a specified Markov CHAIN on a given INSTANCE.\n";
    std::cout << "The total mixing time is the number of steps the Markov CHAIN is requires to\n";
    std::cout << "reach a probability distribution whose total variation distance to the uniform\n";
    std::cout << "distribution is at most EPSILON." << std::endl;
    std::cout << std::endl;
    std::cout << "The parameter CHAIN must be one of the following: " << std::endl;
    std::cout << "  'classical-switch':" << std::endl;
    std::cout << "       Markov chain defined by 'Kannan et al. Simple Markov-chain  algorithms" << std::endl;
    std::cout << "       for generating bipartite graphs and tournaments. Random Structures and" << std::endl;
    std::cout << "       Algorithms 14 (1997), 293–308'." << std::endl;
    std::cout << "  'edge-switch':" << std::endl;
    std::cout << "       Variant of the Markov-chain suggested by 'Kannan et al.' based on an" << std::endl;
    std::cout << "       informed edge selection at the cost of a larger memory consumption." << std::endl;
    std::cout << "  'curveball':" << std::endl;
    std::cout << "       Markov chain defined by 'Strona et al. A fast and unbiased procedure to" << std::endl;
    std::cout << "       randomize ecological binary matrices with fixed row and column totals." << std::endl;
    std::cout << "       Nature communications 5 (2014).'" << std::endl;
    std::cout << "  'sys-curveball':" << std::endl;
    std::cout << "       TODO" << std::endl;
    std::cout << "  'quasi-curveball':" << std::endl;
    std::cout << "       TODO" << std::endl;
    std::cout << "  'split-quasi-curveball':" << std::endl;
    std::cout << "       TODO" << std::endl;
    std::cout << std::endl;
    std::cout << "EPSILON must be a floating point number in the open interval (0,1).\n" << std::endl;
    std::cout << "The parameter INSTANCE is a string-encoded input instance of the form \"r*;c*\"." << std::endl;
    std::cout << "While the i-th r defines the sum of row i, the j-th c is the sum of column j." << std::endl;
    std::cout << "For example, the instance \"2,2,2;1,2,1,2\" corresponds to the" << std::endl;
    std::cout << std::endl;
    std::cout << "        row sums:    (2,2,2)" << std::endl;
    std::cout << "        column sums: (1,2,1,2)" << std::endl;
    std::cout << std::endl;
}

#endif //MARATHON_HELPER_H
