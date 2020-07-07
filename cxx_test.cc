#include <stdlib.h>
#include <cxxopts.hpp>

int main(int argc, char** argv) {

    int alpha = 0;
    double beta = 0;
    bool ferro = false;

    cxxopts::Options options("TestProgram", "Tests the cxxopts header");

    options.add_options()
        ("a,alpha", "value of alpha", cxxopts::value<int>())
        ("b,beta", "value of beta", cxxopts::value<double>())
        ("f,ferro", "ferro or not");

    auto result = options.parse(argc, argv);

    alpha = result["a"].as<int>();
    beta = result["b"].as<double>();
    ferro = result["f"].as<bool>();

    printf("alpha = %d, beta = %f, %d\n", alpha, beta, ferro);
    
    }