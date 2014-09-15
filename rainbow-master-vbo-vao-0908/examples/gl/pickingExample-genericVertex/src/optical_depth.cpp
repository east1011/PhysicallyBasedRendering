//
//  optical_depth.h
//  Solar Position Algorithm
//
//  Created by jd on 2014. 2. 1..
//  Copyright (c) 2014 jd. All rights reserved.
//

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdio.h>


#include "optical_depth.h"

void get_waterDensity_and_waterDepth (int &n, double &L) {
    
    std::string input;
    // get numbers
    while (true) {
        std::cout << "Please enter number density (10000 ~ 100000): ";
        getline(std::cin, input);
        
        // This code converts from string to number safely.
        std::stringstream myStream(input);
        if (myStream >> n)
            break;
        std::cout << "Invalid number, please try again" << std::endl;
    }
    std::cout << "You entered: " << n << std::endl;
    
    while (true) {
        std::cout << "Please enter water depth (default: 5(m)): ";
        getline( std::cin, input);
        
        // This code converts from string to number safely.
        std::stringstream myStream(input);
        if (myStream >> L)
            break;
        std::cout << "Invalid number, please try again" << std::endl;
    }
    std::cout << "You entered: " << L << std::endl;
    
}

double calculate_optDepth (int n, double L) {
    
    double tau_N;
    double scatCross = 6.601272e-006;
    
    tau_N = scatCross * n * L;
    
    return tau_N;
}


