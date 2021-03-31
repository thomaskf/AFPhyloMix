//
//  myRand.cpp
//  HaploCount
//
//  Created by Thomas Wong on 15/2/19.
//

#include "myRand.h"

// constructor
MyRand::MyRand(int inSeed) {

    seed = inSeed;
    
#ifdef USE_BOOST
    generator = boost::mt19937(seed);
    distribution = new uniform_real<>(0.0,1.0);
#else
    generator = default_random_engine(seed);
    distribution = new uniform_real_distribution<double>(0.0,1.0);
#endif

}

MyRand::MyRand() {
    seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
    // seed = 7777777;

#ifdef USE_BOOST
    generator = boost::mt19937(seed);
    distribution = new uniform_real<>(0.0,1.0);
#else
    generator = default_random_engine(seed);
    distribution = new uniform_real_distribution<double>(0.0,1.0);
#endif
}

// destructor
MyRand::~MyRand() {
    delete distribution;
}

// uniform distribution
// report a double between minValue and maxValue
double MyRand::runif(double minValue, double maxValue) {
#ifdef USE_BOOST
    uniform_real<> unif(minValue,maxValue);
#else
    uniform_real_distribution<double> unif(minValue,maxValue);
#endif
    return unif(generator);
}

// uniform distribution between 0.0 and 1.0
double MyRand::runif() {
    return distribution->operator()(generator);
}

// report an integer
int MyRand::iunif(int minValue, int maxValue) {
    if (minValue == maxValue)
        return minValue;
#ifdef USE_BOOST
    uniform_int<> unii(minValue,maxValue);
#else
    uniform_int_distribution<double> unii(minValue,maxValue);
#endif
    return unii(generator);
}
