//
//  myRand.h
//  HaploCount
//
//  Created by Thomas Wong on 15/2/19.
//

#ifndef myRand_h
#define myRand_h

// comment the following line if boost is not installed in the machine
#define USE_BOOST

#ifdef USE_BOOST
#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
using namespace boost;
#else
#include <random>
#endif

#include <stdio.h>
#include <chrono>
#include <ctime>

using namespace std;

class MyRand {
public:

    // uniform random number generator
    unsigned seed;

#ifdef USE_BOOST
    boost::mt19937 generator;
    uniform_real<>* distribution;
#else
    std::default_random_engine generator;
    uniform_real_distribution<double>* distribution;
#endif
    
    // constructor
    MyRand();
    MyRand(int inSeed);
    
    // destructor
    ~MyRand();

    // uniform distribution
    // report a double between minValue and maxValue
    double runif(double minValue, double maxValue);
    
    // report an integer
    int iunif(int minValue, int maxValue);
    
    // uniform distribution between 0.0 and 1.0
    double runif();
};

#endif /* myRand_h */
