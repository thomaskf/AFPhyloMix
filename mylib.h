#ifndef _MYLIB_
#define _MYLIB_

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <cstring>
#include <time.h>
#include <float.h>

#define MAXDIM 10000000

using namespace std;

static string i2s[10] = {"0","1","2","3","4","5","6","7","8","9"};

static long double cumlogI[MAXDIM];

static bool initialCumLogI = false;

string int2str(int i);

void tokenizer(string seq, string separators, vector<string>* result);

void randomizePos(vector<int>& arr, int n);

void randomizePos(vector<int>& arr);

void statistics(long double* arr, int n, long double& min, long double& max, long double& avg, long double& median, long double& stdev);

void statistics(vector<long double>& arr, long double& min, long double& max, long double& avg, long double& median, long double& stdev);

// return the value of topK
int getTopK(vector<int>& arr, int topK); 

long combin(int n, int k);

long fact(int n);

void initializeCumLogI();

long double logfact(int n);

long double logCombin(int n, int k);

long double logBinomial(int x, int n, long double p, bool& isUndefined);

// (1) if p=0 and x>0, return undefine
// (2) if p=1 and x<n, return undefine
long double logBinomial_noCoeff(int x, int n, long double p, bool& isUndefined);

void ratios(vector<long double>& in_logf, vector<long double>& out_ratiof);

// remove all spaces inside the string
void removeSpaces(string& s);

// check whethe the string contains space
bool containSpace(string& s);

// is included each other
// 0: no; 1: s1 = s2; 2: s1 is included by s2; 3: s2 is included by s1
int isIncluded(string& s1, string& s2);

// compute the BIC value
long double BIC(long double logL, int degFree, int dataSize);

// compute the AIC value
long double AIC(long double logL, int degFree, int dataSize);

// compute the Adjusted BIC value
long double ABIC(long double logL, int degFree, int dataSize);

// compute the CAIC value
long double CAIC(long double logL, int degFree, int dataSize);

// compute the AICc value
long double AICc(long double logL, int degFree, int dataSize);

// double to string
string doublToStr(double d, double decimalPlace);

// round the double to a certain number of decimal place
double roundNumber(double d, int decimalPlace);

// compute log(x+y), given log(x) and log(y)
long double log_x_plus_y(long double logx, long double logy);

// compute the log value of dirichlet multinomial distribution
// i.e. Log ( Prob(beta | alpha) )
// the sizes of alpha and beta are the same
// all elements in alpha has to be greater than one
// need to have invoked initializeCumLogI()
long double logDiriMultDist(vector<int>& alpha, vector<int>& beta);

// compute the log value of multinomial distribution
// i.e. log ( prob(freq | p) )
// the sizes of freq and p are the same
long double logMultiDist(vector<int>& freq, vector<double>& p);

// compute the log value of multinomial coefficient
long double logMultiCoeff(vector<int>& freq);

// compute the log value of multinomial coefficient
long double logMultiCoeff(int* freq, int n);

int minInt(int x, int y);
int maxInt(int x, int y);
// get maximum of three integers
int maxInt(int x, int y, int z);
// get minimum item out of the first K itemss inside doublelist
// numZero - number of items are zeros
// minNoZero - the minimum non-zerio item
void minDouble(vector<double>& doubleList, int& numZero, double& minNoZero, int K);

// use continue function to extimate the maximum of a set of long double
// larger value of k_factor, closer to the maximum value
// long double maxLongFun(vector<long double>& values, int k_factor);

// same as the previous function
// use continue function to extimate the maximum of a set of long double
// larger value of k_factor, closer to the maximum value
// except that all the input values are log values
long double logmaxLongFun(vector<long double>& logvalues, int k_factor);

// get indices with maximum values
void getMaxIndices(vector<long double>& values, vector<int>& maxInd);

// build cigar string
// input: alignS - query alignment
//        alignT - genome alignment
// if the whole read cannot be aligned, return false
bool buildCigar(string& alignS, string& alignT, string& cigarStr);

// if x1 > x2, swap between them
// else no change
void reorder(int& x1, int& x2);

// edit distance between two strings
int editDist(string& s1, string& s2);

// hamming distance between two strings
int hamDist(string& s1, string& s2);

// select samples according to weights
void randSelect(vector<double>& weights, vector<int>& selects, int num);

void timeElapsed(clock_t fr_t, clock_t to_t, string msg);

// for Kahan summation algorithm
class Summation {
private:
    long double s; // sum
    long double y, t, c;
    bool islocked;
public:
    Summation(); // constructor
    void reset();
    void add(long double input);
    long double sum();
    void lock(); // lock it so that no other function can use it at the same time
    void unlock(); // release it so that other function can use it now.
};

// get the double after the position p
// return the double, and the starting and the ending position of the double
double getDouble(string& str, int p, int& fr_p, int& to_p);

// generate a set of random numbers
void genRandom(double fr_value, double to_value, int n, vector<double>& rand_nums);
void genRandom(double fr_value, double to_value, double sum_value, int n, vector<double>& rand_nums);

#endif
