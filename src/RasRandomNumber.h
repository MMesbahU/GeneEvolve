#include <ctime>
#include <cstdio>
#include <cmath>
#include <random>
#include "RasMatrix.h"
#include <vector>


namespace RasRandomNumber
{
    unsigned ras_now_nanoseconds(void);
    std::vector<std::vector<double> > ras_mvnorm(unsigned long int n, std::vector<double> &mu, std::vector<std::vector<double> > &corr, unsigned seed);
    std::vector<int> ras_rpois(unsigned long int n, double lam, unsigned seed);
    int ras_rpois(double lam, unsigned seed);
    double ras_uniform(unsigned seed);
    std::vector<unsigned long int> ras_SampleWithoutReplacement(unsigned long int populationSize, unsigned long int sampleSize, unsigned seed);
}


