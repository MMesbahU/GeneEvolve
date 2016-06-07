#include <ctime>
#include <cstdio>
#include <cmath>
#include <random>
#include "RasMatrix.h"
#include <vector>


namespace RasRandomNumber
{
    unsigned ras_now_nanoseconds(void);
    std::vector<std::vector<double> > ras_mvnorm(int n, std::vector<double> mu, std::vector<std::vector<double> > corr);
    std::vector<int> ras_rpois(int n, double lam);
    int ras_rpois(double lam);
    double ras_uniform(void);
    std::vector<int> ras_SampleWithoutReplacement(int populationSize, int sampleSize);
}
