#pragma once
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>


double Rastrigin(std::vector<double> v);
double Michalewicz(std::vector<double> v);
double DeJong(std::vector<double> v);
double Schwefel(std::vector<double> v);