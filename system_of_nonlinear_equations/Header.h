#pragma once
#include <vector>
#include<iostream>
#include <cstdlib> 
#include <ctime>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <windows.h>
#include <shellapi.h>
#include <stdio.h> 
#include <iostream> 
#include <vector> 
#include <fstream> 
#include <strstream> 
#include <sstream> 
#include <cassert> 
#include <set> 
#include <map> 
#include <iomanip> 
#include <iterator> 
#include <algorithm> 
#include <list> 
#include <string> 
#include <cmath> 
#include <random>
using namespace std;

double f1(double x, double y);
double f2(double x, double y);
double f_x(double y);
double f_y(double x);
double f1_deriv_x(double x, double y);
double f2_deriv_x(double x, double y);
double f1_deriv_y(double x, double y);
double f2_deriv_y(double x, double y);
double l2(vector<double> x);
double l1(vector<double> x);
double cheb(vector<double> x);
vector<double> simple_it(double x0, double y0);
vector<double> grad(vector <vector<double>> a, vector<double> b);
vector<double> matrix_mult(vector<vector<double>> a, vector<double> b);
vector<double> newton(double x0, double y0);
double scalar_mult(vector<double> a, vector<double> b);
vector<double> gauss2(vector<vector<double>> a, vector<double> y);
void decompose(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u);
vector<double> lu(vector <vector<double>> a, vector<double> b);
void show(vector<vector<double>> a);
int sign(double x);