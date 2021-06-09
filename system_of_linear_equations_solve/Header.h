#pragma once

#include <iostream>
#include <vector>
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
#include<unordered_map>
using namespace std;
double LP(double x, vector<double> xv, vector<double> yv);
double f(double x);
void LocLag(vector<long double>* res, double a, double b, vector<double> x, vector<double> y, int n);
double l2(vector<double> x);
double l1(vector<double> x);
double cheb(vector<double> x);
vector<double> rand(double a, double b, int l);
double Lagr_basis(vector<double> x, double t, int alpha);
double scalar_mult(vector<double> a, vector<double> b);
vector<double> gauss(vector<vector<double>> a, vector<double> y);
void show(vector<vector<double>> a);
vector<double> lu(vector <vector<double>> a, vector<double> b);
void decompose(vector <vector <double>> a, vector <vector <double>>* l,vector <vector <double>>* u);
vector<double> relax(vector <vector<double>> a, vector<double> b, double w);
void relax2(vector <vector<double>> a, vector<double> b);
vector<double> grad(vector <vector<double>> a, vector<double> b);
vector<double> matrix_mult(vector<vector<double>> a, vector<double> b);
void proisv(vector <vector <double>> A, vector <vector <double>> B);
vector<double> lu_for_lent(vector <vector<double>> a, vector<double> b);
void decompose_lent(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u);
vector<double> relax_lent(vector <vector<double>> a, vector<double> b,double w);
vector<double> grad_lent(vector <vector<double>> a, vector<double> b);
void decompose3(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u);
void decompose_lent2(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u);
vector<double> gauss_lent(vector<vector<double>> a, vector<double> y);
vector<double> gauss2(vector<vector<double>> a, vector<double> y);
void decompose_lent3(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u);