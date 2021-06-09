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

double f1(double y1);
double f2(double x, double y, double y1);
double f2_temp(double x, double y, double y1);
double p(double x);
double q(double x);
double r(double x);
double f2_dy(double x, double y, double y1);
double f2_dy1(double x, double y, double y1);
double l2(vector<double> x);
double l1(vector<double> x);
double cheb(vector<double> x);
vector<double> grad(vector <vector<double>> a, vector<double> b);
vector<double> matrix_mult(vector<vector<double>> a, vector<double> b);
double scalar_mult(vector<double> a, vector<double> b);
void show(vector<vector<double>> a);
vector<double> sub(vector<double> a, vector<double> b);
void runge_kutt_4(vector<double> x, vector<double>* y, vector<double>* y1);
void runge_kutt_4_temp(vector<double> x, vector<double>* y, vector<double>* y1);
void runge_kutt_4_temp2(vector<double> x, vector<double>* y, vector<double>* y1);
void ballistic(vector<double> x, vector<double>* y);
void newton(vector<double> x, vector<double>* y);
void chordes(vector<double> x, vector<double>* y);
vector<double> gauss_lent(vector<vector<double>> a, vector<double> y);
void finite_subs(vector<double> x, vector<double>* y);
void decompose_lent3(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u);
vector<double> lu_for_lent(vector <vector<double>> a, vector<double> b);
double f2_nonlinear(double x, double y, double y1);
void runge_kutt_4_nonlinear(vector<double> x, vector<double>* y, vector<double>* y1);
void solveMatrix(int n, vector<double> a, vector<double> c, vector<double> b, vector<double> f, vector<double> x);
vector<double> lu(vector <vector<double>> a, vector<double> b);
void decompose(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u);