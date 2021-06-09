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

double f(double x, double t);
double psi0(double t);
double psi1(double t);
double phi(double x);
double l2(vector<double> x);
double l1(vector<double> x);
double cheb(vector<double> x);
vector<double> grad(vector <vector<double>> a, vector<double> b);
vector<double> matrix_mult(vector<vector<double>> a, vector<double> b);
double scalar_mult(vector<double> a, vector<double> b);
void show(vector<vector<double>> a);
vector<double> sub(vector<double> a, vector<double> b);
vector<vector<double>> explicit_schema(vector<double> x, vector<double> t);
vector<vector<double>> weighted_schema(vector<double> x, vector<double> t, double sigma);
vector<double> grad_lent(vector <vector<double>> a, vector<double> b);
vector<double> matrix_mult_lent(vector<vector<double>> a, vector<double> b);
vector<double> relax_lent(vector <vector<double>> a, vector<double> b, double w);
vector<double> plain_matrix(vector<vector<double>> a);
vector<double> gauss_lent(vector<vector<double>> a, vector<double> y);
void decompose(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u);
vector<double> lu(vector <vector<double>> a, vector<double> b);
vector<double> progonka(vector<vector<double>> a, vector<double> b);
