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
double l2(vector<double> x);
double l1(vector<double> x);
double cheb(vector<double> x);
double scalar_mult(vector<double> a, vector<double> b);
void show(vector<vector<double>> a); 
void euler(vector<double> x, vector<double>* y, vector<double>* y1);
vector<double> sub(vector<double> a, vector<double> b);
void euler_predict(vector<double> x, vector<double>* y, vector<double>* y1);
void runge_kutt_2(vector<double> x, vector<double>* y, vector<double>* y1);
void runge_kutt_4(vector<double> x, vector<double>* y, vector<double>* y1);
void adams_3(vector<double> x, vector<double>* y, vector<double>* y1);
void runge_kutt_4_small(vector<double> x, vector<double>* y, vector<double>* y1);
