#include "Header.h"

double f1(double y1)
{
	return y1;
}
double f2(double x, double y, double y1)
{
	return y1 * x * x + y * tan(x) + (x - 1) * (x - 1);
}
void euler(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	for (int i = 1; i < n; ++i)
	{
		(*y)[i] = (*y)[i - 1] + h * f1((*y1)[i - 1]);
		(*y1)[i] = (*y1)[i - 1] + h * f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
	}
}
void euler_predict(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double y_pred, y1_pred;
	double h2 = h / 2.0;
	for (int i = 1; i < n; ++i)
	{
		y_pred = (*y)[i - 1] + h * f1((*y1)[i - 1]);
		y1_pred = (*y1)[i - 1] + h * f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
		(*y)[i] = (*y)[i - 1] + h2 * (f1((*y1)[i - 1]) + f1(y1_pred)) ;
		(*y1)[i] = (*y1)[i - 1] + h2 * (f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]) + f2(x[i],y_pred,y1_pred));
	}
}
void runge_kutt_2(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double y_pred, y1_pred;
	double h2 = h / 2.0;
	for (int i = 1; i < n; ++i)
	{
		//predictor
		y_pred = (*y)[i - 1] + h * f1((*y1)[i - 1]);
		y1_pred = (*y1)[i - 1] + h * f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
		//corrector
		(*y)[i] = (*y)[i - 1] + h2 * (f1((*y1)[i - 1]) + f1(y1_pred));
		(*y1)[i] = (*y1)[i - 1] + h2 * (f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]) + f2(x[i], y_pred, y1_pred));
		//result
		(*y)[i] = (*y)[i - 1] + h2 * (f1((*y1)[i - 1]) + f1(y1_pred));
		(*y1)[i] = (*y1)[i - 1] + h2 * (f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]) + f2(x[i], (*y)[i], (*y1)[i]));
	}
}
void runge_kutt_4(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h / 2.0, h6 = h / 6.0;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (int i = 1; i < n; ++i)
	{
		k1 = f1((*y1)[i - 1]);
		k11 = f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
		k2 = f1((*y1)[i - 1] + h2 * k11);
		k21 = f2(x[i - 1] + h2, (*y)[i - 1] + h2 * k1, (*y1)[i - 1] + h2 * k11);
		k3 = f1((*y1)[i - 1] + h2 * k21);
		k31 = f2(x[i - 1] + h2, (*y)[i - 1] + h2 * k2, (*y1)[i - 1] + h2 * k21);
		k4 = f1((*y1)[i - 1] + h * k31);
		k41 = f2(x[i - 1] + h, (*y)[i - 1] + h * k3, (*y1)[i - 1] + h * k31);
		(*y)[i] = (*y)[i - 1] + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		(*y1)[i] = (*y1)[i - 1] + h6 * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
	}
}
void runge_kutt_4_small(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h / 2.0, h6 = h / 6.0;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (int i = 1; i < 3; ++i)
	{
		k1 = f1((*y1)[i - 1]);
		k11 = f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
		k2 = f1((*y1)[i - 1] + h2 * k11);
		k21 = f2(x[i - 1] + h2, (*y)[i - 1] + h2 * k1, (*y1)[i - 1] + h2 * k11);
		k3 = f1((*y1)[i - 1] + h2 * k21);
		k31 = f2(x[i - 1] + h2, (*y)[i - 1] + h2 * k2, (*y1)[i - 1] + h2 * k21);
		k4 = f1((*y1)[i - 1] + h * k31);
		k41 = f2(x[i - 1] + h, (*y)[i - 1] + h * k3, (*y1)[i - 1] + h * k31);
		(*y)[i] = (*y)[i - 1] + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		(*y1)[i] = (*y1)[i - 1] + h6 * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
	}
}
void adams_3(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double k1 = 0, k2 = 0, k3 = 0;
	double k11 = 0, k21 = 0, k31 = 0;
	double delta_y, delta_y1;
	runge_kutt_4_small(x, y, y1);
	for (int i = 3; i < n; ++i)
	{
		k1 = h * f1((*y1)[i - 1]);
		k11 = h * f2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
		k2 = h * f1((*y1)[i - 2]);
		k21 = h * f2(x[i - 2], (*y)[i - 2], (*y1)[i - 2]);
		k3 = h * f1((*y1)[i - 3]);
		k31 = h * f2(x[i - 3], (*y)[i - 3], (*y1)[i - 3]);
		delta_y = (23.0 * k1 - 16.0 * k2 + 5.0 * k3) / 12.0;
		delta_y1 = (23.0 * k11 - 16.0 * k21 + 5.0 * k31) / 12.0;
		(*y)[i] = (*y)[i - 1] + delta_y;
		(*y1)[i] = (*y1)[i - 1] + delta_y1;
	}
}

double l2(vector<double> x)
{
	double sum = 0;
	for (size_t i = 0; i < x.size(); ++i)
		sum += pow(abs(x[i]), 2);
	return sqrt(sum);
}
double l1(vector<double> x)
{
	double sum = 0;
	for (size_t i = 0; i < x.size(); ++i)
		sum += abs(x[i]);
	return sum;
}
double cheb(vector<double> x)
{
	double max = 0;
	for (size_t i = 0; i < x.size(); ++i)
		if (max <= abs(x[i]))
			max = abs(x[i]);
	return max;
}
double scalar_mult(vector<double> a, vector<double> b)
{
	double res = 0.0;
	for (int i = 0; i < a.size(); ++i) res += a[i] * b[i];

	return res;
}
void show(vector<vector<double>> a)
{
	cout << setprecision(2);
	for (size_t i = 0; i < a.size(); ++i)
	{
		for (size_t j = 0; j < a[0].size(); ++j)
			cout << "\t" << a[i][j] << "\t";
		cout << endl;
	}

}
vector<double> sub(vector<double> a, vector<double> b)
{
	int n = min(a.size(), b.size());
	vector<double> c(n);
	if (a.size() > b.size())
		for (size_t i = 0; i < n; ++i)
			c[i] = b[i] - a[2 * i];
	else if (a.size() < b.size())
		for (size_t i = 0; i < n; ++i)
			c[i] = b[2 * i] - a[i];
	else 
		for (size_t i = 0; i < n; ++i)
			c[i] = b[i] - a[i];
	return c;
}