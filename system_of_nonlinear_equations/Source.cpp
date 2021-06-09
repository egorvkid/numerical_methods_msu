#include "Header.h"

double f1(double x, double y)
{
	return pow(x, 2) - y;
}
double f2(double x, double y)
{
	return pow(x, 2) + pow(y, 2) - 1;
}
double f2_y(double x)
{
	return sqrt(fabs(1 - x * x));
}
double f1_x(double y)
{
	return sqrt(fabs(y));
}

double f1_deriv_x(double x, double y)
{
	return 2 * x;
}
double f1_deriv_y(double x, double y)
{
	return -1;
}
double f2_deriv_x(double x, double y)
{
	return 2 * x;
}
double f2_deriv_y(double x, double y)
{
	return 2 * y;
}
int sign(double x)
{
	if (x >= 0)
		return 1;
	else if (x < 0)
		return -1;
	else return 0;
}
vector<double> simple_it(double x0, double y0)
{
	int num_it = 0;
	double eps = 1e-3;
	vector<double> err(2);
	vector<double> iter_xy(2);
	vector<double> ans(3);
	iter_xy[0] = x0, iter_xy[1] = y0;
	while (1)
	{ 
		num_it++;
		if (num_it > 100)
		{
			ans[0] = 0;
			ans[1] = 0;
			break;
		}


		ans[0] = sign(iter_xy[0]) * f1_x(iter_xy[1]);
		ans[1] = sign(iter_xy[1]) * f2_y(iter_xy[0]);

		err[0] = iter_xy[0] - ans[0];
		err[1] = iter_xy[1] - ans[1];
		if (l2(err)/l2(ans) < eps)
			break;
		iter_xy[0] = ans[0];
		iter_xy[1] = ans[1];
	}
	ans[2] = num_it;
	return ans;
}
vector<double> newton(double x0, double y0)
{
	vector<vector<double>> w(2, vector<double>(2));
	vector<double> gf(2), delta(2), err(2);
	vector<double> ans(3,0);
	int num_it = 0;
	double eps = 1e-3;
	ans[0] = x0, ans[1] = y0;
	while (1)
	{
		num_it++;
		if (num_it > 100)
			break;
		w[0][0] = f1_deriv_x(ans[0], ans[1]);
		w[0][1] = f1_deriv_y(ans[0], ans[1]);
		w[1][0] = f2_deriv_x(ans[0], ans[1]);
		w[1][1] = f2_deriv_y(ans[0], ans[1]);
		gf[0] = f1(ans[0], ans[1]);
		gf[1] = f2(ans[0], ans[1]);
		delta = lu(w, gf);
		ans[0] -= delta[0];
		ans[1] -= delta[1];
		if (l2(delta) / l2(ans) < eps)
			break;
	}
	ans[2] = num_it;
	//cout << "num_it=" << num_it << endl;
	return ans;
}
void decompose(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u)
{
	*u = a;
	size_t n = a.size();
	for (int i = 0; i < n; ++i)
		for (int j = i; j < n; ++j)
			(*l)[j][i] = (*u)[j][i] / (*u)[i][i];

	for (int k = 1; k < n; ++k)
	{
		for (int i = k - 1; i < n; ++i)
			for (int j = i; j < n; ++j)
				(*l)[j][i] = (*u)[j][i] / (*u)[i][i];

		for (int i = k; i < n; ++i)
			for (int j = k - 1; j < n; ++j)
				(*u)[i][j] -= (*l)[i][k - 1] * (*u)[k - 1][j];
		//break;
	}
	//proisv(*l, *u);

}
vector<double> lu(vector <vector<double>> a, vector<double> b)
{
	size_t n = a.size();
	size_t p = a[0].size();
	vector<vector<double>> l(n, vector<double>(n)), u(n, vector<double>(n));
	decompose(a, &l, &u);
	//cout <<"l="<< endl;
	//show(l);
	//cout <<"u="<< endl;
	//show(u);

	double sum = 0;
	vector<double> x(n, 0), y(n, 0);
	for (int i = 0; i < n; ++i)
	{
		sum = 0;
		for (int j = 0; j < i; ++j)
			sum += l[i][j] * y[j];
		y[i] = b[i] - sum;

	}
	for (int i = n - 1; i >= 0; --i)
	{
		sum = 0;
		for (int j = n - 1; j >= i; --j)
			sum += u[i][j] * x[j];
		x[i] = (y[i] - sum) / u[i][i];
	}
	return x;
}
vector<double> matrix_mult(vector<vector<double>> a, vector<double> b)
{
	size_t n = a.size();
	vector<double> x(n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			x[i] += a[i][j] * b[j];
	return x;
}
vector<double> grad(vector <vector<double>> a, vector<double> b)
{
	int it = 0;
	size_t n = a.size();
	double b_norm = l2(b);
	vector<double> x(n, 0);
	vector<double> r_old(b), r_new(n);
	vector<double> p(r_old), p_new(n);
	vector<double> mult;
	double r_mult;
	double alpha, beta, eps = 0.001, err = 1000.0;
	while (err > eps) {
		++it;
		mult = matrix_mult(a, p);
		r_mult = scalar_mult(r_old, r_old);
		if (it > 1000)break;
		alpha = r_mult / scalar_mult(mult, p);
		for (int i = 0; i < n; ++i)
			x[i] += alpha * p[i];
		for (int i = 0; i < n; ++i)
			r_new[i] = r_old[i] - alpha * mult[i];
		beta = scalar_mult(r_new, r_new) / r_mult;
		for (int i = 0; i < n; ++i)
			p_new[i] = r_new[i] + beta * p[i];
		r_old = r_new;
		p = p_new;
		err = l2(r_old) / b_norm;
	}
	return x;
}
vector<double> gauss2(vector<vector<double>> a, vector<double> y)
{
	int n = a.size();
	double b;
	vector<double> x(n, 0);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			if (i != j) {

				b = a[i][j] / a[j][j];
				for (int k = 0; k < n; k++) {
					a[i][k] = a[i][k] - b * a[j][k];
				}
				y[i] -= b * y[j];
			}

		}
		//cout << a[j][j] << endl;
	}
	//cout << "a_res=";
	//show(a);
	for (int i = 0; i < n; i++) {
		//cout << y[i] << endl;
		x[i] = y[i] / a[i][i];
		//cout << "x" << i << "=" << x[i] << " ";
	}
	return x;
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