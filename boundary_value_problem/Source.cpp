#include "Header.h"

double f1(double y1)
{
	return y1;
}
double f2(double x, double y, double y1)
{
	return y1 * x * x + y * tan(x) + (x - 1) * (x - 1);
}
double f2_nonlinear(double x, double y, double y1)
{
	return y1 * y1 * x * x + y * y * tan(x) + (x - 1) * (x - 1);
}

double f2_temp(double x, double y, double y1)
{
	return y1 * x * x + y * tan(x);
}
double f2_temp2(double x, double y, double y1)
{
	return f2_dy(x,y,y1) * y + f2_dy1(x,y,y1);
}
double p(double x)
{
	return x * x;
}
double q(double x)
{
	return tan(x);
}
double r(double x)
{
	return (x - 1) * (x - 1);
}
double f2_dy(double x, double y, double y1)
{
	return 2 * y * tan(x);
}
double f2_dy1(double x, double y, double y1)
{
	return 2 * y1 * x * x;
}
void ballistic(vector<double> x, vector<double>* y)
{
	int n = x.size();
	vector<double> coshi1(n), coshi2(n);
	vector<double> temp_y1(n);
	temp_y1[0] = 0; coshi1[0] = (*y)[0];
	runge_kutt_4(x, &coshi1, &temp_y1);
	temp_y1[0] = 1; coshi2[0] = 0;
	runge_kutt_4_temp(x, &coshi2, &temp_y1);
	for (int i = 1; i < n - 1; ++i)
		(*y)[i] = coshi1[i] + ((*y)[n - 1] - coshi1[n - 1]) / coshi2[n - 1] * coshi2[i];
}

void finite_subs(vector<double> x, vector<double>* y)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	//c1 * y(a) + c2 * y'(a) = c;
	//d1 * y(b) + d2 * y'(b) = d
	double c1 = 1, c2 = 0, c = (*y)[0];
	double d1 = 1, d2 = 0, d = (*y)[n - 1];
	vector<double> gf(n);
	vector <vector<double> > W(n);
	for (int i = 0; i < n; i++) W[i] = vector<double>(n, 0);

	for (int i = 0; i < n; i++)
		W[i][i] = 4e0 + h * h * 2 * q(x[i]);
	for (int i = 1; i < n; i++)
		W[i][i - 1] = -(2e0 + h * p(x[i]));;
	for (int i = 1; i < n; i++)
		W[i - 1][i] = -(2e0 - h * p(x[i]));

	W[0][0] = c1 * h - c2; W[n - 1][n - 1] = d1 * h + d2;
	W[0][1] = c2; W[n - 1][n - 2] = -d2;
	gf[0] = c * h;
	for (int i = 2; i < n - 1; ++i)
	{
		gf[i] = - 2 * r(x[i]) * h * h;
	}

	gf[n - 1] = d * h;

	(*y) = lu(W, gf);
}
void newton(vector<double> x, vector<double>* y)
{
	int n = x.size();
	double s = 0, delta;
	double q1 = (*y)[n - 1];
	double eps = 1e-4;
	int num_it = 0;
	vector<double> y1(n),dyds(n),duds(n);
	while (1) {
		num_it++;
		if (num_it > 500)
			break;
		y1[0] = s;
		runge_kutt_4_nonlinear(x, y, &y1);
		dyds[0] = 0;
		duds[0] = 1;
		runge_kutt_4_temp2(x, &dyds, &duds);
		delta = ((*y)[n - 1] - q1) / dyds[n - 1];
		if (fabs(delta) < eps)
			break;
		s -= delta;
	}
	cout << "newton_num_it=" << num_it << endl;
}
void chordes(vector<double> x, vector<double>* y)
{
	int n = x.size();
	double alpha1 = 0.0, alpha2 = 0.5, delta;
	vector<double> y1(n), y2(n);
	vector<double> temp_y(n);
	double q1 = (*y)[n - 1];
	double eps = 1e-4, temp;
	temp_y[0] = (*y)[0];
	int num_it = 0;
	while (1)
	{
		num_it++;
		if (num_it > 500)
			break;
		y1[0] = alpha1;
		y2[0] = alpha2;
		runge_kutt_4_nonlinear(x, y, &y2);
		runge_kutt_4_nonlinear(x, &temp_y, &y1);
		if (fabs(((*y)[n - 1] - q1) / q1) < eps)
		{
			break;
		}
		temp = alpha2;
		alpha2 -= (alpha2 - alpha1) * ((*y)[n - 1] - q1) / ((*y)[n - 1] - temp_y[n - 1]);
		alpha1 = temp;
	}
	cout << "secant_num_it=" << num_it << endl;
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
void runge_kutt_4_temp(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h / 2.0, h6 = h / 6.0;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (int i = 1; i < n; ++i)
	{
		k1 = f1((*y1)[i - 1]);
		k11 = f2_temp(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
		k2 = f1((*y1)[i - 1] + h2 * k11);
		k21 = f2_temp(x[i - 1] + h2, (*y)[i - 1] + h2 * k1, (*y1)[i - 1] + h2 * k11);
		k3 = f1((*y1)[i - 1] + h2 * k21);
		k31 = f2_temp(x[i - 1] + h2, (*y)[i - 1] + h2 * k2, (*y1)[i - 1] + h2 * k21);
		k4 = f1((*y1)[i - 1] + h * k31);
		k41 = f2_temp(x[i - 1] + h, (*y)[i - 1] + h * k3, (*y1)[i - 1] + h * k31);
		(*y)[i] = (*y)[i - 1] + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		(*y1)[i] = (*y1)[i - 1] + h6 * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
	}
}
void runge_kutt_4_nonlinear(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h / 2.0, h6 = h / 6.0;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (int i = 1; i < n; ++i)
	{
		k1 = f1((*y1)[i - 1]);
		k11 = f2_nonlinear(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
		k2 = f1((*y1)[i - 1] + h2 * k11);
		k21 = f2_nonlinear(x[i - 1] + h2, (*y)[i - 1] + h2 * k1, (*y1)[i - 1] + h2 * k11);
		k3 = f1((*y1)[i - 1] + h2 * k21);
		k31 = f2_nonlinear(x[i - 1] + h2, (*y)[i - 1] + h2 * k2, (*y1)[i - 1] + h2 * k21);
		k4 = f1((*y1)[i - 1] + h * k31);
		k41 = f2_nonlinear(x[i - 1] + h, (*y)[i - 1] + h * k3, (*y1)[i - 1] + h * k31);
		(*y)[i] = (*y)[i - 1] + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		(*y1)[i] = (*y1)[i - 1] + h6 * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
	}
}
void runge_kutt_4_temp2(vector<double> x, vector<double>* y, vector<double>* y1)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double h2 = h / 2.0, h6 = h / 6.0;
	double k1, k2, k3, k4;
	double k11, k21, k31, k41;
	for (int i = 1; i < n; ++i)
	{
		k1 = f1((*y1)[i - 1]);
		k11 = f2_temp2(x[i - 1], (*y)[i - 1], (*y1)[i - 1]);
		k2 = f1((*y1)[i - 1] + h2 * k11);
		k21 = f2_temp2(x[i - 1] + h2, (*y)[i - 1] + h2 * k1, (*y1)[i - 1] + h2 * k11);
		k3 = f1((*y1)[i - 1] + h2 * k21);
		k31 = f2_temp2(x[i - 1] + h2, (*y)[i - 1] + h2 * k2, (*y1)[i - 1] + h2 * k21);
		k4 = f1((*y1)[i - 1] + h * k31);
		k41 = f2_temp2(x[i - 1] + h, (*y)[i - 1] + h * k3, (*y1)[i - 1] + h * k31);
		(*y)[i] = (*y)[i - 1] + h6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		(*y1)[i] = (*y1)[i - 1] + h6 * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
	}
}
void decompose_lent3(vector <vector <double>> a, vector <vector <double>>* l, vector <vector <double>>* u)
{
	*u = a;
	size_t n = a.size();
	size_t p = a[0].size();
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < p; ++j)
			(*l)[i][j] = (*u)[i][j] / (*u)[i][0];


	for (int k = 1; k < n; ++k)
	{
		for (int i = k - 1; i < n; ++i)
			for (int j = 0; j < p; ++j)
				(*l)[i][j] = (*u)[i][j] / (*u)[i][0];

		for (int i = 1; i < p; ++i)
			for (int j = 1; j < p; ++j)
				if (i + k - 1 < n)
					if (j - i >= 0)
						(*u)[i + k - 1][j - i] -= (*l)[k - 1][i] * (*u)[k - 1][j];
		//break;
	}
	//proisv(*l, *u);

}
vector<double> lu_for_lent(vector <vector<double>> a, vector<double> b)
{
	size_t n = a.size();
	size_t p = a[0].size();
	vector<vector<double>> l(n, vector<double>(p)), u(n, vector<double>(p));
	decompose_lent3(a, &l, &u);

	double sum = 0;
	vector<double> x(n, 0), y(n, 0);
	for (int i = 0; i < n; ++i)
	{
		sum = 0;
		for (int j = 0; j < p; ++j)
			if (i - j >= 0)
				sum += l[i - j][j] * y[i - j];
		y[i] = b[i] - sum;

	}
	for (int i = n - 1; i >= 0; --i)
	{
		sum = 0;
		for (int j = 0; j < p; ++j)
			if (i + j < n)
				sum += u[i][j] * x[i + j];
		x[i] = (y[i] - sum) / u[i][0];
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
vector<double> gauss_lent(vector<vector<double>> a, vector<double> y)
{
	size_t n = a.size();
	size_t p = (a[0].size() - 1) / 2;
	//vector<vector<double>> a(n, vector<double>(2 * p - 1));
	/*for (int i = 0; i < n; ++i)
		for (int j = 0; j < p; ++j)
		{

			a[i][j + p - 1] = a_sym[i][j];
			a[i][j] = a_sym[i][p - j - 1];
		}*/
	double b = 0;
	//cout << "a_old=" << endl;
	//show(a);
	vector<double> x(n, 0);
	for (int j = 0; j < n; ++j)
	{
		for (int i = 0; i < p - 1; ++i)

		{
			if ((i != p - 1) & (j + p - 1 - i < n))
			{
				b = a[j][i] / a[j][p - 1];
				//cout << "b2=" << b << endl;
				for (int k = 0; k < p; ++k)
				{
					if ((i + k <= p - 1) & ((j + k) < n))
						a[j + k][i + k] -= b * a[j][k + p - 1];
					else if ((i + k > p - 1)& ((j + p - 1 - i) < n))
						a[j + p - 1 - i][i + k] -= b * a[j][k + p - 1];
				}
				y[j + p - 1 - i] -= b * y[j];
			}

		}
		//if (j ==2)break;
	}
	for (int i = n - 1; i >= 0; --i)
	{
		double sum = 0;
		for (int j = 0; j < p; ++j)
		{
			if (i + j < n)
				sum += x[i + j] * a[i][p - 1 + j];
		}
		x[i] = (y[i] - sum) / a[i][p - 1];
	}
	return x;
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