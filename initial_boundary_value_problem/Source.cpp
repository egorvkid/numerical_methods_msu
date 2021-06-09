#include "Header.h"
#define PI 3.14159265 

double f(double x, double t)
{
	return (1 - t) * sin(PI * x);
}
double psi1(double t)
{
	return 1 - t * t;
}
double psi0(double t)
{
	return t * t - 1;
}
double phi(double x)
{
	return -cos(PI * x);
}
vector<vector<double>> explicit_schema(vector<double> x, vector<double> t)
{
	int n_x = x.size(), n_t = t.size();
	double h_x = (x[n_x - 1] - x[0]) / (n_x - 1);
	double h_t = (t[n_t - 1] - t[0]) / (n_t - 1);
	double th = h_t / (h_x * h_x);
	vector<vector<double>> u(n_x, vector<double> (n_t, 0));
	for (int i = 0; i < n_x; ++i)
		u[i][0] = phi(x[i]);
	for (int j = 0; j < n_t - 1; ++j)
	{
		for (int i = 1; i < n_x - 1; ++i)
		{
			u[i][j + 1] = u[i][j] + th * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) + h_t * f(x[i], t[j] +h_t * 0.5);
			//cout << u[i][j] << ' ';
		}
		u[0][j + 1] = psi0(t[j + 1]);
		u[n_x - 1][j + 1] = psi1(t[j + 1]);
	}
	//cout << endl;
	return u;
}
vector<vector<double>> weighted_schema(vector<double> x, vector<double> t, double sigma)
{
	int n_x = x.size(), n_t = t.size();
	double h_x = (x[n_x - 1] - x[0]) / (n_x - 1);
	double h_t = (t[n_t - 1] - t[0]) / (n_t - 1);
	double th = h_t / (h_x * h_x), a = sigma * th, c = -(1 + 2 * sigma * th);

	
	vector<vector<double>> u(n_x, vector<double>(n_t, 0));
	for (int i = 0; i < n_x; ++i)
		u[i][0] = phi(x[i]);
	vector<double> gf(n_x, 0);
	vector<double> temp(n_x, 0);
	int n = n_x;
	vector <vector<double> > W(n);
	for (int i = 0; i < n; i++) W[i] = vector<double>(n, 0);

	for (int i = 0; i < n; i++)
		W[i][i] = c;
	for (int i = 1; i < n; i++)
		W[i][i - 1] = a;
	for (int i = 1; i < n; i++)
		W[i - 1][i] = a;

	W[0][0] = 1; W[n - 1][n - 1] = 1;
	W[0][1] = 0; W[n - 1][n - 2] = 0;
	for (int j = 0; j < n_t - 1; ++j)
	{
		for (int i = 1; i < n_x - 1; ++i)
			gf[i] = -(u[i][j] + (1 - sigma) * th * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) + h_t * f(x[i], t[j]));
		gf[0] = u[0][j];
		gf[n_x - 1] = u[n_x - 1][j];
		temp = lu(W, gf);
		for (int i = 1; i < n_x - 1; ++i)
		{
			u[i][j + 1] =temp[i];
			//cout << u[i] << ' ';
		}
		u[0][j + 1] = psi0(t[j + 1]);
		u[n_x - 1][j + 1] = psi1(t[j + 1]);
	}
	//cout << endl;
	return u;
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
// Метод прогонки
	// прямой ход
vector<double> progonka(vector<vector<double>> a, vector<double> b)
{
	int n = a.size();
	vector<double> alpha(n - 1), betta(n - 1), gamma(n - 1);
	vector<double> y(n, 0);
	gamma[0] = a[0][0];
	betta[0] = b[0] / gamma[0];
	alpha[0] = -a[0][1] / gamma[0];
	for (int i = 1; i < n - 2; i++) {
		gamma[i] = a[i][i] + a[i][i - 1] * alpha[i - 1];
		betta[i] = (b[i] - a[i][i - 1] * betta[i - 1]) / gamma[i];
		alpha[i] = -a[i][i + 1] / gamma[i];
	}
	gamma[n - 2] = a[n - 2][n - 2] + a[n - 2][n - 3] * alpha[n - 3];
	betta[n - 2] = (b[n - 2] - a[n - 2][n - 3] * betta[n - 3]) / gamma[n - 2];

	//Обратный ход
	y[n - 1] = betta[n - 2];
	for (int i = n - 2; i > 0; i--) 
		y[i] = alpha[i - 1] * y[i + 1] + betta[i - 1];
	
	return y;
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
	//cout << "a_new=" << endl;
	//show(a);
	/*for (int i = 0; i < n; ++i) {
		cout << y[i] << endl;
		x[i] = y[i] / a[i][p - 1];
		//cout << "x" << i << "=" << x[i] << " ";
	}*/
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
vector<double> matrix_mult_lent(vector<vector<double>> a, vector<double> b)
{
	size_t n = a.size();
	size_t p = a[0].size();
	vector<double> x(n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			if (i + j < n)
				x[i] += a[i][j] * b[i + j];
			if (i - j >= 0)
				x[i] += a[i - j][j] * b[i - j];

		}
		x[i] -= a[i][0] * b[i];
	}
	return x;
}
vector<double> grad_lent(vector <vector<double>> a, vector<double> b)
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
		mult = matrix_mult_lent(a, p);
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
vector<double> relax_lent(vector <vector<double>> a, vector<double> b, double w)
{
	size_t n = a.size();
	size_t p = a[0].size();
	vector<double> x_new(n);
	vector<double> x_old(n, 0);
	vector<double> diff(n);
	double err = 1000000;
	double eps = 1e-6;
	double sum = 0;
	double max;
	int it = 0;
	for (int i = 0; i < n; i++)
		x_old[i] = b[i] / a[i][0];
	while (err > eps) {
		++it;
		for (int i = 0; i < n; ++i)
		{
			sum = 0;
			x_new[i] = w * b[i] / a[i][0] - x_old[i] * (w - 1);

			for (int j = 1; j < p; ++j) {

				if (i - j >= 0)
					sum += a[i - j][j] * x_new[i - j];
			}
			for (int j = 0; j < p; ++j)
			{
				if (i + j < n)
					sum += a[i][j] * x_old[j + i];
			}

			sum -= a[i][0] * x_old[i];
			x_new[i] -= w * sum / a[i][0];
		}

		for (int i = 0; i < n; ++i)
			diff[i] = x_new[i] - x_old[i];
		err = cheb(diff);
		x_old = x_new;
	}
	cout << "num_it = " << it << endl;;
	return x_old;
}
vector<double> plain_matrix(vector<vector<double>> a)
{
	int n = a.size();
	int m = a[0].size();
	vector<double> res(n * m);
	for (int j = 0; j < m; ++j)
		for (int i = 0; i < n; ++i)
		{
			res[j * n + i] = a[i][j];
		}
			
	return res;
}