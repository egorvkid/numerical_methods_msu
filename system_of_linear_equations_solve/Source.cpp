#include "Header.h"
#include <random>
#include <cmath>


double f(double x)                            // определяет интерполируемую функцию
{
	return sin(x);
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
	vector<double> x(n, 0), start(n);
	vector<double> ax = matrix_mult_lent(a, x);
	for (int i = 0; i < n; ++i)
		start[i] = b[i] - ax[i];
	vector<double> r(start);//антиградиент
	vector<double> z(r);//базисный вектор
	vector<double> mult;
	double r_mult_new, r_mult_old;
	double alpha, beta, eps = 1e-6, err = 1000.0;
	r_mult_old = scalar_mult(r, r);
	while (err > eps) {
		++it;
		mult = matrix_mult_lent(a, z);
		if (it > 1000)break;
		alpha = r_mult_old / scalar_mult(mult, z);
		for (int i = 0; i < n; ++i)
			x[i] += alpha * z[i];
		for (int i = 0; i < n; ++i)
			r[i] -= alpha * mult[i];
		r_mult_new = scalar_mult(r, r);
		beta = r_mult_new / r_mult_old;
		for (int i = 0; i < n; ++i)
			z[i] = r[i] + beta * z[i];
		r_mult_old = r_mult_new;
		err = l2(r) / b_norm;
	}
	cout << "grad_lent_num_it = " << it << endl;
	return x;
}
vector<double> grad(vector <vector<double>> a, vector<double> b)
{
	int it = 0;
	size_t n = a.size();
	double b_norm = l2(b);
	vector<double> x(n, 0), start(n);
	vector<double> ax = matrix_mult(a, x);
	for (int i = 0; i < n; ++i)
		start[i] = b[i] - ax[i];
	vector<double> r(start);//антиградиент
	vector<double> z(r);//базисный вектор
	vector<double> mult;
	double r_mult_new, r_mult_old;
	double alpha, beta,eps = 1e-6,err=1000.0;
	r_mult_old = scalar_mult(r, r);
	while (err>eps) {
		++it;
		mult = matrix_mult(a, z);
		if (it > 1000)break;
		alpha =  r_mult_old / scalar_mult(mult, z);
		for (int i = 0; i < n; ++i)
			x[i] += alpha * z[i];
		for (int i = 0; i < n; ++i)
			r[i] -= alpha * mult[i];
		r_mult_new = scalar_mult(r, r);
		beta = r_mult_new / r_mult_old;
		for (int i = 0; i < n; ++i)
			z[i] = r[i] + beta * z[i];
		r_mult_old = r_mult_new;
		err = l2(r) / b_norm;
	}
	cout << "grad_num_it = " << it << endl;
	return x;
}

//подобрать w
vector<double> relax(vector <vector<double>> a, vector<double> b, double w)
{
	size_t n = a.size();
	vector<double> x(n,0);
	vector<double> diff(n);
	double err = 10000000;
	double eps = 1e-6;
	double sum = 0;
	int  it = 0;
	for (int i = 0; i < n; i++)
		x[i] = b[i] / a[i][i];
	while (err>eps) {
		++it;
		for (int i = 0; i < n; ++i)
		{
			sum = 0;
			
			for (int k = 0; k < n; ++k)
				sum += a[i][k] * x[k];
			diff[i] = w * b[i] / a[i][i] - w * sum / a[i][i];
			x[i] += diff[i];
		}

		err = cheb(diff);
	}
	cout << "relax_num_it = " << it << endl;
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
			diff[i] = w * b[i] / a[i][0] - w * sum / a[i][0];
			x_new[i] -= w * sum / a[i][0];
		}

		for (int i = 0; i < n; ++i)
			diff[i] = x_new[i] - x_old[i];
		err = cheb(diff);
		x_old = x_new;
	}
	cout << "relax_lent_num_it = " << it << endl;;
	return x_old;
}

//для ленточных матриц
void decompose(vector <vector <double>> a, vector <vector <double>> *l,vector <vector <double>> *u)
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
			if (i - j >= 0 )
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
vector<double> lu(vector <vector<double>> a, vector<double> b)
{
	size_t n = a.size();
	size_t p = a[0].size();
	vector<vector<double>> l(n, vector<double>(n)) , u(n, vector<double>(n));
	decompose(a, &l, &u);

	double sum = 0;
	vector<double> x(n,0), y(n,0);
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
	}


	for (int i = 0; i < n; i++) 
		x[i] = y[i] / a[i][i];
	
	return x;
	}
//гаусс для ленточных матриц
vector<double> gauss_lent(vector<vector<double>> a_sym, vector<double> y)
{
	size_t n = a_sym.size();
	size_t p = a_sym[0].size();
	vector<vector<double>> a(n, vector<double>(2 * p - 1));
	for(int i = 0; i < n; ++i)
		for (int j = 0; j < p; ++j)
		{

			a[i][j + p - 1] = a_sym[i][j];
			a[i][j] = a_sym[i][p - j - 1];
		}
	double b = 0;

	vector<double> x(n,0);
	for(int j = 0; j < n; ++j)
	{
		for (int i = 0; i < p - 1; ++i)
		
			{
				if ((i != p - 1) & (j + p - 1 - i < n))
				{
					b = a[j][i] / a[j][p - 1];
					for (int k = 0; k < p; ++k)
					{
						if ((i + k <= p - 1) & ((j + k) < n))
							a[j + k][i + k] -= b * a[j][k + p - 1];
						else if ((i + k > p - 1)  & ((j + p - 1 - i) < n))
							a[j + p - 1 - i][i + k] -= b * a[j][k + p - 1];
					}
					y[j + p - 1 - i] -= b * y[j];
				}

			}
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
vector<double> gauss(vector<vector<double>> a, vector<double> y)
{
	int n = y.size();
	vector<double> x(n, 0);
	double	max;
	int k, index;
	double eps = 0.00001;  // точность
	k = 0;
	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "error ";
			cout << index << endl;
			return x;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	return x;
}

double get_scalar_value(vector<double> a, vector<double> b)
{
	double val = 0.0;
	for (int i = 0; i < a.size(); ++i) val += a[i] * b[i];

	return val;
}
double LP(double x, vector<double> xv, vector<double> yv) { //Lagrange polynomial
	int size = xv.size(); //Количество точек (для удобства)
	long double sum = 0; //Значение функции
	//double eps = 1e-1;
	for (int i = 0; i < size; i++) {
		double mul = 1; //Произведение
		for (int j = 0; j < size; j++) {

			if (i != j) {
				//if (abs(xv[i] - xv[j]) < eps) mul*=1;
				mul *= (x - xv[j]) / (xv[i] - xv[j]);
			}
		}
		sum += yv[i] * mul;
	}
	return sum;
}
double Lagr_basis(vector<double> x, double t, int alpha) {
	double mul = 1;
	for (size_t j = 0; j < x.size(); ++j)
		if (j != alpha)
			mul *= (t - x[alpha]) / (x[alpha] - x[j]);
	return mul;
}

void LocLag(vector<long double>* res, double a, double b, vector<double> x, vector<double> y, int n)
{
	//double h = (b - a) / n;
	vector<double> arr(n);
	for (size_t i = 0; i < arr.size(); ++i)
	{
		arr[i] = a + (b - a) / n * i;
		(*res).push_back(LP(arr[i], x, y)); //cout << arr[i] << ' ' << (*res)[i]<< endl;
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
vector<double> rand(double a, double b, int l)
{
	vector<double> tmp(l);
	double EPS = 0.01;
	for (size_t i = 0; i < tmp.size(); ++i)
	{
		tmp[i] = (static_cast <double> (rand()) / static_cast <double> (RAND_MAX))* (b - a) + a;
	}
	//sort(tmp.begin(), tmp.end());
	return tmp;
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
			cout <<"\t"<< a[i][j] <<"\t";
		cout << endl;
	}
	
}
