#include "Header.h"
int main()
{
	double a = 0.0, b = 1.0;
	int n = 100;
	double q0 = -0.5, q1 = 1.0;
	double h = (b - a) / (n - 1);
	vector<double> x(n);
	for (int i = 0; i < n; ++i)
		x[i] = a + i * h;
	vector<double> y(n), y_1(n), y_2(n), y_3(n,0);
	vector<double> y1(n), y1_1(n), y1_2(n), y1_3(n);
	y[0] = q0, y[n - 1] = q1;
	y_1[0] = q0, y_1[n - 1] = q1;
	y_2[0] = q0, y_2[n - 1] = q1;
	y_3[0] = q0, y_3[n - 1] = q1;
	ballistic(x, &y);
	newton(x, &y_1);
	chordes(x, &y_2);
	finite_subs(x, &y_3);
	int n2 = 2 * n;
	double h2 = (b - a) / (n2 - 1);
	vector<double> x2(n2);
	for (int i = 0; i < n2; ++i)
		x2[i] = a + i * h2;
	vector<double> y2(n2), y2_1(n2), y2_2(n2), y2_3(n2);
	vector<double> y21(n2), y21_1(n2), y21_2(n2), y21_3(n2);
	y2[0] = q0, y2[n2 - 1] = q1;
	y2_1[0] = q0, y2_1[n2 - 1] = q1;
	y2_2[0] = q0, y2_2[n2 - 1] = q1;
	y2_3[0] = q0, y2_3[n2 - 1] = q1;
	ballistic(x2, &y2);
	newton(x2, &y2_1);
	chordes(x2, &y2_2);
	finite_subs(x2, &y2_3);
	int n3 = n * 0.5;
	double h3 = (b - a) / (n3 - 1);
	vector<double> x3(n3);
	for (int i = 0; i < n3; ++i)
		x3[i] = a + i * h3;
	vector<double> y3(n3), y3_1(n3), y3_2(n3), y3_3(n3);
	vector<double> y31(n3), y31_1(n3), y31_2(n3), y31_3(n3);
	y3[0] = q0, y3[n3 - 1] = q1;
	y3_1[0] = q0, y3_1[n3 - 1] = q1;
	y3_2[0] = q0, y3_2[n3 - 1] = q1;
	y3_3[0] = q0, y3_3[n3 - 1] = q1;
	ballistic(x3, &y3);
	newton(x3, &y3_1);
	chordes(x3, &y3_2);
	finite_subs(x3, &y3_3);

	cout << "absolute error:\t\t" << "l2:\t\t" << "l1:\t\t" << "cheb:" << endl;//выровнять столбцы (одинаковая длина)
	cout << "ballistic:\t\t" << l2(sub(y, y2)) << ", " << l2(sub(y, y3)) << '\t' << l1(sub(y, y2)) << ", " << l1(sub(y, y3)) << '\t' << cheb(sub(y, y2)) << ", " << cheb(sub(y, y3)) << endl;
	cout << "finite:\t\t\t" << l2(sub(y_3, y2_3)) << ", " << l2(sub(y_3, y3_3)) << '\t' << l1(sub(y_3, y2_3)) << ", " << l1(sub(y_3, y3_3)) << '\t' << cheb(sub(y_3, y2_3)) << ", " << cheb(sub(y_3, y3_3)) << endl;
	cout << "newton:\t\t\t" << l2(sub(y_1, y2_1)) << ", " << l2(sub(y_1, y3_1)) << '\t' << l1(sub(y_1, y2_1)) << ", " << l1(sub(y_1, y3_1)) << '\t' << cheb(sub(y_1, y2_1)) << ", " << cheb(sub(y_1, y3_1)) << endl;
	cout << "chordes:\t\t" << l2(sub(y_2, y2_2)) << ", " << l2(sub(y_2, y3_2)) << '\t' << l1(sub(y_2, y2_2)) << ", " << l1(sub(y_2, y3_2)) << '\t' << cheb(sub(y_2, y2_2)) << ", " << cheb(sub(y_2, y3_2)) << endl;

	cout << "relative error:" << endl;
	cout << "ballistic:\t" << '\t' << l2(sub(y, y2)) / l2(y) << ", " << l2(sub(y, y3)) / l2(y) << '\t' << l1(sub(y, y2)) / l1(y) << ", " << l1(sub(y, y3)) / l1(y) << '\t' << cheb(sub(y, y2)) / cheb(y) << ", " << cheb(sub(y, y3)) / cheb(y) << endl;
	cout << "finite:\t\t" << '\t' << l2(sub(y_3, y2_3)) / l2(y_3) << ", " << l2(sub(y_3, y3_3)) / l2(y_3) << '\t' << l1(sub(y_3, y2_3)) / l1(y_3) << ", " << l1(sub(y_3, y3_3)) / l1(y_3) << '\t' << cheb(sub(y_3, y2_3)) / cheb(y_3) << ", " << cheb(sub(y_3, y3_3)) / cheb(y_3) << endl;
	cout << "newton:\t\t" << '\t' << l2(sub(y_1, y2_1)) / l2(y_1) << ", " << l2(sub(y_1, y3_1)) / l2(y_1) << '\t' << l1(sub(y_1, y2_1)) / l1(y_1) << ", " << l1(sub(y_1, y3_1)) / l1(y_1) << '\t' << cheb(sub(y_1, y2_1)) / cheb(y_1) << ", " << cheb(sub(y_1, y3_1)) / cheb(y_1) << endl;
	cout << "chordes:\t" << '\t' << l2(sub(y_2, y2_2)) / l2(y_2) << ", " << l2(sub(y_2, y3_2)) / l2(y_2) << '\t' << l1(sub(y_2, y2_2)) / l1(y_2) << ", " << l1(sub(y_2, y3_2)) / l1(y_2) << '\t' << cheb(sub(y_2, y2_2)) / cheb(y_2) << ", " << cheb(sub(y_2, y3_2)) / cheb(y_2) << endl;

	ofstream fout("linear.txt");
	for (size_t i = 0; i < x.size(); i++)
		fout << x[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < y.size(); i++)
		fout << y[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < y_3.size(); i++)
		fout << y_3[i] << ' ';
	fout << endl;
	fout.close();
	ofstream fout1("nonlinear.txt");
	for (size_t i = 0; i < x.size(); i++)
		fout1 << x[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < y_1.size(); i++)
		fout1 << y_1[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < y_2.size(); i++)
		fout1 << y_2[i] << ' ';
	fout1 << endl;
	fout1.close();
	system("python \"vis.py\" ");
	system("python \"vis2.py\" ");
	return 0;
}