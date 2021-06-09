#include "Header.h"


int main()
{
	double a = 0.0, b = 1.0;
	double sigma = 0.5;
	int n_x0 = 10, n_t = 16 * n_x0 * n_x0;
	int  n_x01 = 10, n_t2 = 200; //* (1 - 2 * sigma);
	double q0 = 0, q1 = -1;
	double h_x = (b - a) / (n_x0 - 1);
	double h_x1 = (b - a) / (n_x01 - 1);
	double h_t = (b - a) / (n_t - 1);
	double h_t2 = (b - a) / (n_t2 - 1);
	

	vector<double> x0(n_x0), x01(n_x01), t(n_t), t2(n_t2);
	for (int i = 0; i < n_x0; ++i)
		x0[i] = a + i * h_x;
	for (int i = 0; i < n_x01; ++i)
		x01[i] = a + i * h_x1;
	for (int i = 0; i < n_t; ++i)
		t[i] = a + i * h_t;
	for (int i = 0; i < n_t2; ++i)
		t2[i] = a + i * h_t2;
	vector<vector<double>> u_1_h1,u_2_h1;

	u_1_h1 = explicit_schema(x0, t);
	u_2_h1 = weighted_schema(x01, t2, sigma);

	int n_x1 = 2 * n_x0;
	int n_x11 = 2 * n_x01;
	h_x = (b - a) / (n_x1 - 1);
	h_x1 = (b - a) / (n_x11 - 1);
	vector<double> x1(n_x1);
	for (int i = 0; i < n_x1; ++i)
		x1[i] = a + i * h_x;
	vector<double> x11(n_x11);
	for (int i = 0; i < n_x1; ++i)
		x11[i] = a + i * h_x1;
	vector<vector<double>> u_1_h2, u_2_h2;
	u_1_h2 = explicit_schema(x1, t);
	u_2_h2 = weighted_schema(x11, t2, sigma);

	int n_x2 = 0.5 * n_x0;
	h_x = (b - a) / (n_x2 - 1);
	vector<double> x2(n_x2);
	for (int i = 0; i < n_x2; ++i)
		x2[i] = a + i * h_x;
	int n_x21 = 0.5 * n_x01;
	h_x1 = (b - a) / (n_x21 - 1);
	vector<double> x21(n_x21);
	for (int i = 0; i < n_x21; ++i)
		x21[i] = a + i * h_x1;
	vector<vector<double>> u_1_h3, u_2_h3;
	u_1_h3 = explicit_schema(x2, t);
	u_2_h3 = weighted_schema(x21, t2, sigma);
	vector<double> err1h1, err1h2, err2h1, err2h2;
	err1h1 = sub(plain_matrix(u_1_h1), plain_matrix(u_1_h2));
	err1h2 = sub(plain_matrix(u_1_h1), plain_matrix(u_1_h3));
	err2h1 = sub(plain_matrix(u_2_h1), plain_matrix(u_2_h2));
	err2h2 = sub(plain_matrix(u_2_h1), plain_matrix(u_2_h3));

	cout << "absolute error:\t\t\t" << "l2:\t\t" << "l1:\t\t" << "cheb:" << endl;//выровнять столбцы (одинаковая длина)
	cout << "explicit_schema:\t\t" << l2(err1h1) << ", " << l2(err1h2) << '\t' << l1(err1h1) << ", " << l1(err1h2) << '\t' << cheb(err1h1) << ", " << cheb(err1h2) << endl;
	cout << "weighted_schema:\t\t" << l2(err2h1) << ", " << l2(err2h2) << '\t' << l1(err2h1) << ", " << l1(err2h2) << '\t' << cheb(err2h1) << ", " << cheb(err2h2) << endl;
	cout << "relative error:" << endl;
	cout << "explicit_schema:\t\t" << l2(err1h1) / l2(plain_matrix(u_1_h1)) << ", " << l2(err1h2) / l2(plain_matrix(u_1_h1)) << '\t' << l1(err1h1) / l1(plain_matrix(u_1_h1)) << ", " << l1(err1h2) / l1(plain_matrix(u_1_h1)) << '\t' << cheb(err1h1) / cheb(plain_matrix(u_1_h1)) << ", " << cheb(err1h2) / cheb(plain_matrix(u_1_h1)) << endl;
	cout << "weighted_schema:\t\t" << l2(err2h1) / l2(plain_matrix(u_2_h1)) << ", " << l2(err2h2) / l2(plain_matrix(u_2_h1)) << '\t' << l1(err2h1) / l1(plain_matrix(u_2_h1)) << ", " << l1(err2h2) / l1(plain_matrix(u_2_h1)) << '\t' << cheb(err2h1) / cheb(plain_matrix(u_2_h1)) << ", " << cheb(err2h2) / cheb(plain_matrix(u_2_h1)) << endl;
	ofstream fout("explicit.txt");
	for (size_t i = 0; i < x0.size(); i++)
		fout << x0[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < t.size(); i++)
		fout << t[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < u_1_h1.size(); i++)
	{
		for (size_t j = 0; j < u_1_h1[0].size(); j++)
		{
			fout << u_1_h1[i][j] << ' ';
		}
		fout << endl;
	}
	fout.close();
	ofstream fout1("weighted.txt");
	for (size_t i = 0; i < x01.size(); i++)
		fout1 << x01[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < t2.size(); i++)
		fout1 << t2[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < u_2_h1.size(); i++)
	{
		for (size_t j = 0; j < u_2_h1[0].size(); j++)
			fout1 << u_2_h1[i][j] << ' ';
		fout1 << endl;
	}
	fout1.close();
	system("python \"vis.py\" \"explicit.txt\"");
	system("python \"vis.py\" \"weighted.txt\"");
	return 0;
}