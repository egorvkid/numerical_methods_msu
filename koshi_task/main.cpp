#include "Header.h"
int main()
{
	double a = 0, b = 1;
	int n = 200;
	double q0 = 0, q1 = -1;
	double h = (b - a) / (n - 1);
	vector<double> x(n);
	for (int i = 0; i < n; ++i)
		x[i] = a + i * h;
	vector<double> y(n), y_1(n), y_2(n),y_3(n),y_4(n);
	vector<double> y1(n), y1_1(n), y1_2(n),y1_3(n), y1_4(n);
	y[0] = q0, y1[0] = q1;
	y_1[0] = q0, y1_1[0] = q1;
	y_2[0] = q0, y1_2[0] = q1;
	y_3[0] = q0, y1_3[0] = q1;
	y_4[0] = q0, y1_4[0] = q1;
	euler(x, &y, &y1);
	euler_predict(x, &y_1, &y1_1);
	runge_kutt_2(x, &y_2, &y1_2);
	runge_kutt_4(x, &y_3, &y1_3);
	adams_3(x, &y_4, &y1_4);
	int n2 = 2 * n;
	double h2 = (b - a) / (n2 - 1);
	vector<double> x2(n2);
	for (int i = 0; i < n2; ++i)
		x2[i] = a + i * h2;
	// _ euler, _1 euler_predict
	vector<double> y2(n2), y2_1(n2), y2_2(n2),y2_3(n2),y2_4(n2);
	vector<double> y21(n2), y21_1(n2), y21_2(n2), y21_3(n2), y21_4(n2);
	y2[0] = q0, y21[0] = q1;
	y2_1[0] = q0, y21_1[0] = q1;
	y2_2[0] = q0, y21_2[0] = q1;
	y2_3[0] = q0, y21_3[0] = q1;
	y2_4[0] = q0, y21_4[0] = q1;
	euler(x2, &y2, &y21);
	euler_predict(x2, &y2_1, &y21_1);
	runge_kutt_2(x2, &y2_2, &y21_2);
	runge_kutt_4(x2, &y2_3, &y21_3);
	adams_3(x2, &y2_4, &y21_4);
	int n3 = n * 0.5;
	double h3 = (b - a) / (n3 - 1);
	vector<double> x3(n3);
	for (int i = 0; i < n3; ++i)
		x3[i] = a + i * h3;
	vector<double> y3(n3), y3_1(n3), y3_2(n3), y3_3(n3), y3_4(n3);
	vector<double> y31(n3), y31_1(n3), y31_2(n3), y31_3(n3), y31_4(n3);
	y3[0] = q0, y31[0] = q1;
	y3_1[0] = q0, y31_1[0] = q1;
	y3_2[0] = q0, y31_2[0] = q1;
	y3_3[0] = q0, y31_3[0] = q1;
	y3_4[0] = q0, y31_4[0] = q1;
	euler(x3, &y3, &y31);
	euler_predict(x3, &y3_1, &y31_1);
	runge_kutt_2(x3, &y3_2, &y31_2);
	runge_kutt_4(x3, &y3_3, &y31_3);
	adams_3(x3, &y3_4, &y31_4);
	cout << "absolute error:\t\t" << "l2:\t\t" << "l1:\t\t" << "cheb:" << endl;
	cout << "euler:\t\t" << l2(sub(y, y2)) << ", " << l2(sub(y, y3)) << '\t' << l1(sub(y, y2)) << ", " << l1(sub(y, y3)) << '\t' << cheb(sub(y, y2)) << ", " << cheb(sub(y, y3)) << endl;
	cout << "euler_predict:\t" << l2(sub(y_1, y2_1)) << ", " << l2(sub(y_1, y3_1)) << '\t' << l1(sub(y_1, y2_1)) << ", " << l1(sub(y_1, y3_1)) << '\t' << cheb(sub(y_1, y2_1)) << ", " << cheb(sub(y_1, y3_1)) << endl;
	cout << "runge_kutt_2:\t" << l2(sub(y_2, y2_2)) << ", " << l2(sub(y_2, y3_2)) << '\t' << l1(sub(y_2, y2_2)) << ", " << l1(sub(y_2, y3_2)) << '\t' << cheb(sub(y_2, y2_2)) << ", " << cheb(sub(y_2, y3_2)) << endl;
	cout << "runge_kutt_4:\t" << l2(sub(y_3, y2_3)) << ", " << l2(sub(y_3, y3_3)) << '\t' << l1(sub(y_3, y2_3)) << ", " << l1(sub(y_3, y3_3)) << '\t' << cheb(sub(y_3, y2_3)) << ", " << cheb(sub(y_3, y3_3)) << endl;
	cout << "adams_3:\t" << l2(sub(y_4, y2_4)) << ", " << l2(sub(y_4, y3_4)) << '\t' << l1(sub(y_4, y2_4)) << ", " << l1(sub(y_4, y3_4)) << '\t' << cheb(sub(y_4, y2_4)) << ", " << cheb(sub(y_4, y3_4)) << endl;
	cout << "relative error:" << endl;
	cout << "euler:\t" << '\t' << l2(sub(y, y2)) / l2(y) << ", " << l2(sub(y, y3)) / l2(y) << '\t' << l1(sub(y, y2)) / l1(y) << ", " << l1(sub(y, y3)) / l1(y) << '\t' << cheb(sub(y, y2)) / cheb(y) << ", " << cheb(sub(y, y3)) / cheb(y) << endl;
	cout << "euler_predict:" << '\t' << l2(sub(y_1, y2_1)) / l2(y_1) << ", " << l2(sub(y_1, y3_1)) / l2(y_1) << '\t' << l1(sub(y_1, y2_1)) / l1(y_1) << ", " << l1(sub(y_1, y3_1)) / l1(y_1) << '\t' << cheb(sub(y_1, y2_1)) / cheb(y_1) << ", " << cheb(sub(y_1, y3_1)) / cheb(y_1) << endl;
	cout << "runge_kutt_2:" << '\t' << l2(sub(y_2, y2_2)) / l2(y_2) << ", " << l2(sub(y_2, y3_2)) / l2(y_2) << '\t' << l1(sub(y_2, y2_2)) / l1(y_2) << ", " << l1(sub(y_2, y3_2)) / l1(y_2) << '\t' << cheb(sub(y_2, y2_2)) / cheb(y_2) << ", " << cheb(sub(y_2, y3_2)) / cheb(y_2) << endl;
	cout << "runge_kutt_4:" << '\t' << l2(sub(y_3, y2_3)) / l2(y_3) << ", " << l2(sub(y_3, y3_3)) / l2(y_3) << '\t' << l1(sub(y_3, y2_3)) / l1(y_3) << ", " << l1(sub(y_3, y3_3)) / l1(y_3) << '\t' << cheb(sub(y_3, y2_3)) / cheb(y_3) << ", " << cheb(sub(y_3, y3_3)) / cheb(y_3) << endl;
	cout << "adams_3:" << '\t' << l2(sub(y_4, y2_4)) / l2(y_4) << ", " << l2(sub(y_4, y3_4)) / l2(y_4) << '\t' << l1(sub(y_4, y2_4)) / l1(y_4) << ", " << l1(sub(y_4, y3_4)) / l1(y_4) << '\t' << cheb(sub(y_4, y2_4)) / cheb(y_4) << ", " << cheb(sub(y_4, y3_4)) / cheb(y_4) << endl;
	ofstream fout("y.txt");
	for (size_t i = 0; i < x.size(); i++)
		fout << x[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < y.size(); i++)
		fout << y[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < y_1.size(); i++)
		fout << y_1[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < y_2.size(); i++)
		fout << y_2[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < y_3.size(); i++)
		fout << y_3[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < y_3.size(); i++)
		fout << y_3[i] << ' ';
	fout << endl;
	for (size_t i = 0; i < y_4.size(); i++)
		fout << y_4[i] << ' ';
	fout << endl;

	fout.close();
	ofstream fout1("y1.txt");
	for (size_t i = 0; i < x.size(); i++)
		fout1 << x[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < y1.size(); i++)
		fout1 << y1[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < y1_1.size(); i++)
		fout1 << y1_1[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < y1_2.size(); i++)
		fout1 << y1_2[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < y1_3.size(); i++)
		fout1 << y1_3[i] << ' ';
	fout1 << endl;
	for (size_t i = 0; i < y1_4.size(); i++)
		fout1 << y1_4[i] << ' ';
	fout1 << endl;
	fout1.close();
	system("python \"vis.py\" ");
	system("python \"vis2.py\" ");
	return 0;
}