#include "Header.h"
int main()
{
	int n = 50;
	//≈сли решение не сошлось за 100 итераций при начальном приближении график рисует (0,0);
	double x_start = 0.1, y_start = 0.5, delta = 0.02;
	x_start -= delta * n / 2;
	y_start -= delta * n / 2;
	vector<double> x0(n), y0(n);
	for (int i = 0; i < n; ++i)
	{
		x0[i] = x_start + delta * i;
		y0[i] = y_start + delta * i;
	}
	vector<double> temp(3,0);
	vector<vector<double>> ans_x(n, vector<double>(n)), ans_y(n, vector<double>(n)), num_it(n,vector<double>(n));
	for(int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			temp = simple_it(x0[i], y0[j]);
			ans_x[i][j] = temp[0];
			ans_y[i][j] = temp[1];
			num_it[i][j] = temp[2];

		}
	ofstream fout; 
	fout.open("iter_num.txt"); 
	for (int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
		{
		fout << num_it[i][j] << ' ';
		}
	fout.close();
	fout.open("x0.txt"); 	
	for (int i = 0; i < n; ++i) {
		fout << x0[i] << ' ';
	}
	fout.close();
	fout.open("y0.txt"); 
	for (int i = 0; i < n; ++i) {
		fout << y0[i] << ' ';
	}
	fout.close();
	fout.open("ans_x.txt");
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
	{
		fout << ans_x[i][j] << ' ';
	}
	fout.close();
	fout.open("ans_y.txt");
	for (int i = 0; i < n; ++i) 
	for(int j = 0; j < n; ++j){
		fout << ans_y[i][j] << ' ';
	}
	fout.close();
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			temp = newton(x0[i], y0[j]);
			ans_x[i][j] = temp[0];
			ans_y[i][j] = temp[1];
			num_it[i][j] = temp[2];
		}
	//ofstream fout;
	fout.open("iter_num1.txt");
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			fout << num_it[i][j] << ' ';
		}
	fout.close();
	fout.open("x01.txt");
	for (int i = 0; i < n; i++) {
		fout << x0[i] << ' ';
	}
	fout.close();
	fout.open("y01.txt");
	for (int i = 0; i < n; i++) {
		fout << y0[i] << ' ';
	}
	fout.close();
	fout.open("ans_x1.txt");
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			fout << ans_x[i][j] << ' ';
		}
	fout.close();
	fout.open("ans_y1.txt");
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			fout << ans_y[i][j] << ' ';
		}
	fout.close();
	system("python viz.py");
	return 0;
}
