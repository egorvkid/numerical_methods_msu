
#include "Header.h"
using namespace std;

#pragma warning(disable : 4996) 

int main()
{
	double a = 0.0;
	double b = 3.14159;
	int K = 3;
	int N = 6;
	int L = 10;
	int wide = N - 1;//Полуширина ленты

	int M = K * (N - 1) + 1;
	double h_lp = (b - a) / (double)K;
	double h = (b - a) / (double)(M - 1);


	vector<vector<double>> g(M, vector<double>(K * L, 0));

	int cur_index = 0;
	int start = 0;
	int cur = 0;
	double left = a; double right = a + h_lp;
	vector<double> tmp_rand;
	vector<double> tmp_x(N);
	vector<double> tmp_y(N, 0);
	vector<double> x_rand(K * L);
	vector<double> temp(1);
	vector<double> f_rand(K * L);

	for (int k = 0; k < K; ++k)
	{
		tmp_rand = rand(left, right, L);
		for (int j = 0; j < L; ++j)
			x_rand[k * L + j] = tmp_rand[j];
		for (int j = 0; j < L; ++j)
			f_rand[k * L + j] = f(tmp_rand[j]);
		for (int j = 0; j < N; ++j)
			tmp_x[j] = a + h_lp * k + h * j;

		cur = start;
		for (int i = 0; i < N; ++i)
		{

			tmp_y[i] = 1;
			for (int j = 0; j < L; ++j)
			{
				g[cur_index][cur] = LP(tmp_rand[j], tmp_x, tmp_y);
				++cur;
			}
			tmp_y[i] = 0;
			cur = start;
			++cur_index;
		}
		start += L;
		--cur_index;
		left += h_lp; right += h_lp;
	}



	vector<vector<double>> G(M, vector<double>(M));
	vector<vector<double>> lent(M, vector<double>(wide + 1));
	vector<double> gf(M);

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < M; ++j)
			for (int k = 0; k < K * L; ++k)
				G[i][j] += g[i][k] * g[j][k];
		
		for (int j = 0; j < K * L; ++j)
			gf[i] += g[i][j] * f_rand[j];
	}
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < wide + 1; ++j)
		{
			if (i + j >= M) { lent[i][j] = 0; continue; }
			for (int k = 0; k < K * L; ++k)
				lent[i][j] += g[i][k] * g[j + i][k];
		}

		for (int j = 0; j < K * L; ++j)
			gf[i] += g[i][j] * f_rand[j];
		//cout << gf[i] << endl;
	}

	vector<double> c0(M),c1(M), c2(M), c3(M), c4(M), c11(M),c21(M),c31(M),c41(M);

	double w = 1.1;

	c0 = gauss(G, gf);
	c1 = gauss2(G, gf);
	c11 = gauss_lent(lent, gf);
	c2 = lu(G, gf);
	c21 = lu_for_lent(lent, gf);
	c3 = relax(G, gf, w);
	c31 = relax_lent(lent, gf, w);
	c4 = grad(G, gf);
	c41 = grad_lent(lent, gf);
	
	cout << endl << "internet:" << endl;
	for (size_t i = 0; i < c0.size(); ++i)
		cout << c0[i] << ' ';
	cout << endl<<"gauss:" << endl;
	for (size_t i = 0; i < c1.size(); ++i)
		cout << c1[i] << ' ';
	cout <<endl<<"lu:"<< endl;
	for (size_t i = 0; i < c1.size(); ++i)
		cout << c2[i] << ' ';
	cout <<endl<< "relax:"<<endl;
	for (size_t i = 0; i < c1.size(); ++i)
		cout << c3[i] << ' ';
	cout << endl<< "grad:" << endl;
	for (size_t i = 0; i < c1.size(); ++i)
		cout << c4[i] << ' ';
	cout << endl << "relax_lent:" << endl;
	for (size_t i = 0; i < c31.size(); ++i)
		cout << c31[i] << ' ';
	cout << endl << "grad_lent:" << endl;
	for (size_t i = 0; i < c41.size(); ++i)
		cout << c41[i] << ' ';
	cout << endl << "gauss_lent:" << endl;
	for (size_t i = 0; i < c11.size(); ++i)
		cout << c11[i] << ' ';
	cout << endl << "lu_lent:" << endl;
	for (size_t i = 0; i < c21.size(); ++i)
		cout << c21[i] << ' ';
	cout << endl;
	vector<double> nev1(M), nev2(M), nev3(M), nev4(M), nev5(M), nev6(M), nev7(M), nev8(M), nev9(M);
	vector<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7, tmp8, tmp9;
	tmp1 = matrix_mult(G, c1);
	tmp2 = matrix_mult(G, c2);
	tmp3 = matrix_mult(G, c3);
	tmp4 = matrix_mult(G, c4);
	tmp5 = matrix_mult(G, c0);
	tmp6 = matrix_mult(G, c11);
	tmp7 = matrix_mult(G, c21);
	tmp8 = matrix_mult(G, c31);
	tmp9 = matrix_mult(G, c41);
	for (int i = 0; i < M; ++i)
	{
		nev1[i] = tmp1[i] - gf[i];
		nev2[i] = tmp2[i] - gf[i];
		nev3[i] = tmp3[i] - gf[i];
		nev4[i] = tmp4[i] - gf[i];
		nev5[i] = tmp5[i] - gf[i];
		nev6[i] = tmp6[i] - gf[i];
		nev7[i] = tmp7[i] - gf[i];
		nev8[i] = tmp8[i] - gf[i];
		nev9[i] = tmp9[i] - gf[i];
	}
	double gf_l1 = l1(gf), gf_l2 = l2(gf), gf_cheb = cheb(gf);
	cout << setprecision(5);
	cout << "  \tl1 \t\t" << "  l2 \t\t" << "  cheb \t\t" << endl;
	cout << "gauss:\t"<<l1(nev1) << '\t' << l2(nev1) << '\t' << cheb(nev1) << endl;
	cout << "lowup:\t" << l1(nev2) << '\t' << l2(nev2) << '\t' << cheb(nev2) << endl;
	cout << "relax:\t" << l1(nev3) << '\t' << l2(nev3) << '\t' << cheb(nev3) << endl;
	cout << "sgrad:\t" << l1(nev4) << '\t' << l2(nev4) << '\t' << cheb(nev4) << endl;
	cout << "NotMy:\t" << l1(nev5) << '\t' << l2(nev5) << '\t' << cheb(nev5) << endl;
	cout << "glent:\t" << l1(nev6) << '\t' << l2(nev6) << '\t' << cheb(nev6) << endl;
	cout << "llent:\t" << l1(nev7) << '\t' << l2(nev7) << '\t' << cheb(nev7) << endl;
	cout << "rlent:\t" << l1(nev8) << '\t' << l2(nev8) << '\t' << cheb(nev8) << endl;
	cout << "slent:\t" << l1(nev9) << '\t' << l2(nev9) << '\t' << cheb(nev9) << endl;
	cout << "relative:" << endl;
	cout << "\t l1 \t\t" << "  l2 \t\t" << "  cheb \t" << endl;
	cout << "gauss:\t" << l1(nev1)/gf_l1 << '\t' << l2(nev1)/gf_l2 << '\t' << cheb(nev1)/gf_cheb << endl;
	cout << "lowup:\t" << l1(nev2) / gf_l1 << '\t' << l2(nev2) / gf_l2 << '\t' << cheb(nev2) / gf_cheb << endl;
	cout << "relax:\t" << l1(nev3) / gf_l1 << '\t' << l2(nev3) / gf_l2 << '\t' << cheb(nev3) / gf_cheb << endl;
	cout << "sgrad:\t" << l1(nev4) / gf_l1 << '\t' << l2(nev4) / gf_l2 << '\t' << cheb(nev4) / gf_cheb << endl;
	cout << "NotMy:\t" << l1(nev5) / gf_l1 << '\t' << l2(nev5) / gf_l2 << '\t' << cheb(nev5) / gf_cheb << endl;
	cout << "glent:\t" << l1(nev6) / gf_l1 << '\t' << l2(nev6) / gf_l2 << '\t' << cheb(nev6) / gf_cheb << endl;
	cout << "llent:\t" << l1(nev7) / gf_l1 << '\t' << l2(nev7) / gf_l2 << '\t' << cheb(nev7) / gf_cheb << endl;
	cout << "rlent:\t" << l1(nev8) / gf_l1 << '\t' << l2(nev8) / gf_l2 << '\t' << cheb(nev8) / gf_cheb << endl;
	cout << "slent:\t" << l1(nev9) / gf_l1 << '\t' << l2(nev9) / gf_l2 << '\t' << cheb(nev9) / gf_cheb << endl;

	return 0;

}

