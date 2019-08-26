#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>
#include"INA.h"
using namespace std;
int cycnum = 1;
vector<bool> change_mark;
vector<double> initial_degree;
vector<double> initial_single_degree;
vector<vector<vector<double>>> Knn_all;
void main()
{	
	int Agen = 5, k = 10, layernum = LayerN;
	vector<int> network_size;
	int N = 100;
	vector<double> average_ini(Agen, 0);
	vector<double> average_res(Agen, 0);
	vector<double> average_R3(Agen, 0);
	vector<double> R_res;
	vector<vector<vector<double>>> Knn_average(Agen);
	ofstream result;
	ofstream Knnresult;
	result.open("result.txt");
	Knnresult.open("Knnresult.txt");
	for (int i = 0; i < Agen; i++)
	{
		for (int j = 0; j < k; j++)
		{
			network_size.push_back(N);
			R_res.push_back(mainloop_MA(N));
			cycnum++;
		}
		N += 100;
	}
	N -= 100;
	for (int j = 0,t=0; j < Agen; j++)
	{
		for (int i = 0; i < k; i++,t++)
		{
			average_ini[j] += initial_degree[t];
			average_res[j] += R_res[t];
		}
	}
	for (int j = 0; j < Agen; j++)
	{
		average_ini[j] = average_ini[j] / k;
		average_res[j] = average_res[j] / k;
	}
	for (int i = 0; i < Agen; i++) {
		Knn_average[i] = Knn_all[i*k];
		for (int j = 0; j < Knn_average[i].size(); j++) {
			for (int y = 0; y < Knn_average[i][y].size(); y++) {
				Knn_average[i][j][y] = 0;
			}
		}
	}
	for (int i = 0, b = 0; i < Agen; i++) {
		for (int j = 0; j < k; j++) {
			for (int y = 0; y < network_size[b]; y++) {
				Knn_average[i][0][y] += Knn_all[b][0][y];
				Knn_average[i][1][y] += Knn_all[b][1][y];
			}
			b++;
		}
		for (int y = 0; y < network_size[i*k]; y++) {
			Knn_average[i][0][y] = Knn_average[i][0][y] / k;
			Knn_average[i][1][y] = Knn_average[i][1][y] / k;
		}
	}

	for (int i = 0; i < Agen*k; i++)
	{
		result << "the cycnum is " << i <<" layer="<<layernum<< ", the size is " << network_size[i] << endl;
		if (change_mark[0] == false )
		{
			result << "the degree has change" << endl;
		}
		result << "the initial degree is " << initial_degree[i] << endl;
		result << "the interdependent network MA--edgecrossover is " << R_res[i] << endl;
	}
	result << "网络大小:" << N - 100 * (Agen - 1) << "-" << N << endl;
	result << "平均初始度" << "\t多层网络优化结果 " << endl;
	for (int j = 0; j < Agen; j++)
	{	
		result  << average_ini[j] << "\t" << average_res[j]  << endl;
	}
	for (int i = 0; i < Agen; i++) {
		Knnresult << "N: " << network_size[i*k] << endl;
		Knnresult << "k		" << "Knn" <<"初始"<<endl;
		for (int j = 0; j < network_size[i*k]; j++) {
			Knnresult << j << "		" << Knn_average[i][0][j] << endl;
		}
		Knnresult << "k		" << "Knn" << "结束" << endl;
		for (int j = 0; j < network_size[i*k]; j++) {
			Knnresult << j << "		" << Knn_average[i][1][j] << endl;
		}
	}
	Knnresult.close();
	result << "MA-SF-2网络" << endl;
	result.close();
}



double mainloop_MA(int N)
{
	srand(cycnum);
	int gen = 10000, falsenum = 0;
	bool searchoperator = true;//true-localsearch,false-randomswap
	int layer = 0, layer_ops = 1,maxnode=0, pop_size = 10;
	double R_max = 0, R_transmit = 0, R_best = 0;
	int begin_layer = 0, end_layer = LayerN;
	srand(cycnum);
	double R2[LayerN], R2_mid[LayerN], R2_comb = 0;
	int countnum = 0;
	vector<vector<double>> Knn(2, vector<double>(N,0));
	vector<double> R(pop_size, 0);
	vector<double> Rc(pop_size, 0);
	vector<vector<int>> choice(pop_size / 2, vector<int>(2, -1));
	vector<vector<int>> choice2(pop_size, vector<int>(2, -1));
	//初始化网络
	vector<vector<node>> rnode(pop_size, vector<node>(N));
	vector<vector<node>> rnodec(pop_size, vector<node>(N));
	vector<vector<node>> rnode_backup( pop_size, vector<node>(N));
	vector<vector<node>> rnodec_backup(pop_size, vector<node>(N));
	vector<vector<node>> rnode_2(pop_size, vector<node>(N));
	vector<node> rnode_best;
	

	for (int i = 0; i < pop_size; i++)
	{
		initialization_degree(rnode[i]);//初始化度
		activate_node(rnode[i]);//激活节点
	}

	for (int i = 0; i < LayerN; i++)
	{
		initialization_SF(rnode[0], i);//创建网络
	}

//	initialization_ER(rnode[0], 0);
//	initialization_SF(rnode[0], 1);

/*	int k_degree = 0; double k_res;
	for (int i = 0; i < rnode[0].size(); i++)
	{
		k_degree += rnode[0][i].neibor[1].size();
	}
	k_res = (double)k_degree / rnode[0].size();*/
	
	initialSA(rnode[0]);//初始化SA
	rnode_2[0] = rnode[0];
	mark_degree(rnode[0]);
	mark_degree(rnode_2[0]);
	R[0]= computeR(rnode[0], 0,LayerN);
	initial_degree.push_back(R[0]);
	for (int i = 0; i < rnode[0].size(); i++) {
		Knn[0][i]=computeKnn(rnode[0], 0, i);
	}
//	initial_single_degree.push_back(computeR_single(rnode[0], layer));
	print_network(rnode[0],0, 0);

	//预训练
/*	for (int i = 0; i < N * 30; i++) {
		for (int j = 0; j < LayerN; j++) {
			rnode_backup[0] = rnode[0];
			swap(rnode[0], j, falsenum);
			iterator_operator(rnode[0], rnode_backup[0]);
		}
		
	}*/

	//初始化种群
	for (int i = 1; i < pop_size; i++)
	{
		rnode[i] = rnode[0];
		for (int j = 0; j < N*2; j++)
		{
			for (int s = 0; s < LayerN; s++)
			{
				swap(rnode[i], s, falsenum);
			}
		}
	}
	

	//开始迭代,0-pop_size属于父代，pop_size-2*pop_size属于子代

	for (int i = 0; i < gen; i++)
	{
		choice_in(choice, pop_size);
		for (int j = 0; j < pop_size / 2; j++)
		{
			crossover_edge(rnode[choice[j][0]], rnode[choice[j][1]], rnodec[j], rnodec[j + pop_size / 2]);
		}
		rnode_backup = rnode;
		rnodec_backup = rnodec;
		choice_in(choice2, pop_size * 2);
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < LayerN; k++)//随机交换边--局部搜索算子
			{
					if ((double)rand() / RAND_MAX < Pm)
					{
						if (choice2[j][0] >= pop_size)
						{
							calculate_eSA(rnodec[choice2[j][0] - pop_size], k);
							if (searchoperator) {
								calculate_SA2(rnodec[choice2[j][0] - pop_size], begin_layer, end_layer);
								local_search_operator2(rnodec[choice2[j][0] - pop_size], k);
							}
							else while (!swap(rnodec[choice2[j][0] - pop_size], k, falsenum));
						}
						else
						{
							calculate_eSA(rnode[choice2[j][0]], k);
							if (searchoperator) {
								calculate_SA2(rnode[choice2[j][0]], begin_layer, end_layer);
								local_search_operator2(rnode[choice2[j][0]], k);
							}
							else while (!swap(rnode[choice2[j][0]], k, falsenum));

						}
						if (choice2[j][1] >= pop_size)
						{
							calculate_eSA(rnodec[choice2[j][1] - pop_size], k);
							if (searchoperator) {
								calculate_SA2(rnodec[choice2[j][1] - pop_size], begin_layer, end_layer);
								local_search_operator2(rnodec[choice2[j][1] - pop_size], k);//
							}
							else while (!swap(rnodec[choice2[j][1] - pop_size], k, falsenum));

						}
						else
						{
							calculate_eSA(rnode[choice2[j][1]], k);
							if (searchoperator) {
								calculate_SA2(rnode[choice2[j][1]], begin_layer, end_layer);
								local_search_operator2(rnode[choice2[j][1]], k);
							}
							else while (!swap(rnode[choice2[j][1]], k, falsenum));

						}
					}
			}
		}	
		for (int j = 0; j < pop_size; j++)//计算所有个体R
		{
			R[j] = computeR(rnode[j], 0, LayerN);
			Rc[j] = computeR(rnodec[j], 0, LayerN);
		}
		for (int j = 0, selected = 0; j < pop_size; j++)//选择算子,只选两个个体
		{
			choice_in(choice2, pop_size * 2);
			selected = select_operator(choice2, pop_size, R, Rc);
			if (selected >= pop_size) {
				rnode_2[j] = rnodec[selected - pop_size];
			//	iterator_operator(rnodec[selected - pop_size], rnodec_backup[selected - pop_size]);
			}
			else {
			//	iterator_operator(rnode[selected], rnode_backup[selected]);
				rnode_2[j] = rnode[selected];
			}
		}
		for (int j = 0; j < pop_size; j++)//产生下一代
		{
			rnode[j] = rnode_2[j];
		}
		//保留最优个体
		for (int j = 0; j < pop_size; j++)
		{
			if (R[j] > R_best)
			{
				R_best = R[j];
			}
			if (Rc[j] > R_best)
			{
				R_best = Rc[j];
			}
		}
		countnum++;
		printf("gen:%d, cycnum:%d , falsenum:%d, R: %f \n", countnum, cycnum, falsenum, R_best);
	}
	change_mark.push_back(check_degree(rnode[0]));
	for (int i = 0; i < rnode[0].size(); i++) {
		Knn[1][i] = computeKnn(rnode[0], 0, i);
	}
	Knn_all.push_back(Knn);
	print_network(rnode[0], 0, 0);
	return R_best;
}

/*	for (int i = 0; i < gen; i++)
{
choice_in(choice,pop_size);
for (int k = 0; k < LayerN; k++)
{
for (int j = 0; j < pop_size / 2; j++)
{
crossover(rnode[choice[j][0]], rnode[choice[j][1]], rnodec[j], rnodec[j + pop_size / 2], k);
activate_node(rnodec[j]);
activate_node(rnodec[j+pop_size/2]);
}
}
choice_in(choice2, pop_size * 2);
for (int k = 0; k < LayerN; k++)//随机交换边--局部搜索算子
{
for (int j = 0; j < pop_size / 2; j++)
{
if (choice2[j][0] >= pop_size)
{
while (!swap(rnodec[choice2[j][0] - pop_size], k, falsenum));
}
else
{
while (!swap(rnode[choice2[j][0]], k, falsenum));
}
if (choice2[j][1] >= pop_size)
{
while (!swap(rnodec[choice2[j][1] - pop_size], k, falsenum));
}
else
{
while (!swap(rnode[choice2[j][1]], k, falsenum));
}
}
}
for (int j = 0; j < pop_size; j++)//计算所有个体R
{
R[j] = computeR(rnode[j], 0, LayerN);
Rc[j] = computeR(rnodec[j], 0, LayerN);
}
for (int j = 0,selected=0; j < pop_size; j++)//选择算子
{
choice_in(choice2, pop_size * 2);
selected = select_operator(choice2, pop_size, R, Rc);
if (selected >= pop_size) {
copy_node(rnodec[selected - pop_size], rnode_2[j]);
}
else {
copy_node(rnode[selected], rnode_2[j]);
}
}
for (int j = 0; j < pop_size; j++)//产生下一代
{
copy_node(rnode_2[j], rnode[j]);
}
//保留最优个体
for (int j = 0; j < pop_size; j++)
{
if (R[j] > R_best)
{
R_best = R[j];
}
if (Rc[j] > R_best)
{
R_best = Rc[j];
}
}

countnum++;
printf("gen:%d, cycnum:%d , falsenum:%d, R: %f \n", countnum, cycnum, falsenum, R_best);
}*/