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

bool check_degree(vector<node> rnode)//检查度是否发生变化
{
	for (int i = 0; i < rnode.size(); i++)
	{
		for (int j = 0; j < LayerN; j++)
		{
			if (rnode[i].neibor[j].size() != rnode[i].degree_record[j][0])
			{
				return false;
			}
		}
	}
	return true;
}

void mark_degree(vector<node> &rnode)//记录度
{
	for (int i = 0; i < rnode.size(); i++)
	{
		for (int j = 0; j < LayerN; j++)
		{
			rnode[i].degree_record[j].push_back(rnode[i].neibor[j].size());
		}
	}
}

/*void print_network(vector<node> rnode, int layer)//打印单层网络
{
	ofstream PrintNetwork;
	vector<vector<int>> G;
	G.resize(rnode.size(), vector<int>(rnode.size(), 0));
	PrintNetwork.open("singlenetwork.txt");
	for (int i = 0; i < rnode.size(); i++)
	{
		for (int j = 0; j < rnode.size(); j++)
		{
			G[i][j] = 0;
		}
	}
	for (int i = 0; i < rnode.size(); i++)
	{
		for (int j = 0; j < rnode[i].neibor[layer].size(); j++)
		{
			G[i][rnode[i].neibor[layer][j]] = 1;
		}
	}
	for (int i = 0; i < rnode.size(); i++)
	{
		for (int j = 0; j < rnode.size(); j++)
		{
			PrintNetwork << G[i][j];
		}
		PrintNetwork << "\n";
	}
	PrintNetwork.close();
}*/

void choice_in(vector<vector<int>> &choice,int pop_size)//将种群两两配对
{
	vector<int> popnum;
	int ppi;
	for (int i = 0; i < pop_size; i++)
	{
		popnum.push_back(i);
	}
	for (int i = 0; i < pop_size / 2; i++)
	{
		ppi = ((double)rand() / RAND_MAX)*(popnum.size() - 1);
		choice[i][0] = popnum[ppi];
		popnum.erase(popnum.begin()+ppi);
		ppi = ((double)rand() / RAND_MAX)*(popnum.size() - 1);
		choice[i][1] = popnum[ppi];
		popnum.erase(popnum.begin() + ppi);
	}
}

void choice_in2(vector<vector<int>> &choice, int pop_size)//选出两个个体
{
	vector<int> popnum;
	int ppi;
	for (int i = 0; i < pop_size; i++)
	{
		popnum.push_back(i);
	}
	for (int i = 0; i < 1; i++)
	{
		ppi = ((double)rand() / RAND_MAX)*(popnum.size() - 1);
		choice[i][0] = popnum[ppi];
		popnum.erase(popnum.begin() + ppi);
		ppi = ((double)rand() / RAND_MAX)*(popnum.size() - 1);
		choice[i][1] = popnum[ppi];
		popnum.erase(popnum.begin() + ppi);
	}
}

int select_operator(vector<vector<int>> &choice,int pop_size,vector<double> R,vector<double> Rc)//找到choice中一半个体中R最大的个体
{
	double R_max = 0, R_transmit = 0;
	int maxpos=choice[0][0];
	for (int k = 0; k < pop_size/2; k++)//找一半的个体
	{
		if (choice[k][0] >= pop_size) {
			R_transmit = Rc[choice[k][0] - pop_size];
		}
		else {
			R_transmit = R[choice[k][0]];
		}
		if (R_transmit > R_max)
		{
			R_max = R_transmit;
			maxpos = choice[k][0];
		}
		if (choice[k][1] >= pop_size) {
			R_transmit = Rc[choice[k][1] - pop_size];
		}
		else {
			R_transmit = R[choice[k][1]];
		}
		if (R_transmit > R_max)
		{
			R_max = R_transmit;
			maxpos = choice[k][1];
		}
	}
	return maxpos;
}

int select_operator2(vector<vector<int>> &choice, int pop_size, vector<double> R, vector<double> Rc)//找到choice个体中R最大的个体
{
	double R_max = 0, R_transmit = 0;
	int maxpos = choice[0][0];
	for (int k = 0; k < 1; k++)//找一半的个体
	{
		if (choice[k][0] >= pop_size) {
			R_transmit = Rc[choice[k][0] - pop_size];
		}
		else {
			R_transmit = R[choice[k][0]];
		}
		if (R_transmit > R_max)
		{
			R_max = R_transmit;
			maxpos = choice[k][0];
		}
		if (choice[k][1] >= pop_size) {
			R_transmit = Rc[choice[k][1] - pop_size];
		}
		else {
			R_transmit = R[choice[k][1]];
		}
		if (R_transmit > R_max)
		{
			R_max = R_transmit;
			maxpos = choice[k][1];
		}
	}
	return maxpos;
}

void copy_node_layer(vector<node> &rnode1, vector<node> &rnode2, int layer)//复制某一层的网络
{
	for (int i = 0; i < rnode1.size(); i++)
	{
		rnode2[i].degree[layer] = rnode1[i].degree[layer];
		rnode2[i].neibor[layer] = rnode1[i].neibor[layer];
	}
}