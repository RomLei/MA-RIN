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
void initialization_ER(vector<node> &rnode,int layer)//初始化ER网络
{
	int edgenumber;//存储边的数量
	int rn, rn2;
	edgenumber = rnode.size() * 2;
	for (int i = 0; i < rnode.size(); i++)
	{
		rn = i;
		rn2 = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
		if (find(rnode[rn].neibor[layer].begin(), rnode[rn].neibor[layer].end(), rn2) == rnode[rn].neibor[layer].end())//若节点rn没有和rn2相连接
		{
			rnode[rn].neibor[layer].push_back(rn2);
			rnode[rn2].neibor[layer].push_back(rn);
			rnode[rn].degree[layer][0]++;
			rnode[rn2].degree[layer][0]++;
			edgenumber--;
		}
	}
	while (edgenumber != 0)
	{
		rn = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
		rn2 = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
		if (find(rnode[rn].neibor[layer].begin(), rnode[rn].neibor[layer].end(), rn2) == rnode[rn].neibor[layer].end())//若节点rn没有和rn2相连接
		{
			rnode[rn].neibor[layer].push_back(rn2);
			rnode[rn2].neibor[layer].push_back(rn);
			rnode[rn].degree[layer][0]++;
			rnode[rn2].degree[layer][0]++;
			edgenumber--;
		}
	}
}

void initialization_SF(vector<node> &rnode, int layer)//初始化SF网络
{
	int degree_ALL=0,degree_MID=0;
	double rn, rn2;
	//创建初始的三个点
	rnode[0].neibor[layer].push_back(1);
	rnode[0].neibor[layer].push_back(2);
	rnode[1].neibor[layer].push_back(0);
	rnode[1].neibor[layer].push_back(2);
	rnode[2].neibor[layer].push_back(0);
	rnode[2].neibor[layer].push_back(1);
	rnode[0].degree[layer][0]+=2;
	rnode[1].degree[layer][0]+=2;
	rnode[2].degree[layer][0]+=2;
	//加入其他点
	for (int i = 3; i < rnode.size(); i++)
	{
		degree_ALL = calculate_degree(rnode, layer);
		rn = (double)rand() / RAND_MAX;
		degree_MID = 0;
		for (int j = 0; j < i; j++)
		{
			degree_MID += rnode[j].degree[layer][0];
			if ((rn <= ((double)degree_MID / degree_ALL)) && (find(rnode[j].neibor[layer].begin(), rnode[j].neibor[layer].end(), i) == rnode[j].neibor[layer].end()))
			{
				rnode[i].neibor[layer].push_back(j);
				rnode[j].neibor[layer].push_back(i);
				rnode[i].degree[layer][0]++;
				rnode[j].degree[layer][0]++;
				rn = (double)rand() / RAND_MAX;
				degree_ALL++;
				j = -1;
				degree_MID = 0;
			}
			if (rnode[i].degree[layer][0] == 2)
			{
				break;
			}
		}
	}
}

void initialization_WS(vector<node> &rnode, int layer)//k=4,pc=0.5
{
	int k = 4, rn = 0;
	double pc = 0.5;
	for (int i = 0; i < rnode.size()-2; i++)
	{
			if ((double)rand() / RAND_MAX < pc)
			{
				do
				{
					rn = (double)rand() / RAND_MAX*(rnode.size() - 1);
				} while (rn==i||find(rnode[i].neibor[layer].begin(),rnode[i].neibor[layer].end(),rn)!=rnode[i].neibor[layer].end());
				connect_node(rnode, layer, i, rn);
			}
			else
			{
				connect_node(rnode, layer, i, i + 2);
			}
			if ((double)rand() / RAND_MAX < pc)
			{
				do
				{
					rn = (double)rand() / RAND_MAX*(rnode.size() - 1);
				} while (rn == i || find(rnode[i].neibor[layer].begin(), rnode[i].neibor[layer].end(), rn) != rnode[i].neibor[layer].end());
				connect_node(rnode, layer, i, rn);
			}
			else
			{
				connect_node(rnode, layer, i, i + 1);
			}
	}
	if ((double)rand() / RAND_MAX < pc)
	{
		do
		{
			rn = (double)rand() / RAND_MAX*(rnode.size() - 1);
		} while (rn == (rnode.size()-1) || find(rnode[rnode.size() - 1].neibor[layer].begin(), rnode[rnode.size() - 1].neibor[layer].end(), rn) != rnode[rnode.size() - 1].neibor[layer].end());
		connect_node(rnode, layer, rnode.size()-1, rn);
	}
	else
	{
		connect_node(rnode, layer, rnode.size() - 1, 1);
	}
	if ((double)rand() / RAND_MAX < pc)
	{
		do
		{
			rn = (double)rand() / RAND_MAX*(rnode.size() - 1);
		} while (rn == (rnode.size() - 1) || find(rnode[rnode.size() - 1].neibor[layer].begin(), rnode[rnode.size() - 1].neibor[layer].end(), rn) != rnode[rnode.size() - 1].neibor[layer].end());
		connect_node(rnode, layer, rnode.size() - 1, rn);
	}
	else
	{
		connect_node(rnode, layer, rnode.size() - 1, 0);
	}
	if ((double)rand() / RAND_MAX < pc)
	{
		do
		{
			rn = (double)rand() / RAND_MAX*(rnode.size() - 1);
		} while (rn == (rnode.size()-2) || find(rnode[rnode.size()-2].neibor[layer].begin(), rnode[rnode.size()-2].neibor[layer].end(), rn) != rnode[rnode.size()-2].neibor[layer].end());
		connect_node(rnode, layer, rnode.size() - 2, rn);
	}
	else
	{
		connect_node(rnode, layer, rnode.size()-2, 0);
	}
	if ((double)rand() / RAND_MAX < pc)
	{
		do
		{
			rn = (double)rand() / RAND_MAX*(rnode.size() - 1);
		} while (rn == (rnode.size() - 2) || find(rnode[rnode.size() - 2].neibor[layer].begin(), rnode[rnode.size() - 2].neibor[layer].end(), rn) != rnode[rnode.size() - 2].neibor[layer].end());
		connect_node(rnode, layer, rnode.size() - 2, rn);
	}
	else
	{
		connect_node(rnode, layer, rnode.size() - 2, rnode.size()-1);
	}
	
}

void initialization_NW(vector<node> &rnode, int layer)//k=4,pa=0.01
{
	double pa = 0.01;
	for (int i = 0; i < rnode.size() - 2; i++)
	{
	
			connect_node(rnode, layer, i, i + 2);
			connect_node(rnode, layer, i, i + 1);
	}
	connect_node(rnode, layer, rnode.size() - 1, 1);
	connect_node(rnode, layer, rnode.size() - 1, 0);
	connect_node(rnode, layer, rnode.size() - 2, 0);
	connect_node(rnode, layer, rnode.size() - 2, rnode.size() - 1);
	for (int i = 0,x=rnode.size()-1; i < rnode.size(); i++,x--)
	{
		for (int j = 1; j <= x; j++)
		{
			if (((double)rand() / RAND_MAX < pa) && (find(rnode[i].neibor[layer].begin(), rnode[i].neibor[layer].end(), i + j) == rnode[i].neibor[layer].end()))
			{
				connect_node(rnode, layer, i, i + j);
			}
		}
	}
}

void initialization_degree(vector<node> &rnode)//初始化度
{
	for (int i = 0; i < rnode.size(); i++)
	{
		for (int j = 0; j < LayerN; j++)
		{
			rnode[i].degree[j].push_back(0);
		}
	}
}

int calculate_degree(vector<node> &rnode, int layer)//计算总的度
{
	int degree_ALL = 0;
	for (int i = 0; i < rnode.size(); i++)
	{
		degree_ALL += rnode[i].degree[layer][0];
	}
	return degree_ALL;
}

void activate_node(vector<node> &rnode)//激活节点
{
	for (int i = 0; i < rnode.size(); i++)
	{
		rnode[i].activation = 1;
	}
}

void connect_node(vector<node> &rnode, int layer,int node1,int node2)//连接两个节点
{
	rnode[node1].neibor[layer].push_back(node2);
	rnode[node2].neibor[layer].push_back(node1);
	rnode[node1].degree[layer][0]++;
	rnode[node2].degree[layer][0]++;
}

void initialSA(vector<node> &rnode) 
{
	for (int i = 0; i < rnode.size(); i++) {
		for (int j = 0; j < LayerN; j++) {
			rnode[i].SA.push_back(1.0 / rnode.size());
			rnode[i].eSA[j].resize(rnode[i].neibor[j].size(), 0);
		}	
	}
}