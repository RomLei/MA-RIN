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
bool swap(vector<node> &rnode, int layer,int &falsenum)//随机交换一层中的两条边,返回false说明交换失败，true成功
{
	int node1, node1_neibor, node2, node2_neibor;
	node1 = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
	node1_neibor = ((double)rand() / RAND_MAX)*(rnode[node1].neibor[layer].size()-1);
	node1_neibor = rnode[node1].neibor[layer][node1_neibor];
	node2 = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
	node2_neibor = ((double)rand() / RAND_MAX)*(rnode[node2].neibor[layer].size()-1);
	node2_neibor = rnode[node2].neibor[layer][node2_neibor];
	if (find(rnode[node1].neibor[layer].begin(), rnode[node1].neibor[layer].end(), node2_neibor) != rnode[node1].neibor[layer].end())
	{
		falsenum++;
		return false;
	}
	else if(find(rnode[node2].neibor[layer].begin(), rnode[node2].neibor[layer].end(), node1_neibor) != rnode[node2].neibor[layer].end())
	{
		falsenum++;
		return false;
	}
	else if (node1 == node2 || node1_neibor == node2 || node1 == node2_neibor || node2 == node1_neibor)
	{
		falsenum++;
		return false;
	}
	else
	{
		rnode[node1].neibor[layer].erase(find(rnode[node1].neibor[layer].begin(), rnode[node1].neibor[layer].end(), node1_neibor));
		rnode[node1_neibor].neibor[layer].erase(find(rnode[node1_neibor].neibor[layer].begin(), rnode[node1_neibor].neibor[layer].end(), node1));
		rnode[node2].neibor[layer].erase(find(rnode[node2].neibor[layer].begin(), rnode[node2].neibor[layer].end(), node2_neibor));
		rnode[node2_neibor].neibor[layer].erase(find(rnode[node2_neibor].neibor[layer].begin(), rnode[node2_neibor].neibor[layer].end(), node2));
		rnode[node1].neibor[layer].push_back(node2_neibor);
		rnode[node2_neibor].neibor[layer].push_back(node1);
		rnode[node2].neibor[layer].push_back(node1_neibor);
		rnode[node1_neibor].neibor[layer].push_back(node2);
	}
	return true;
}