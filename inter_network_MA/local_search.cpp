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
int lcc_code_local = 0;
int largest_lcc_code_local = 0;
void calculate_eSA(vector<node> &rnode, int layer)//计算eSA
{
	for (int i = 0; i < rnode.size(); i++) {
		for (int j = 0; j < rnode[i].neibor[layer].size(); j++) {
			rnode[i].eSA[layer][j] = 1 / (1 + (double)rnode.size()*fabs(rnode[i].SA[layer] - rnode[rnode[i].neibor[layer][j]].SA[layer])*exp(-(rnode[i].SA[layer] + rnode[rnode[i].neibor[layer][j]].SA[layer])));
		}
	}
}

void mark_lived_node(vector<node> &rnode_back,vector<node> &rnode)//记录存活节点，并+1
{
	for (int i = 0; i < rnode.size(); i++) {
		if (rnode_back[i].activation == 1) {
			rnode[i].SA2++;
		}
	}
}

void calculate_SA2(vector<node> &rnode,int begin_layer, int end_layer)
{
	int maxdegreeP, lcc_size = 0, layer_ops;
	double R = 0;
	int round = 1;
	vector<node> rnodeB(rnode.size());//创建备份节点
	rnodeB = rnode;
	vector<int> res_node(LayerN, 0);
	while (calculate_res_node(rnodeB) != 0)
	{
		mark_lived_node(rnodeB, rnode);
		maxdegreeP = max_degree_P(rnodeB, begin_layer);
		for (int j = begin_layer; j < end_layer; ++j)
		{
			attack_NODE(rnodeB, j, maxdegreeP);
		}
		do
		{
			res_node[0] = calculate_res_node(rnodeB);
			for (int i = begin_layer, layers = begin_layer; i < end_layer; ++i, ++layers)
			{
				lcc_size = find_largestconnectcluster(rnodeB, layers);
				for (int j = begin_layer; j < end_layer; ++j)
				{
					delete_node(rnodeB, j);
				}
				clear_lcc_flag(rnodeB);
			}
			res_node[1] = calculate_res_node(rnodeB);
		} while (res_node[0] != res_node[1]);
		round++;
	}
	for (int i = 0; i < rnode.size(); i++) {
		rnode[i].SA2 = rnode[i].SA2 / round;
	}
}

double min_SA2(vector<node> rnode, int node1, int node2)
{
	if (rnode[node1].SA2 < rnode[node2].SA2)
	{
		return rnode[node1].SA2;
	}
	else {
		return rnode[node2].SA2;
	}
}

void local_search_operator2(vector<node>& rnode, int layer)//采用后向传播算法2的局部搜索算子
{
	double sumSA = 0, ppi, accumulaterate = 0;
	vector<node> rnode_backup;
	vector<double> Sijkl(2, 0);
	vector<double> Sikjl(2, 0);
	vector<int> n;
	int aa = 0;
	bool flag = true;
	while (flag) {
	
		// 随机选两个边
		n.push_back(((double)rand() / RAND_MAX)*(rnode.size() - 1));
		n.push_back(((double)rand() / RAND_MAX)*(rnode[n[0]].neibor[layer].size() - 1));
		n[1] = rnode[n[0]].neibor[layer][n[1]];
		n.push_back(((double)rand() / RAND_MAX)*(rnode.size() - 1));
		n.push_back(((double)rand() / RAND_MAX)*(rnode[n[2]].neibor[layer].size() - 1));
		n[3] = rnode[n[2]].neibor[layer][n[3]];

		Sikjl[0] = min_SA2(rnode, n[0], n[1]);
		Sikjl[1] = min_SA2(rnode, n[2], n[3]);
		Sijkl[0] = min_SA2(rnode, n[0], n[2]);
		Sijkl[1] = min_SA2(rnode, n[1], n[3]);
		flag = false;
		if (find(rnode[n[0]].neibor[layer].begin(), rnode[n[0]].neibor[layer].end(), n[2]) != rnode[n[0]].neibor[layer].end())
		{
			flag = true;
			n.clear();
		}
		else if (find(rnode[n[3]].neibor[layer].begin(), rnode[n[3]].neibor[layer].end(), n[1]) != rnode[n[3]].neibor[layer].end())
		{
			n.clear();
			flag = true;
		}
		else if (n[0] == n[2] || n[1] == n[3] || n[0] == n[3] || n[2] == n[1])
		{
			n.clear();
			flag = true;
		}
		if (Sijkl[0] + Sijkl[1]>Sikjl[0] + Sikjl[1]) {
			n.clear();
			flag = true;
		}
	}
	/*	rnode_backup = rnode;
	rnode_backup[n[0]].neibor[layer].erase(find(rnode_backup[n[0]].neibor[layer].begin(), rnode_backup[n[0]].neibor[layer].end(), n[1]));
	rnode_backup[n[1]].neibor[layer].erase(find(rnode_backup[n[1]].neibor[layer].begin(), rnode_backup[n[1]].neibor[layer].end(), n[0]));
	rnode_backup[n[2]].neibor[layer].erase(find(rnode_backup[n[2]].neibor[layer].begin(), rnode_backup[n[2]].neibor[layer].end(), n[3]));
	rnode_backup[n[3]].neibor[layer].erase(find(rnode_backup[n[3]].neibor[layer].begin(), rnode_backup[n[3]].neibor[layer].end(), n[2]));
	rnode_backup[n[0]].neibor[layer].push_back(n[2]);
	rnode_backup[n[2]].neibor[layer].push_back(n[0]);
	rnode_backup[n[3]].neibor[layer].push_back(n[1]);
	rnode_backup[n[1]].neibor[layer].push_back(n[3]);*/
	//	if (computeR(rnode, 0, LayerN) <= computeR(rnode_backup, 0, LayerN)) {
	rnode[n[0]].neibor[layer].erase(find(rnode[n[0]].neibor[layer].begin(), rnode[n[0]].neibor[layer].end(), n[1]));
	rnode[n[1]].neibor[layer].erase(find(rnode[n[1]].neibor[layer].begin(), rnode[n[1]].neibor[layer].end(), n[0]));
	rnode[n[2]].neibor[layer].erase(find(rnode[n[2]].neibor[layer].begin(), rnode[n[2]].neibor[layer].end(), n[3]));
	rnode[n[3]].neibor[layer].erase(find(rnode[n[3]].neibor[layer].begin(), rnode[n[3]].neibor[layer].end(), n[2]));
	rnode[n[0]].neibor[layer].push_back(n[2]);
	rnode[n[2]].neibor[layer].push_back(n[0]);
	rnode[n[3]].neibor[layer].push_back(n[1]);
	rnode[n[1]].neibor[layer].push_back(n[3]);
	//	}

}

void local_search_operator(vector<node>& rnode,int layer)//采用后向传播算法的局部搜索算子
{
	double sumSA=0,ppi,accumulaterate=0;
	vector<node> rnode_backup;
	vector<double> Sijkl(2,0);
	vector<double> Sikjl(2,0);
	vector<int> n;
	int aa=0;
	bool flag = true;
	for (int i = 0; i < rnode.size(); i++) {
		for (int j = 0; j < rnode[i].eSA[layer].size(); j++) {
			sumSA += rnode[i].eSA[layer][j];
		}
	}
	if (sumSA == 0) {
		while (!swap(rnode, layer, aa));
			return;
	}
	while (flag) {
	/*	for (int i = 0; i < 2; i++) {
			ppi = ((double)rand() / (RAND_MAX+1));
			accumulaterate = 0;
			for (int j = 0; j < rnode.size(); j++) {
				for (int k = 0; k < rnode[j].eSA[layer].size(); k++) {
					accumulaterate += (rnode[j].eSA[layer][k] / sumSA);
					if (accumulaterate >= ppi) {
						Sijkl[i]=(rnode[j].eSA[layer][k]);
						n.push_back(j); n.push_back(rnode[j].neibor[layer][k]);
						j = rnode.size();
						break;
					}
					if (j == rnode.size() - 1 && k == rnode[j].eSA[layer].size() - 1) {
						Sijkl[i]=(rnode[j].eSA[layer][k]);
						n.push_back(j); n.push_back(rnode[j].neibor[layer][k]);
						j = rnode.size();
						break;
					}
				}
			}
		}*/
		// 随机选两个边
		n .push_back(((double)rand() / RAND_MAX)*(rnode.size() - 1));
		n.push_back( ((double)rand() / RAND_MAX)*(rnode[n[0]].neibor[layer].size() - 1));
		n[1] = rnode[n[0]].neibor[layer][n[1]];
		n.push_back( ((double)rand() / RAND_MAX)*(rnode.size() - 1));
		n.push_back(((double)rand() / RAND_MAX)*(rnode[n[2]].neibor[layer].size() - 1));
		n[3] = rnode[n[2]].neibor[layer][n[3]];
		
		Sikjl[0]=1/(1+ (double)rnode.size()*fabs(rnode[n[0]].SA[layer] - rnode[n[2]].SA[layer]) *exp(-rnode[n[0]].SA[layer] - rnode[n[2]].SA[layer]));
		Sikjl[1]=1/(1+ (double)rnode.size()*fabs(rnode[n[1]].SA[layer] - rnode[n[3]].SA[layer]) *exp (-rnode[n[1]].SA[layer] - rnode[n[3]].SA[layer]));
		Sijkl[0]= 1 / (1 + (double)rnode.size()*fabs(rnode[n[0]].SA[layer] - rnode[n[1]].SA[layer]) *exp(-rnode[n[0]].SA[layer] - rnode[n[1]].SA[layer]));
		Sijkl[1]= 1 / (1 + (double)rnode.size()*fabs(rnode[n[2]].SA[layer] - rnode[n[3]].SA[layer]) *exp(-rnode[n[2]].SA[layer] - rnode[n[3]].SA[layer]));
		flag = false;
		if (find(rnode[n[0]].neibor[layer].begin(), rnode[n[0]].neibor[layer].end(), n[2]) != rnode[n[0]].neibor[layer].end())
		{
			flag=true;
			n.clear();
		}
		else if (find(rnode[n[3]].neibor[layer].begin(), rnode[n[3]].neibor[layer].end(), n[1]) != rnode[n[3]].neibor[layer].end())
		{
			n.clear();
			flag = true;
		}
		else if (n[0] == n[2] || n[1] == n[3] || n[0] == n[3] || n[2] == n[1])
		{
			n.clear();
			flag = true;
		}
		if (Sijkl[0]+Sijkl[1]>Sikjl[0]+Sikjl[1]) {
			n.clear();
			flag = true;
		}
	}
/*	rnode_backup = rnode;
	rnode_backup[n[0]].neibor[layer].erase(find(rnode_backup[n[0]].neibor[layer].begin(), rnode_backup[n[0]].neibor[layer].end(), n[1]));
	rnode_backup[n[1]].neibor[layer].erase(find(rnode_backup[n[1]].neibor[layer].begin(), rnode_backup[n[1]].neibor[layer].end(), n[0]));
	rnode_backup[n[2]].neibor[layer].erase(find(rnode_backup[n[2]].neibor[layer].begin(), rnode_backup[n[2]].neibor[layer].end(), n[3]));
	rnode_backup[n[3]].neibor[layer].erase(find(rnode_backup[n[3]].neibor[layer].begin(), rnode_backup[n[3]].neibor[layer].end(), n[2]));
	rnode_backup[n[0]].neibor[layer].push_back(n[2]);
	rnode_backup[n[2]].neibor[layer].push_back(n[0]);
	rnode_backup[n[3]].neibor[layer].push_back(n[1]);
	rnode_backup[n[1]].neibor[layer].push_back(n[3]);*/
//	if (computeR(rnode, 0, LayerN) <= computeR(rnode_backup, 0, LayerN)) {
		rnode[n[0]].neibor[layer].erase(find(rnode[n[0]].neibor[layer].begin(), rnode[n[0]].neibor[layer].end(), n[1]));
		rnode[n[1]].neibor[layer].erase(find(rnode[n[1]].neibor[layer].begin(), rnode[n[1]].neibor[layer].end(), n[0]));
		rnode[n[2]].neibor[layer].erase(find(rnode[n[2]].neibor[layer].begin(), rnode[n[2]].neibor[layer].end(), n[3]));
		rnode[n[3]].neibor[layer].erase(find(rnode[n[3]].neibor[layer].begin(), rnode[n[3]].neibor[layer].end(), n[2]));
		rnode[n[0]].neibor[layer].push_back(n[2]);
		rnode[n[2]].neibor[layer].push_back(n[0]);
		rnode[n[3]].neibor[layer].push_back(n[1]);
		rnode[n[1]].neibor[layer].push_back(n[3]);
//	}

}

inline int max_degree_P_local(vector<node> &rnode, int layer)//找到度最大的节点
{
	int P = 0, maxdegree = 0;
	for (int i = 0; i < rnode.size(); i++)
	{
		if ((rnode[i].degree[layer][0] > maxdegree) && (rnode[i].activation == 1))
		{
			P = i;
			maxdegree = rnode[i].degree[layer][0];
		}
	}
	if (maxdegree == 0)
	{
		for (int i = 0; i < rnode.size(); i++)
		{
			if (rnode[i].activation == 1)
			{
				P = i;
				return P;
			}
		}
	}
	return P;
}

inline void attack_NODE_local(vector<node> &rnode, int layer, int attackP)//攻击节点
{
	rnode[attackP].activation = 0;
	int numbernode = 0;
	for (int i = 0; i < rnode[attackP].neibor[layer].size(); i++)
	{
		numbernode = rnode[attackP].neibor[layer][i];
		rnode[numbernode].neibor[layer].erase(find(rnode[numbernode].neibor[layer].begin(), rnode[numbernode].neibor[layer].end(), attackP));
		rnode[numbernode].degree[layer][0]--;
	}
	rnode[attackP].neibor[layer].clear();
	rnode[attackP].degree[layer][0] = 0;
}
inline bool whether_search_end_local(vector<int> &search, vector<node> &rnode)//返回false说明全部点都被检索过了，返回true则还有点未被检索
{
	int N = rnode.size();
	bool u = false;
	for (int i = 0; i < N; i++)
	{
		if (search[i] != -1 && rnode[i].activation != 0)
		{
			u = true;
		}
	}
	return u;
}

inline int Include_Nodes_local(int p, int layer, vector<int> &search, vector<node> &rnode)//输入节点位置p，检索该节点连接的未被检索过的节点数，返回其值
{
	int N = rnode.size();
	int res;
	int neibornode;
	if (search[p] != -1 && rnode[p].activation != 0)
	{
		search[p] = -1;
		rnode[p].largest_connect_cluster = lcc_code_local;
		int maxconnectnum = 1;
		for (int i = 0; i < rnode[p].neibor[layer].size(); i++)
		{
			neibornode = rnode[p].neibor[layer][i];
			if (rnode[neibornode].activation == 1 && search[neibornode] != -1)
			{
				maxconnectnum += Include_Nodes_local(neibornode, layer, search, rnode);
			}
		}
		res = maxconnectnum;
	}
	else
	{
		res = 0;
	}
	return res;
}

inline int find_largestconnectcluster_local(vector<node> &rnode, int layer)
{
	int N = rnode.size();
	vector<int> search(N);//存储节点是否被检索过的信息
	int maxconnectnum = 1;
	int p = 0;
	int IN = 0;
	for (int i = 0; i < N; i++)
	{
		search[i] = i;
	}
	while (whether_search_end_local(search, rnode))
	{
		while (search[p] == -1 || rnode[p].activation == 0)
		{
			p++;
		}
		IN = Include_Nodes_local(p, layer, search, rnode);
		if (maxconnectnum < IN)
		{
			maxconnectnum = IN;
			largest_lcc_code_local = lcc_code_local;
		}
		lcc_code_local++;
	}
	return maxconnectnum;
}





inline void clear_lcc_flag_local(vector<node> &rnode)//清除最大连接簇标记
{
	for (int i = 0; i < rnode.size(); i++)
	{
		rnode[i].largest_connect_cluster = -1;
	}
	largest_lcc_code_local = 0;
	lcc_code_local = 0;
}

inline void delete_node_local(vector<node> &rnode, int layer)//删除不是最大连接簇中的节点
{
	for (int i = 0; i < rnode.size(); i++)
	{
		if (rnode[i].largest_connect_cluster != largest_lcc_code_local)
		{
			attack_NODE_local(rnode, layer, i);
		}
	}
}

inline int calculate_res_node_local(vector<node> &rnode)
{
	int res_node = rnode.size();
	for (int i = 0; i < rnode.size(); i++)
	{
		if (rnode[i].activation == 0)
		{
			res_node--;
		}
	}
	return res_node;
}

void iterator_operator(vector<node> &rnode, vector<node> rnode_before)//迭代算子
{
	int maxdegreeP, lcc_size = 0, layer_ops;
	int begin_layer = 0;//攻击起始层
	int end_layer = LayerN;//攻击结束层
	vector<node> rnode_backup=rnode;
	vector<vector<double>> deta_SAi;
	deta_SAi.resize(LayerN, vector<double>(rnode.size(), 0));
	int res_node[2] = { 0,0 };
	size_t length = rnode.size();
	int sum_E0,sum_E1;
	double EE0, EE1;
	for (int i = 0; i < length; i++) {
		for (size_t layer = 0; layer < LayerN; layer++)
		{
			sum_E0 = 0; sum_E1 = 0;
			for (size_t j = 0; j < length; j++)
			{
				sum_E0 += rnode_backup[j].degree[layer][0];
				sum_E1 += rnode_before[j].degree[layer][0];
			}
			for (size_t j = 0; j < length; j++)
			{
				if (sum_E0 == 0) {
					EE0 = 0;
				}
				else {
					EE0 = (double)rnode_backup[j].degree[layer][0] / sum_E0;
				}
				if (sum_E1 == 0) {
					EE1 = 0;
				}
				else {
					EE1 = (double)rnode_before[j].degree[layer][0] / sum_E1;
				}
				deta_SAi[layer][j]+=(EE0 - EE1);
			}
		}
		if (calculate_res_node_local(rnode_backup)) {
			maxdegreeP = max_degree_P_local(rnode_backup, begin_layer);
			for (int j = begin_layer; j < end_layer; j++)//攻击
			{
				attack_NODE_local(rnode_backup, j, maxdegreeP);
			}
			do
			{

				res_node[0] = calculate_res_node_local(rnode_backup);
				for (int i = begin_layer, layers = begin_layer; i < end_layer; i++, layers++)
				{
					lcc_size = find_largestconnectcluster_local(rnode_backup, layers);
					for (int j = begin_layer; j < end_layer; j++)
					{
						delete_node_local(rnode_backup, j);
					}
					clear_lcc_flag_local(rnode_backup);
				}
				res_node[1] = calculate_res_node_local(rnode_backup);
			} while (res_node[0] != res_node[1]);
		}
		if (calculate_res_node_local(rnode_before)) {
			maxdegreeP = max_degree_P_local(rnode_before, begin_layer);
			for (int j = begin_layer; j < end_layer; j++)//攻击
			{
				attack_NODE_local(rnode_before, j, maxdegreeP);
			}
			do
			{

				res_node[0] = calculate_res_node_local(rnode_before);
				for (int i = begin_layer, layers = begin_layer; i < end_layer; i++, layers++)
				{
					lcc_size = find_largestconnectcluster_local(rnode_before, layers);
					for (int j = begin_layer; j < end_layer; j++)
					{
						delete_node_local(rnode_before, j);
					}
					clear_lcc_flag_local(rnode_before);
				}
				res_node[1] = calculate_res_node_local(rnode_before);
			} while (res_node[0] != res_node[1]);
		}
		else if (!calculate_res_node_local(rnode_backup)) {
			break;
		}

	}
	for (size_t i = 0; i < length; i++)
	{
		for (size_t j = 0; j < LayerN; j++)
		{
			rnode[i].SA[j] += 1/(1+exp(-namuta*deta_SAi[j][i]))-0.5;
		}
	}
}