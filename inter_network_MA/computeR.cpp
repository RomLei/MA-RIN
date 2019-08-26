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
int lcc_code=0;
int largest_lcc_code = 0;
double computeR(vector<node> &rnode,int begin_layer,int end_layer)//��������³���ԣ�ֻ��������һ������,��begin��end��
{
	int maxdegreeP,lcc_size=0, layer_ops;
	double R = 0;
	vector<node> rnodeB(rnode.size());//�������ݽڵ�
	rnodeB = rnode;
	vector<int> res_node(LayerN, 0);
	while (calculate_res_node(rnodeB)!=0)
	{
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
		} while (res_node[0]!=res_node[1]);
		R += (double)lcc_size / rnode.size();
	}
	R = R / rnode.size();
	return R;
}

inline int max_degree_P(vector<node> &rnode, int layer)//�ҵ������Ľڵ�
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

inline void attack_NODE(vector<node> &rnode,int layer, int attackP)//�����ڵ�
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

void copy_node(vector<node> &rnode, vector<node> &rnodeB)//����rnode��Ϣ��rnodeB��
{
	for (int i = 0; i < rnode.size(); i++)
	{
		rnodeB[i].activation = rnode[i].activation;
		rnodeB[i].largest_connect_cluster = rnode[i].largest_connect_cluster;
		for (int j = 0; j < LayerN; j++)
		{
			rnodeB[i].SA[j] = rnode[i].SA[j];
			rnodeB[i].eSA[j] = rnode[i].eSA[j];
			rnodeB[i].degree[j] = rnode[i].degree[j];
			rnodeB[i].neibor[j] = rnode[i].neibor[j];
		}
	}
}

inline int find_largestconnectcluster(vector<node> &rnode,int layer)
{
	int N = rnode.size();
	vector<int> search(N);//�洢�ڵ��Ƿ񱻼���������Ϣ
	int maxconnectnum = 1;
	int p = 0;
	int IN = 0;
	for (int i = 0; i < N; i++)
	{
		search[i] = i;
	}
	while (whether_search_end(search, rnode))
	{
		while (search[p] == -1 || rnode[p].activation == 0)
		{
			p++;
		}
		IN = Include_Nodes(p,layer, search, rnode);
		if (maxconnectnum < IN)
		{
			maxconnectnum = IN;
			largest_lcc_code = lcc_code;
		}
		lcc_code++;
	}
	return maxconnectnum;
}

inline bool whether_search_end(vector<int> &search, vector<node> &rnode)//����false˵��ȫ���㶼���������ˣ�����true���е�δ������
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

inline int Include_Nodes(int p, int layer, vector<int> &search, vector<node> &rnode)//����ڵ�λ��p�������ýڵ����ӵ�δ���������Ľڵ�����������ֵ
{
	int N = rnode.size();
	int res;
	int neibornode;
	if (search[p] != -1 && rnode[p].activation != 0)
	{
		search[p] = -1;
		rnode[p].largest_connect_cluster = lcc_code;
		int maxconnectnum = 1;
		for (int i = 0; i < rnode[p].neibor[layer].size(); i++)
		{
			neibornode = rnode[p].neibor[layer][i];
			if (rnode[neibornode].activation==1 && search[neibornode] != -1)
			{
				maxconnectnum += Include_Nodes(neibornode, layer, search, rnode);
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

inline void clear_lcc_flag(vector<node> &rnode)//���������Ӵر��
{
	for (int i = 0; i < rnode.size(); i++)
	{
		rnode[i].largest_connect_cluster = -1;
	}
	largest_lcc_code = 0;
	lcc_code = 0;
}

inline void delete_node(vector<node> &rnode,int layer)//ɾ������������Ӵ��еĽڵ�
{
	for (int i = 0; i < rnode.size(); i++)
	{
		if (rnode[i].largest_connect_cluster != largest_lcc_code)
		{
			attack_NODE(rnode, layer,i);
		}
	}
}

inline int calculate_res_node(vector<node> &rnode)
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

double computeR_single(vector<node> &rnode, int layer)//���㵥������³����
{
	int maxdegreeP, lcc_size = 0, layer_ops;
	double R = 0;
	vector<node> rnodeB;//�������ݽڵ�
	rnodeB.resize(rnode.size());
	rnodeB = rnode;
	int res_node1 = 0, res_node2 = 0;
	while (calculate_res_node(rnodeB) != 0)
	{
		maxdegreeP = max_degree_P(rnodeB, layer);
		attack_NODE(rnodeB, layer, maxdegreeP);
		lcc_size = find_largestconnectcluster(rnodeB, layer);
		R += (double)lcc_size / rnode.size();
	}
	R = R / rnode.size();
	return R;
}