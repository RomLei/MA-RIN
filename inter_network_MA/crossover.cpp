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

void crossover(vector<node> &rnode_p1, vector<node> &rnode_p2, vector<node> &rnode_c1, vector<node> &rnode_c2,int layer)//交叉算子
{
	int ppi = 0, ppi2 = 0, t_node = 0, s_node = 0;
	vector<int> target_nonshare_node;
	vector<int> target_nonshare_node2;
	vector<int> target_share_node;
	copy_node_layer(rnode_p1, rnode_c1,layer);
	copy_node_layer(rnode_p2, rnode_c2,layer);
	for (int i = 0; i < rnode_p1.size(); i++)
	{
		if ((double)rand() / RAND_MAX < Pc)//小于交叉概率才交叉
		{
			find_targetnode(rnode_c1, rnode_c2, i, target_nonshare_node, target_share_node, layer);// 找到rnode2的邻居中与rnode1的非共有邻节点与共有节点
			find_targetnode(rnode_c2, rnode_c1, i, target_nonshare_node2, target_share_node, layer);//找到rnode1的邻居中与rnode2的非共有邻节点与共有节点
			//rnode_c1
			delete_neibor(rnode_c1, i, target_share_node, layer);//删除共有节点
			for (int j = 0,js=0,jj= rnode_c1[i].neibor[layer].size(); j < jj; j++)//删除非共有节点
			{
				t_node = rnode_c1[i].neibor[layer][js];
				ppi = ((double)rand() / RAND_MAX)*(target_nonshare_node.size() - 1);
				ppi2 = ((double)rand() / RAND_MAX)*(rnode_c1[target_nonshare_node[ppi]].neibor[layer].size() - 1);
				s_node = rnode_c1[target_nonshare_node[ppi]].neibor[layer][ppi2];
				if ((find(rnode_c1[t_node].neibor[layer].begin(), rnode_c1[t_node].neibor[layer].end(), s_node) == rnode_c1[t_node].neibor[layer].end())&&(t_node!=s_node))
				{
					rnode_c1[i].neibor[layer].erase(rnode_c1[i].neibor[layer].begin()+js);//删除非共有节点
					rnode_c1[t_node].neibor[layer].erase(find(rnode_c1[t_node].neibor[layer].begin(), rnode_c1[t_node].neibor[layer].end(), i));
					rnode_c1[t_node].neibor[layer].push_back(s_node);
					rnode_c1[s_node].neibor[layer].push_back(t_node);
					rnode_c1[i].neibor[layer].push_back(target_nonshare_node[ppi]);//添加非共有节点
					rnode_c1[target_nonshare_node[ppi]].neibor[layer].push_back(i);
					rnode_c1[s_node].neibor[layer].erase(find(rnode_c1[s_node].neibor[layer].begin(), rnode_c1[s_node].neibor[layer].end(), target_nonshare_node[ppi]));
					rnode_c1[target_nonshare_node[ppi]].neibor[layer].erase(rnode_c1[target_nonshare_node[ppi]].neibor[layer].begin() + ppi2);
					target_nonshare_node.erase(target_nonshare_node.begin() + ppi);
				}
				else
				{
					js++;
				}
			}
			add_neibor(rnode_c1, i, target_share_node, layer);//添加共有节点

			//rnode_c2
			delete_neibor(rnode_c2, i, target_share_node, layer);//删除共有节点
			for (int j = 0, js = 0, jj = rnode_c2[i].neibor[layer].size(); j < jj; j++)//删除非共有节点
			{
				t_node = rnode_c2[i].neibor[layer][js];
				ppi = ((double)rand() / RAND_MAX)*(target_nonshare_node2.size() - 1);
				ppi2 = ((double)rand() / RAND_MAX)*(rnode_c2[target_nonshare_node2[ppi]].neibor[layer].size() - 1);
				s_node = rnode_c2[target_nonshare_node2[ppi]].neibor[layer][ppi2];
				if ((find(rnode_c2[t_node].neibor[layer].begin(), rnode_c2[t_node].neibor[layer].end(), s_node) == rnode_c2[t_node].neibor[layer].end()) && (t_node != s_node))
				{
					rnode_c2[i].neibor[layer].erase(rnode_c2[i].neibor[layer].begin() + js);//删除非共有节点
					rnode_c2[t_node].neibor[layer].erase(find(rnode_c2[t_node].neibor[layer].begin(), rnode_c2[t_node].neibor[layer].end(), i));
					rnode_c2[t_node].neibor[layer].push_back(s_node);
					rnode_c2[s_node].neibor[layer].push_back(t_node);
					rnode_c2[i].neibor[layer].push_back(target_nonshare_node2[ppi]);//添加非共有节点
					rnode_c2[target_nonshare_node2[ppi]].neibor[layer].push_back(i);
					rnode_c2[s_node].neibor[layer].erase(find(rnode_c2[s_node].neibor[layer].begin(), rnode_c2[s_node].neibor[layer].end(), target_nonshare_node2[ppi]));
					rnode_c2[target_nonshare_node2[ppi]].neibor[layer].erase(rnode_c2[target_nonshare_node2[ppi]].neibor[layer].begin() + ppi2);
					target_nonshare_node2.erase(target_nonshare_node2.begin() + ppi);
				}
				else
				{
					js++;
				}
			}
			add_neibor(rnode_c2, i, target_share_node, layer);//添加共有节点

		}
	}
}

void find_targetnode(vector<node> rnode1, vector<node> rnode2,int nodenumber,vector<int> &target_nonshare_ndoe, vector<int> &target_share_ndoe,int layer)//找到rnode2的邻居中与rnode1的非共有邻节点
{
	target_nonshare_ndoe.clear();
	target_share_ndoe.clear();
	for (int i = 0; i < rnode1[nodenumber].neibor[layer].size(); i++)
	{
		if (find(rnode1[nodenumber].neibor[layer].begin(), rnode1[nodenumber].neibor[layer].end(), rnode2[nodenumber].neibor[layer][i]) == rnode1[nodenumber].neibor[layer].end())//node2的邻居节点在node1没有
		{
			target_nonshare_ndoe.push_back(rnode2[nodenumber].neibor[layer][i]);
		}
		else
		{
			target_share_ndoe.push_back(rnode2[nodenumber].neibor[layer][i]);
		}
	}
}

void delete_neibor(vector<node> &rnode,int nodenum, vector<int> delete_set, int layer)//删除delete_set中的节点
{

	for (int i = 0; i < delete_set.size(); i++)
	{
		rnode[nodenum].neibor[layer].erase(find(rnode[nodenum].neibor[layer].begin(), rnode[nodenum].neibor[layer].end(), delete_set[i]));
	}
}

void add_neibor(vector<node> &rnode, int nodenum, vector<int> add_set, int layer)//添加add_set中的节点
{
	for (int i = 0; i < add_set.size(); i++)
	{
		rnode[nodenum].neibor[layer].push_back(add_set[i]);
	}
}

void crossover_edge(vector<node> &rnode_p1, vector<node> &rnode_p2, vector<node> &rnode_c1, vector<node> &rnode_c2)//交换两个个体某一层的网络产生下一代
{
	rnode_c1 = rnode_p1;
	rnode_c2 = rnode_p2;
	vector<node> rnodeB(rnode_p1.size());
	for (int j = 0; j < LayerN; j++)
	{
		if ((double)rand() / RAND_MAX < Pc2)
		{
			for (int i = 0; i < rnode_c1.size(); i++)
			{
				rnodeB[i].degree[j] = rnode_c1[i].degree[j];
				rnodeB[i].neibor[j] = rnode_c1[i].neibor[j];
			}
			for (int i = 0; i < rnode_c1.size(); i++)
			{
				rnode_c1[i].degree[j] = rnode_c2[i].degree[j];
				rnode_c1[i].neibor[j] = rnode_c2[i].neibor[j];
			}
			for (int i = 0; i < rnode_c1.size(); i++)
			{
				rnode_c2[i].degree[j] = rnodeB[i].degree[j];
				rnode_c2[i].neibor[j] = rnodeB[i].neibor[j];
			}
		}
	}
}