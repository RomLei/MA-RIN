#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>
#include"INA.h"
#include"numeric"
using namespace std;
void print_network(vector<node> &rnode, int nettype, int layer)//���������Ϊgml��ʽ,ֻ��һ���ļ���һ�����磡
{
	vector<vector<int>> G;
	G.resize(rnode.size(), vector<int>(rnode.size(), 0));
	ofstream PrintNetwork;
	if (nettype == 0)
	{
		PrintNetwork.open("InitialNetwork.gml");
	}
	else
	{
		PrintNetwork.open("FinalNetwork.gml");
	}
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
	PrintNetwork << "Creator \" N=" << G.size() << " layer=" << layer << " \"" << endl;
	PrintNetwork << "graph" << endl;
	PrintNetwork << "[" << endl;
	PrintNetwork << "  directed 0" << endl;
	for (int i = 0; i < G.size(); i++)
	{
		PrintNetwork << "  node" << endl;
		PrintNetwork << "  [" << endl;
		PrintNetwork << "    id " << i << endl;
		PrintNetwork << "  ]" << endl;
	}
	for (int i = 0; i < G.size(); i++)
	{
		for (int j = 0; j < G.size(); j++)
		{
			if (G[i][j] == 1)
			{
				PrintNetwork << "  edge" << endl;
				PrintNetwork << "  [" << endl;
				PrintNetwork << "    source " << i << endl;
				PrintNetwork << "    target " << j << endl;
				PrintNetwork << "  ]" << endl;
				G[j][i] = 0;
			}
		}
	}
	PrintNetwork << "]" << endl;
	PrintNetwork.close();
}