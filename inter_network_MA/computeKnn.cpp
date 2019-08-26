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

double computeKnn(vector<node> &rnode, int layer, int k)
{
	double knn1=0;
	int Ekki = 0, Nk = 0;
	Nk=computeNk(rnode, k, layer);
	if (Nk == 0) {
		return 0;
	}
	for (int i = 0; i < rnode.size(); i++) {
		Ekki=computeEkki(rnode, i, k, layer);
		knn1 += (double)rnode[i].neibor[layer].size()*Ekki;
	}
	knn1 = knn1 / ((double)k*Nk);
	return knn1;
}

int computeEkki(vector<node> &rnode, int targetnode,int k,int layer)
{
	int E_num = 0;
	for (int i = 0; i < rnode[targetnode].neibor[layer].size(); i++) {
		if (rnode[rnode[targetnode].neibor[layer][i]].neibor[layer].size() == k) {
			E_num++;
		}
	}
	if (rnode[targetnode].neibor[layer].size() == k) {
		E_num = E_num * 2;
	}
	return E_num;
}

int computeNk(vector<node> &rnode, int k,int layer)
{
	int Nk_num = 0;
	for (int i = 0; i < rnode.size(); i++) {
		if (rnode[i].neibor[layer].size() == k) {
			Nk_num++;
		}
	}
	return Nk_num;
}