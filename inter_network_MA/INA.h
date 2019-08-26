#ifndef _INA_H
#define _INA_H
#define alpha 0.9
#define LayerN 2
#define Pc 0.5
#define Pc2 0.3//crossover 
#define Pm 0.9
#define namuta 3
//3-0.240
//2-0.240
//swap 0.209
#include"vector"
using namespace std;
typedef struct node
{
	vector<int> degree[LayerN];//存储节点的度
	vector<int> neibor[LayerN];//存储结点的邻居
	vector<int> crruentdegree[LayerN];//存储当前的度
	vector<int> degree_record[LayerN];//记录度
	vector<double> SA;//生存能力SA
	double SA2;//生产能力2，以节点存在回合为标准
	vector<double> eSA[LayerN];//边的SA
	int activation = 0;//激活开关,1是激活
	int largest_connect_cluster = -1;//是否位于最大连接簇
}node;


double mainloop_MA(int N);//源.cpp

void initialization_ER(vector<node> &rnode, int layer);//initialization.cpp
void initialization_SF(vector<node> &rnode, int layer);
void initialization_WS(vector<node> &rnode, int layer);
void initialization_NW(vector<node> &rnode, int layer);
void initialization_degree(vector<node> &rnode);
int calculate_degree(vector<node> &rnode, int layer);
void activate_node(vector<node> &rnode);
void connect_node(vector<node> &rnode, int layer, int node1, int node2);
void initialSA(vector<node> &rnode);

double computeR(vector<node> &rnode, int begin_layer, int end_layer);//compute.cpp
int max_degree_P(vector<node> &rnode, int layer);
void copy_node(vector<node> &rnode, vector<node> &rnodeB);
void attack_NODE(vector<node> &rnode, int layer, int attackP);
int find_largestconnectcluster(vector<node> &rnode, int layer);
bool whether_search_end(vector<int> &search, vector<node> &rnode);
int Include_Nodes(int p, int layer, vector<int> &search, vector<node> &rnode);
void clear_lcc_flag(vector<node> &rnode);
void delete_node(vector<node> &rnode, int layer);
int calculate_res_node(vector<node> &rnode);
double computeR_single(vector<node> &rnode, int layer);

bool swap(vector<node> &rnode, int layer, int &falsenum);//swap.cpp

bool check_degree(vector<node> rnode);//tools.cpp
void mark_degree(vector<node> &rnode);
//void print_network(vector<node> rnode, int layer);
void choice_in(vector<vector<int>> &choice, int pop_size);
void choice_in2(vector<vector<int>> &choice, int pop_size);
int select_operator(vector<vector<int>> &choice, int pop_size, vector<double> R, vector<double> Rc);
int select_operator2(vector<vector<int>> &choice, int pop_size, vector<double> R, vector<double> Rc);
void copy_node_layer(vector<node> &rnode1, vector<node> &rnode2, int layer);

void crossover(vector<node> &rnode_p1, vector<node> &rnode_p2, vector<node> &rnode_c1, vector<node> &rnode_c2, int layer);//crossover.cpp
void find_targetnode(vector<node> rnode1, vector<node> rnode2, int nodenumber, vector<int> &target_nonshare_ndoe, vector<int> &target_share_ndoe, int layer);
void delete_neibor(vector<node> &rnode, int nodenum, vector<int> delete_set, int layer);
void add_neibor(vector<node> &rnode, int nodenum, vector<int> add_set, int layer);
void crossover_edge(vector<node> &rnode_p1, vector<node> &rnode_p2, vector<node> &rnode_c1, vector<node> &rnode_c2);

double computeKnn(vector<node> &rnode, int layer, int k);//computeKnn.cpp
int computeEkki(vector<node> &rnode, int targetnode, int k, int layer);
int computeNk(vector<node> &rnode, int k, int layer);

void print_network(vector<node> &rnode, int nettype, int layer);//print_network.cpp

void local_search_operator(vector<node>& rnode, int layer);//local_search.cpp
void calculate_eSA(vector<node> &rnode, int layer);
void iterator_operator(vector<node> &rnode, vector<node> rnode_before);
void calculate_SA2(vector<node> &rnode, int begin_layer, int end_layer);
void local_search_operator2(vector<node>& rnode, int layer);

#endif