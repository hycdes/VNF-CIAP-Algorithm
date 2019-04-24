#include<iostream>
#include<iomanip>
#include<utility>
#include<algorithm>
#include<string>
#include<sstream>
#include<vector>
#include<stack>
#include<queue>
#include<map>
#include<fstream>
#include<cmath>
#include<time.h>

constexpr auto INF = 99999;
using namespace std;
enum sensitive { CPU, cache, memory }; //sensitive resource: 0,1,2
enum status { linear, branch_main, branch_sub }; //status of sfc: linear, have branches and is main-branch or sub-branch
//class: VNF, SFC, vertex, graph
class VNF {
private:
	int type;
	double price, QoSCapacity, CPUpree, calorie;
	sensitive s1, s2;
	vector<double> senDegree;
public:
	VNF(int type, double price, double QoSCapacity, sensitive s1, sensitive s2, vector<double> senDegree, double CPUpree) {
		this->type = type;
		this->price = price;
		this->calorie = 0;
		this->QoSCapacity = QoSCapacity;
		this->s1 = s1;
		this->s2 = s2;
		this->senDegree = senDegree;
		this->CPUpree = CPUpree;
	}
	int get_type() { return type; }
	double get_price() { return price; }
	double get_calorie() { return calorie; }
	double get_QoSCapacity() { return QoSCapacity; }
	sensitive get_major() { return s1; }
	sensitive get_minor() { return s2; }
	vector<double> get_senDegree() { return senDegree; }
	double get_CPUpree() { return CPUpree; }
	void update_calorie(double calorie) { this->calorie = calorie; }
};
class SFC {
private:
	int _num, _start, _end, _branchidx;
	status _status;
	double _avgCalorie;
	vector<VNF> VNFList;
public:
	SFC(int num, int start, int end) {
		_num = num;
		_start = start;
		_end = end;
		_status = linear;
		_branchidx = -1; // -1 means linear
		_avgCalorie = 0;
	}
	int get_num() { return _num; }
	int get_startNum() { return _start; }
	int get_endNum() { return _end; }
	int get_branchidx() { return _branchidx; }
	vector<VNF> get_VNFList() { return VNFList; }
	double get_avgCalorie() { return _avgCalorie; }
	bool is_sub() {
		if (_status == branch_sub) return true;
		else return false;
	}
	void add_VNF(VNF vnf) { VNFList.push_back(vnf); }
	void update_status(status new_status) { _status = new_status; }
	void update_branchidx(int idx) { _branchidx = idx; }
	void update_avgCalorie(double calorie) { _avgCalorie = calorie; }
	void update_vnfCalorie(int vnfPosition, double calorie) { this->VNFList[vnfPosition].update_calorie(calorie); }
	void show_detail() {
		if (_status == linear || _status == branch_main) cout << "sfc_" << _num << "\t";
		else {
			for (int i = 0; i < _branchidx + 1; i++) cout << "\t";
			cout << "   ->\t";
		}
		for (vector<VNF>::iterator it = VNFList.begin(); it != VNFList.end(); it++) {
			cout << it->get_type() << " [" << it->get_QoSCapacity() << "]\t";
		}
		cout << endl;
	}
};
class vertex {
private:
	int _num;
	double _capacity;
	bool _saturated;
	vector<int> adjacentList;
	vector<VNF> VNFPlace;
public:
	vertex(int vertexNum, vector<int> adjacentList, double capacity) {	//init vertex
		this->_num = vertexNum;
		this->_capacity = capacity;
		this->_saturated = false;
		this->adjacentList = adjacentList;
	}
	int get_Num() { return _num; }
	double get_capacity() { return _capacity; }
	bool is_pure() { return VNFPlace.empty(); }
	bool is_saturated() { return _saturated; }
	vector<int> get_AdjList() { return adjacentList; }
	vector<VNF> get_VNFPlacement() { return VNFPlace; }
	void add_Adjacent(vertex v) { this->adjacentList.push_back(v.get_Num()); }
	void add_VNF(VNF vnf) { VNFPlace.push_back(vnf); }
	void update_saturated() { _saturated = true; }
};
class graph {
private:
	vector<vertex> vertexSet;
public:
	graph() {}	//empty graph
	int get_size() { return vertexSet.size(); }
	vector<vertex>::iterator get_vertexIt(int vertexNum) {
		vector<vertex>::iterator vertexIterator;
		for (vector<vertex>::iterator it = vertexSet.begin(); it != vertexSet.end(); it++) {
			if (it->get_Num() == vertexNum) {
				vertexIterator = it;
				break;
			}
		}
		return vertexIterator;
	}
	void add_vertex(vertex v) { this->vertexSet.push_back(v); }
	void show_detail() {
		cout << "[Graph info]\ttotal vertices=" << vertexSet.size() << endl;
		for (vector<vertex>::iterator it = vertexSet.begin(); it != vertexSet.end(); it++) {
			cout << "v_" << it->get_Num() << "\tcapacity=" << it->get_capacity() << "\tadjacent list: ";
			vector<int> cur_adjList = it->get_AdjList();
			for (vector<int>::iterator it_ad = cur_adjList.begin(); it_ad != cur_adjList.end(); it_ad++)
				cout << *it_ad << " ";
			cout << endl;
		}
	}
};

template<typename T> T toNum(string str) {
	T num;
	stringstream ss;
	ss << str;
	ss >> num;
	return num;
}
template<typename T> string toStr(T num) {
	stringstream ss;
	ss << num;
	return ss.str();
}
double vectorPlus(vector<double> v1, vector<double> v2) { return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]; }
void verticesMerge(vector<vertex> &set_vertices, vector<int> list_num, graph GlobalMap) { // 将点号表v2对应的点合并进点集v1中，注意v2和v1数据类型不同
	for (vector<int>::iterator it = list_num.begin(); it != list_num.end(); it++) {
		set_vertices.push_back(*GlobalMap.get_vertexIt(*it));
	}
}
bool isInVertexSet(int vertexNum, vector<vertex> vertexSet) {		//判断一个顶点号对应顶点是否在某个点集中
	for (vector<vertex>::iterator it = vertexSet.begin(); it != vertexSet.end(); it++)
		if (it->get_Num() == vertexNum) return true;
	return false;
}
bool isDeadend(graph G, int vertexNum, vector<vertex> pathSet, vector<vertex> deadendSet) {	//末路：无路可走，即邻接表全是前路或末路点
	vector<int> adjList = (G.get_vertexIt(vertexNum))->get_AdjList();
	for (vector<int>::iterator it = adjList.begin(); it != adjList.end(); it++) {
		int vertexNum = *it;
		bool pass = false;
		for (vector<vertex>::iterator itPath = pathSet.begin(); itPath != pathSet.end(); itPath++) {	//是否为前路
			if (vertexNum == itPath->get_Num()) {
				pass = true;
				break;
			}
		}
		if (pass) continue;
		else {
			for (vector<vertex>::iterator itDE = deadendSet.begin(); itDE != deadendSet.end(); itDE++) {	//是否为末路
				if (vertexNum == itDE->get_Num()) {
					pass = true;
					break;
				}
			}
		}
		if (!pass) return false;	//只要有邻点两者皆非则不是末路
	}
	return true;
}
vector<vector<int>> integerComposition(int parts_max, int size_max, int sum) {
	vector<vector<int>> compositionSet;
	vector<int> temp_solu;
	for (int i = 0; i < parts_max; i++) temp_solu.push_back(0);
	if (parts_max*size_max >= sum) { // possible to make composition, else return empty
		int idx_pick = 0, total = 0;
		bool is_gonext = false;
		while (idx_pick >= 0) {
			if (total < sum) {
				if (temp_solu[idx_pick] < size_max) { // can add at this position
					temp_solu[idx_pick]++;
					total++;
					if (total == sum) {
						compositionSet.push_back(temp_solu);
						is_gonext = true;
					}
				}
				else {
					is_gonext = true;
				}
			}
			if (is_gonext) { // need go next
				if (idx_pick + 1 < parts_max) { // can go next
					if (total == sum) { // find next solution
						temp_solu[idx_pick]--;
						total--;
					}
					idx_pick++; // go next
				}
				else { // backtrack: meet border
					total -= temp_solu[idx_pick];
					temp_solu[idx_pick] = 0; // clear this position
					idx_pick--;
					while (true) { // find last non-zero, if all done, idx_pick=-1
						if (temp_solu[idx_pick] <= 0) idx_pick--;
						if (idx_pick < 0) break;
						if (temp_solu[idx_pick] > 0) break;
					}
					if (idx_pick >= 0) {
						temp_solu[idx_pick]--;
						total--;
						idx_pick++;
					}
				}
				is_gonext = false;
			}
		}
	}
	return compositionSet;
}
bool cmp_placement(pair<double, vector<int>> p1, pair<double, vector<int>> p2) { return p1.first < p2.first; } // compare function to sort placement
bool cmp_sfc(SFC s1, SFC s2) { return s1.get_avgCalorie() > s2.get_avgCalorie(); }
void sort_sfc(vector<SFC> &S) {
	vector<SFC> mainchainSet;
	map<int, SFC> subchainSet;
	for (vector<SFC>::iterator it = S.begin(); it != S.end(); it++) {
		if (it->is_sub()) subchainSet.insert(make_pair(it->get_num(), *it));
		else mainchainSet.push_back(*it);
	}
	sort(mainchainSet.begin(), mainchainSet.end(), cmp_sfc);
	vector<SFC>::iterator it = mainchainSet.begin();
	while (it != mainchainSet.end()) {
		int sfc_num = it->get_num();
		if (subchainSet.find(sfc_num) != subchainSet.end()) { // have subchain
			mainchainSet.insert(it + 1, subchainSet.find(sfc_num)->second);
			subchainSet.erase(sfc_num);
			it = mainchainSet.begin();
		}
		else it++;
	}
	S = mainchainSet;
}

//init 初始化函数
void init_getdata(string datapath_vnf, map<int, VNF> &vnfLib, string datapath_sfc, vector<SFC> &S, string datapath_graph, graph &G) {
	string str_oneline;
	ifstream vnf_data, vertex_data, sfc_data;
	//录入VNF库
	vnf_data.open(datapath_vnf);
	while (!vnf_data.eof()) {
		getline(vnf_data, str_oneline);
		istringstream istr_oneline(str_oneline);
		string contents[9];
		int temp_type;
		double temp_price, temp_capacity, temp_pree;
		sensitive temp_major, temp_minor;
		vector<double> temp_sen;
		istr_oneline >> contents[0] >> contents[1] >> contents[2] >> contents[3] >> contents[4] >> contents[5] >> contents[6] >> contents[7] >> contents[8];
		temp_type = toNum<int>(contents[0]);
		temp_price = toNum<double>(contents[1]);
		temp_capacity = toNum<double>(contents[2]);
		if (toNum<int>(contents[3]) == 0) temp_major = CPU;
		else if (toNum<int>(contents[3]) == 1) temp_major = cache;
		else if (toNum<int>(contents[3]) == 2) temp_major = memory;
		if (toNum<int>(contents[4]) == 0) temp_minor = CPU;
		else if (toNum<int>(contents[4]) == 1) temp_minor = cache;
		else if (toNum<int>(contents[4]) == 2) temp_minor = memory;
		temp_sen.push_back(toNum<double>(contents[5])); temp_sen.push_back(toNum<double>(contents[6])); temp_sen.push_back(toNum<double>(contents[7]));
		temp_pree = toNum<double>(contents[8]);
		VNF temp_vnf(temp_type, temp_price, temp_capacity, temp_major, temp_minor, temp_sen, temp_pree);
		vnfLib.insert(make_pair(temp_type, temp_vnf));
	}
	vnf_data.close();
	//录入SFC信息
	sfc_data.open(datapath_sfc);
	while (!sfc_data.eof()) {
		getline(sfc_data, str_oneline);
		bool is_sub = false;
		istringstream istr_oneline(str_oneline);
		int temp_num, temp_start, temp_end, temp_type, vnfidx = 0;
		double temp_capacity;
		string str_content;
		istr_oneline >> str_content;
		temp_num = toNum<int>(str_content);
		istr_oneline >> str_content;
		temp_start = toNum<int>(str_content);
		istr_oneline >> str_content;
		temp_end = toNum<int>(str_content);
		SFC temp_s(temp_num, temp_start, temp_end);
		while (istr_oneline >> str_content) {
			int temp_value = toNum<int>(str_content);
			if (temp_value < 0) {
				is_sub = true;
				istr_oneline >> str_content;
			}
			else {
				temp_type = temp_value;
				istr_oneline >> str_content;
				temp_capacity = toNum<double>(str_content);
				VNF vnfInfo = vnfLib.find(temp_type)->second;
				VNF temp_vnf(temp_type, vnfInfo.get_price(), temp_capacity, vnfInfo.get_major(), vnfInfo.get_minor(), vnfInfo.get_senDegree(), vnfInfo.get_CPUpree());
				temp_s.add_VNF(temp_vnf);
				if (is_sub) { // mark sub-chain
					is_sub = false;
					S.back().update_branchidx(vnfidx - 1);
					S.back().update_status(branch_main);
					temp_s.update_branchidx(vnfidx - 1);
					temp_s.update_status(branch_sub);
				}
			}
			vnfidx++;
		}
		S.push_back(temp_s);
	}
	sfc_data.close();
	//录入图信息
	vertex_data.open(datapath_graph);
	while (!vertex_data.eof()) {
		getline(vertex_data, str_oneline);
		istringstream istr_oneline(str_oneline);
		vector<int> temp_adj;
		int temp_num, temp_capacity;
		string str_content;
		istr_oneline >> str_content;
		temp_num = toNum<int>(str_content);
		istr_oneline >> str_content;
		temp_capacity = toNum<int>(str_content);
		while (istr_oneline >> str_content) {
			temp_adj.push_back(toNum<int>(str_content));
		}
		vertex temp_v(temp_num, temp_adj, temp_capacity);
		G.add_vertex(temp_v);
	}
	vertex_data.close();
}
void inti_getcalorie(vector<SFC> &S, map<int, double> &vnfSet) {	//SFC平均热值
	int totalNum = 0;
	double totalCapacity = 0;
	for (vector<SFC>::iterator itS = S.begin(); itS != S.end(); itS++) {		//count all VNF
		vector<VNF> serviceChain = itS->get_VNFList();
		for (vector<VNF>::iterator itV = serviceChain.begin(); itV != serviceChain.end(); itV++) {
			if (vnfSet.find(itV->get_type()) == vnfSet.end())
				vnfSet.insert(make_pair(itV->get_type(), itV->get_QoSCapacity()));
			else
				vnfSet.find(itV->get_type())->second += (itV->get_QoSCapacity());
			totalCapacity += itV->get_QoSCapacity();
		}
	}
	for (vector<SFC>::iterator itS = S.begin(); itS != S.end(); itS++) {		//update calorie of VNF and SFC
		int sfcLen = (itS->get_VNFList()).size();
		double sfcTotalCalorie = 0;
		vector<VNF> serviceChain = itS->get_VNFList();
		for (int i = 0; i < sfcLen; i++) {		//update each vnfCalorie in a SFC
			double price = serviceChain[i].get_price();
			double heat = (serviceChain[i].get_QoSCapacity()) / (vnfSet.find(serviceChain[i].get_type())->second);
			double calorie = price * heat;
			sfcTotalCalorie += calorie;
			itS->update_vnfCalorie(i, calorie);
		}
		double temp_result = sfcTotalCalorie / double(sfcLen);
		itS->update_avgCalorie(temp_result);
	}
}

// debug ONLY
void show_vertices(vector<vertex> vSet) {
	cout << endl << "V_SET: ";
	for (vector<vertex>::iterator it = vSet.begin(); it != vSet.end(); it++) cout << it->get_Num() << ", ";
	cout << endl;
}
void show_map(map<int, int> M) {
	cout << "Map value:";
	for (map<int, int>::iterator it = M.begin(); it != M.end(); it++) {
		cout << "<" << it->first << "," << it->second << "> ";
	}
	cout << endl;
}
template<typename T> void show_vector(vector<T> v) {
	cout << endl << "VECTOR: ";
	int size = v.size();
	for (int i = 0; i < size;i++) cout << v[i] << ", ";
	cout << endl;
}
void show_icSet(vector<vector<int>> compositionSet) {
	for (vector<vector<int>>::iterator it = compositionSet.begin(); it != compositionSet.end(); it++) {
		for (vector<int>::iterator it_e = it->begin(); it_e != it->end(); it_e++) {
			cout << *it_e << " ";
		}
		cout << endl;
	}
}
void show_ic(vector<int> composition) {
	cout << "ic: ";
	vector<int>::iterator it = composition.begin();
	cout << *it;
	it++;
	for (; it != composition.end(); it++) cout << " -> " << *it;
	cout << endl;
}
void show_icSetPair(vector<pair<double, vector<int>>> set_placement) {
	for (vector<pair<double, vector<int>>>::iterator it = set_placement.begin(); it != set_placement.end(); it++) {
		cout << "value=" << it->first << "\t";
		show_ic(it->second);
	}
	cout << endl;
}
void show_input(map<int, double> vnfSet_capacity, map<int, VNF> vnfLib, vector<SFC> S, graph G) {
	cout << "[VNF info]\ttotal vnf=" << vnfLib.size() << endl;
	for (map<int, VNF>::iterator it = vnfLib.begin(); it != vnfLib.end(); it++)
		cout << "vnf_" << (it->second).get_type() << "\tprice=" << (it->second).get_price() << "\tI=(" << ((it->second).get_senDegree())[0] << ", " << ((it->second).get_senDegree())[1] << ", " << ((it->second).get_senDegree())[2] << ")\tCPUpree=" << (it->second).get_CPUpree() << endl;
	cout << endl;

	cout << "[SFC info]\ttotal sfc=" << S.back().get_num() + 1 << endl;
	for (vector<SFC>::iterator it = S.begin(); it != S.end(); it++)
		it->show_detail();
	cout << endl;
	G.show_detail();
	cout << endl;

	cout << "[VNF-capacity]" << endl;
	for (map<int, double>::iterator it = vnfSet_capacity.begin(); it != vnfSet_capacity.end(); it++) {
		cout << "vnf_" << it->first << "=" << it->second << "\tprice=" << (vnfLib.find(it->first)->second).get_price() << endl;
	}
	cout << endl;

	cout << "[SFC-calorie]" << endl;
	for (vector<SFC>::iterator it = S.begin(); it != S.end(); it++) {
		cout << "sfc_" << it->get_num();
		if (it->is_sub()) cout << "'";
		cout << "\t";
		vector<VNF> current_vnfList = it->get_VNFList();
		for (vector<VNF>::iterator itv = current_vnfList.begin(); itv != current_vnfList.end(); itv++) {
			cout << "vnf_" << itv->get_type() << "=" << itv->get_calorie() << "  ";
		}
		cout << endl << "\tavgCalorie=" << it->get_avgCalorie() << endl;
	}
	cout << endl;
}
void show_path(vector<vertex> path) {
	for (vector<vertex>::iterator it = path.begin(); it != path.end(); it++) {
		if (it != path.begin()) cout << " -> ";
		cout << it->get_Num();
	}
	cout << endl;
}