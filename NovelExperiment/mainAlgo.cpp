#include"graphADT.h"
/*
< ADT > // need update!
	VNF		get_type/price/calorie/capacity/major/minor/senDegree/CPUpree
			is_proper
			update_calorie/capacity
	SFC		get_num/startNum/endNum/VNFList/avgCalorie
			add_VNF
			update_substatus/branchidx/avgCalorie/vnfCalorie
	vertex	get_Num/AdjList/capacity/VNFPlacement/minCalorie/totalCalorie
			add_Adjacent
	graph	add_vertex
			show_detail
< Basic function >
	toNum<num_type>		>num_type	string->num
	toStr				>string		num->string
	vectorPlus			>double		vector-plus
	verticesMerge		将点号向量的对应点合并进点集
	isInVertexSet		判断点是否在点集中
	isDeadend			判断寻路中的点是否末路(判断其邻接表点是否全是前路或末路)
	integerComposition	返回所有可能的整数组合
< Init function >
	init_getdata		读取数据文件
	inti_getcalorie		SFC集合中所有SFC平均及所有其内VNF热值计算
*/

/*
【算法】
	get_declineRatio	VNF衰退比计算
		get_senSitu		敏感情况，根据I导出e
		get_cpuPree		CPU竞争衰退
		get_declineNetIO	网络I/O竞争衰退

	get_similarity		>double		计算链和点相似度

	find_path			>vector<vertex>	寻路，返回点向量，或空集
		get_similarity			>double		计算两个vector<VNF>余弦相似度
		get_cutoff				>int		找近路点，返回点编号或-1
		get_nextVertex			>int		找下一步，返回点编号或-1

	deploy_sfc			>string			部署sfc，返回部署信息
		cmp_placement			>bool		比较放置方案优劣的比较函数
		get_actualRate			>double		计算实际服务率
		get_declineCalorie		>double		计算损失的热值
		is_fit					>bool		放置是否合适
*/

vector<double> get_senSitu(vector<double> senDegree) {	//敏感情况计算
	vector<double> e;
	for (int i = 0; i < 3; i++) {
		if (senDegree[i] > 0) e.push_back(1);
		else e.push_back(0);
	}
	return e;
}
double get_cpuPree(VNF VNF1, VNF VNF2) {		//CPU竞争衰退：C(Nk,Nl)
	double Ckl;
	bool complementary = false;
	vector<double> e1 = get_senSitu(VNF1.get_senDegree());
	vector<double> e2 = get_senSitu(VNF2.get_senDegree());
	if (e1[0] > 0 && e1[1] > 0 && e2[0] > 0 && e2[1] > 0) complementary = true;
	if (complementary) {		//g(Nk,Nl)
		if (VNF1.get_CPUpree() > VNF2.get_CPUpree())
			Ckl = 0.02*(30 * VNF1.get_CPUpree() + 1 / VNF2.get_CPUpree());
		else
			Ckl = 0.02*(VNF1.get_CPUpree() + 9 * VNF2.get_CPUpree());
	}
	else {						//f(Nk,Nl)
		if (VNF1.get_CPUpree() > VNF2.get_CPUpree()) Ckl = VNF2.get_CPUpree()/VNF1.get_CPUpree();
		else Ckl = VNF1.get_CPUpree()*VNF2.get_CPUpree();
	}
	return Ckl;
}
double get_declineNetIO(VNF VNF1, VNF VNF2) {	//网络I/O带宽竞争衰退：D(Nk,Nl), x=0.157, y=0.06, z=0.145, m,n=0.02
	double Dkl;
	bool complementary = false;
	vector<double> e1 = get_senSitu(VNF1.get_senDegree());
	vector<double> e2 = get_senSitu(VNF2.get_senDegree());
	if (e1[0] > 0 && e1[1] > 0 && e2[0] > 0 && e2[1] > 0) complementary = true;
	if (complementary) {		//(if >) z*f(Nk,Nl)+x*g(Nk,Nl)		if(>=):z*f(Nk,Nl)-y*g(Nk,Nl)
		if (VNF1.get_CPUpree() > VNF2.get_CPUpree()) Dkl = 0.157*0.02*(30 * VNF1.get_CPUpree() + 1 / VNF2.get_CPUpree()) + 0.145*VNF2.get_CPUpree() / VNF1.get_CPUpree();
		else Dkl = 0.145*VNF1.get_CPUpree()*VNF2.get_CPUpree() - 0.06*0.02*(VNF1.get_CPUpree() + 9 * VNF2.get_CPUpree());
	}
	else {						//Δ*f(Nk,Nl)
		if (VNF1.get_CPUpree() > VNF2.get_CPUpree()) Dkl = 0.2*VNF2.get_CPUpree()/VNF1.get_CPUpree();
		else Dkl = 0.2*VNF1.get_CPUpree()*VNF2.get_CPUpree();
	}
	return Dkl;
}
double get_declineRatio(VNF VNF1, VNF VNF2, bool isG) {		//get delta
	vector<double> e1 = get_senSitu(VNF1.get_senDegree());
	vector<double> e2 = get_senSitu(VNF2.get_senDegree());
	vector<double> I1 = VNF1.get_senDegree();
	vector<double> I2 = VNF2.get_senDegree();
	double delta = 0;
	if (VNF1.get_type() == VNF2.get_type())
		delta = 0;
	else {
		if (isG) delta = vectorPlus(I1, e2)*vectorPlus(I1, I2)*get_cpuPree(VNF1, VNF2);
		else delta = vectorPlus(I1, e2)*vectorPlus(I1, I2)*get_cpuPree(VNF1, VNF2) + get_declineNetIO(VNF1, VNF2);
	}
	return delta;
}
double get_similarity(vector<VNF> vlist1,vector<VNF> vlist2) {
	double plus = 0, square = 1, square_1 = 0, square_2 = 0;
	map<int, double> cal1, cal2;
	vector<int> delete1, delete2;
	for (vector<VNF>::iterator it = vlist1.begin(); it != vlist1.end(); it++) { // statistic cal1
		int vnf_type = it->get_type();
		double vnf_calorie = it->get_calorie();
		if (cal1.find(vnf_type) == cal1.end()) cal1.insert(make_pair(vnf_type, vnf_calorie));
		else cal1.find(vnf_type)->second += vnf_calorie;
	}
	for (vector<VNF>::iterator it = vlist2.begin(); it != vlist2.end(); it++) { // statistic cal1
		int vnf_type = it->get_type();
		double vnf_calorie = it->get_calorie();
		if (cal2.find(vnf_type) == cal2.end()) cal2.insert(make_pair(vnf_type, vnf_calorie));
		else cal2.find(vnf_type)->second += vnf_calorie;
	}
	for (map<int, double>::iterator it = cal1.begin(); it != cal1.end(); it++) { // record independent of cal1
		if (cal2.find(it->first) == cal2.end()) delete1.push_back(it->first);
	}
	for (map<int, double>::iterator it = cal2.begin(); it != cal2.end(); it++) { // record independent of cal2
		if (cal1.find(it->first) == cal1.end()) delete2.push_back(it->first);
	}
	for (vector<int>::iterator it = delete1.begin(); it != delete1.end(); it++)
		cal1.erase(*it);
	for (vector<int>::iterator it = delete2.begin(); it != delete2.end(); it++)
		cal2.erase(*it);
	for (map<int, double>::iterator it = cal1.begin(); it != cal1.end(); it++) // vertical
		plus += (it->second)*(cal2.find(it->first)->second);
	for (map<int, double>::iterator it = cal1.begin(); it != cal1.end(); it++) // horizontal
		square_1 += pow(it->second, 2);
	for (map<int, double>::iterator it = cal2.begin(); it != cal2.end(); it++)
		square_2 += pow(it->second, 2);
	square = sqrt(square_1*square_2);
	return plus / square;
}

int get_cutoff(graph GlobalMap, int startNum, int endNum, vector<vertex> pathSet, vector<vertex> deadendSet, vector<vertex> detourSet) {		//返回近路点编号，或-1代表无法抵达终点
	map<int, int> edge;			//edge：		结界点集 <num, distance>
	map<int, int> discover;		//discover：	探索点集 <num, distance>, d_start=0, d_blockade=-1
	map<int, int> detour;		//detour:	绕路点集 <num, distance>, 第一步不能走该集合点
	map<int, stack<int>> trace;	//trace:	路径栈
	bool is_firststep = true;	//是否为第一步
	int currentNum = -1;			//当前点
	int currentDistance = 0;			//当前点起点距
	vector<int> currentAdjList;		//当前点邻接表，上述三个都在遍历结界刷新
	// 结界找最近点 -> 非终点则更新邻接点 -> 更新结界
	edge.insert(make_pair(startNum, 0));	//录入：第一次结界为起点
	stack<int> tr_start;
	tr_start.push(-1);
	trace.insert(make_pair(startNum, tr_start));
	for (vector<vertex>::iterator it = pathSet.begin(); it != pathSet.end(); it++) { // discover path
		discover.insert(make_pair(it->get_Num(), -1));
	}
	for (vector<vertex>::iterator it = deadendSet.begin(); it != deadendSet.end(); it++) { // discover deadend
		discover.insert(make_pair(it->get_Num(), -1));
	}
	for (vector<vertex>::iterator it = detourSet.begin(); it != detourSet.end(); it++) { // include detour vertices
		detour.insert(make_pair(it->get_Num(), -1));
	}
	while (!edge.empty()) {																				//寻找近路点：结界消失（可能末路），或终点为结界最近点
		int min_edgeNum = -1, min_edgeDistance = INF;
		for (map<int, int>::iterator it = edge.begin(); it != edge.end(); it++) {						//遍历结界
			if (it->second < min_edgeDistance || (it->second == min_edgeDistance && it->first == endNum)) {	//找最近点，等距若为终点则选取
				min_edgeNum = it->first;
				min_edgeDistance = it->second;
			}
		}
		if (min_edgeNum == endNum) {		//如果抵达终点
			discover.insert(make_pair(min_edgeNum, min_edgeDistance));
			stack<int> edgeTrace = trace.find(currentNum)->second;
			edgeTrace.push(currentNum);
			trace.insert(make_pair(min_edgeNum, edgeTrace));
			break;
		}
		else {		//尚未抵达终点
			//cout << "edge:"; show_map(edge);
			//cout << "discover:"; show_map(discover);
			//cout << "min=" << min_edgeNum << endl;
			currentNum = min_edgeNum;
			currentDistance = min_edgeDistance;	//更新当前点、当前权重
			discover.insert(make_pair(currentNum, currentDistance));	//探索该点s
			edge.erase(currentNum);		//结界中清除该点
			currentAdjList = (GlobalMap.get_vertexIt(currentNum))->get_AdjList();	//获取邻接表
			for (vector<int>::iterator it = currentAdjList.begin(); it != currentAdjList.end(); it++) {
				if (discover.find(*it) == discover.end()) {		//非探索点
					if (edge.find(*it) == edge.end()) {		//非结界点
						if (is_firststep && detour.find(*it) != detour.end()) {
							is_firststep = false;
							edge.insert(make_pair(*it, currentDistance + INF)); // first step not in detour
						}
						else {
							edge.insert(make_pair(*it, currentDistance + 1));
							stack<int> edgeTrace = trace.find(currentNum)->second;
							edgeTrace.push(currentNum);
							trace.insert(make_pair(*it, edgeTrace));
						}
					}
					else if (edge.find(*it)->second > currentDistance + 1) {	//结界内可以更新
						edge.find(*it)->second = currentDistance + 1;
						stack<int> edgeTrace = trace.find(currentNum)->second;
						edgeTrace.push(currentNum);
						trace.insert(make_pair(*it, edgeTrace));
					}
				}
			}
		}
	}
	if (discover.find(endNum) != discover.end()) { // endnum have been discoverd
		int a = 1;
		stack<int> path;
		if (trace.find(endNum) != trace.end()) { // endnum reachable
			path = trace.find(endNum)->second;
		}
		else {
			return -1; // endnum in deadend
		}
		int cutoffNum = endNum;
		while (true) {
			if (path.top() != startNum) {
				cutoffNum = path.top();
				path.pop();
			}
			else {
				break;
			}
		}
		return cutoffNum;
	}
	else { // endnum not discovered, out of area
		return -1;
	}
}
int get_passby(graph GlobalMap, int startNum, vector<vertex> pathSet, vector<vertex> deadendSet) {
	vector<vertex> detourSet;
	vector<int> currentAdjList = GlobalMap.get_vertexIt(startNum)->get_AdjList();
	map<int, stack<int>> route;
	map<int, int> discover;
	int currentNum, passby = -1;
	stack<int> currentTrace, newTrace;
	queue<int> edge;
	edge.push(startNum);
	currentTrace.push(-1);
	route.insert(make_pair(startNum, currentTrace));
	for (vector<vertex>::iterator it = pathSet.begin(); it != pathSet.end(); it++) { // discover path
		discover.insert(make_pair(it->get_Num(), -1));
	}
	for (vector<vertex>::iterator it = deadendSet.begin(); it != deadendSet.end(); it++) { // discover deadend
		discover.insert(make_pair(it->get_Num(), -1));
	}
	while (true) {
		currentNum = edge.front();
		currentTrace = route.find(currentNum)->second;
		currentAdjList = GlobalMap.get_vertexIt(currentNum)->get_AdjList();
		if (currentNum != startNum && GlobalMap.get_vertexIt(currentNum)->is_pure()) { // find pure vertex
			//passby = get_cutoff(GlobalMap, startNum, currentNum, pathSet, deadendSet, detourSet);
			int this_num = currentTrace.top(), last_num = -1;
			while (this_num != startNum) {
				currentTrace.pop();
				last_num = this_num;
				this_num = currentTrace.top();
			}
			passby = last_num;
			break;
		}
		else {
			for (vector<int>::iterator it = currentAdjList.begin(); it != currentAdjList.end(); it++) {
				if (route.find(*it) == route.end() && discover.find(*it) == discover.end()) {
					edge.push(*it);
					newTrace = currentTrace;
					newTrace.push(currentNum);
					route.insert(make_pair(*it, newTrace));
				}
			}
			edge.pop();
		}
		if (edge.empty()) break;
	}
	return passby;
}
int get_nextVertex(graph GlobalMap, vector<VNF> list_request, int startNum, int endNum, vector<vertex> pathSet, vector<vertex> deadendSet) {	//return vertex_num, no way = -1
	// 当前点邻接状态所属：共享集（相似度>0）、纯净集（空节点）、平行集（相似度=0）、饱和集（饱和路径点）
	// 优先度：共享集 > 纯净集 > 平行集 > 饱和集，等价选取近路点，无路可走回溯
	cout << startNum << " ";
	int nextNum = -1;
	double max_similarity = 0;
	vector<int> currentAdjList = (GlobalMap.get_vertexIt(startNum))->get_AdjList();
	vector<int> set_share, set_pure, set_parallel, set_saturated;
	vector<vertex> detourSet;
	for (vector<int>::iterator it = currentAdjList.begin(); it != currentAdjList.end(); it++) {		//对邻接表所有点
		vertex current_V = *GlobalMap.get_vertexIt(*it);
		vector<VNF> list_placement = current_V.get_VNFPlacement();
		if (!isInVertexSet(*it, pathSet) && !isInVertexSet(*it, deadendSet)) {	//此点可走：不是前路或末路s
			double similarity = get_similarity(list_request, list_placement);
			if (similarity > 0 && !current_V.is_saturated()) {
				if (similarity > max_similarity) {
					max_similarity = similarity;
					set_share.clear();
					set_share.push_back(*it); // set_share: vertex with VNF can be shared
				}
				else if (similarity == max_similarity) {
					set_share.push_back(*it);
				}
			}
			else {
				if (list_placement.size() == 0) set_pure.push_back(*it); // set_pure: empty vertices
				else {
					if (current_V.is_saturated()) set_saturated.push_back(*it); // set_saturated: saturated vertices cannot deploy again, passby
					else set_parallel.push_back(*it); // set_parallel: simi=0 but not empty, lowest priority
				}
			}
		}
	}
	if (!set_share.empty()) { // #1 priority: able to share
		if (set_share.size() == 1) nextNum = set_share.front();
		else {
			verticesMerge(detourSet, set_pure, GlobalMap);
			verticesMerge(detourSet, set_parallel, GlobalMap);
			verticesMerge(detourSet, set_saturated, GlobalMap);
			nextNum = get_cutoff(GlobalMap, startNum, endNum, pathSet, deadendSet, detourSet);
		}
	}
	else if (!set_pure.empty()) { // #2 priority: select empty vertex
		if (set_pure.size() == 1) nextNum = set_pure.front();
		else {
			verticesMerge(detourSet, set_parallel, GlobalMap);
			verticesMerge(detourSet, set_saturated, GlobalMap);
			nextNum = get_cutoff(GlobalMap, startNum, endNum, pathSet, deadendSet, detourSet);
		}
	}
	else if (!set_parallel.empty()) { // #3 priority: select simi=0 vertex (parallel)
		if (set_parallel.size() == 1) nextNum = set_parallel.front();
		else {
			verticesMerge(detourSet, set_saturated, GlobalMap);
			nextNum = get_cutoff(GlobalMap, startNum, endNum, pathSet, deadendSet, detourSet);
		}
	}
	else { // #4 priority(last selection): passby a saturated vertex, and pick which closest to a nearby pure vertex
		if (set_saturated.size() == 1) nextNum = set_saturated.front();
		else {
			int passby = get_passby(GlobalMap, startNum, pathSet, deadendSet);
			if (passby > 0) nextNum = passby;
			else nextNum = get_cutoff(GlobalMap, startNum, endNum, pathSet, deadendSet, detourSet);
		}
	}
	return nextNum;
}
vector<vertex> find_path(graph GlobalMap, SFC s, int sub_startpoint, vector<vertex> frontSet) {		//寻路，返回路径点集pathSet，空集代表无路可走
	vector<vertex> pathSet;
	vector<vertex> deadendSet = frontSet;
	int startNum = s.get_startNum();
	if (sub_startpoint > 0) startNum = sub_startpoint;
	int endNum = s.get_endNum();
	vector<VNF> list_request = s.get_VNFList();
	int currentNum = startNum;
	while (currentNum != endNum && currentNum != -1) {		//找到终点或者判定死局
		int nextNum = get_nextVertex(GlobalMap, list_request, currentNum, endNum, pathSet, deadendSet);	//寻找下一点get_nextVertex
		if (nextNum == -1) {								//无路可走则回溯
			deadendSet.push_back(*GlobalMap.get_vertexIt(currentNum));
			currentNum = pathSet.back().get_Num();
			pathSet.pop_back();
			if (pathSet.empty() && isDeadend(GlobalMap, currentNum, pathSet, deadendSet))		//死局：无路可走且无路可退
				currentNum = -1;
		}
		else {
			pathSet.push_back(*GlobalMap.get_vertexIt(currentNum));
			currentNum = nextNum;

		}
	}
	if (currentNum != -1)
		pathSet.push_back(*GlobalMap.get_vertexIt(endNum));
	return pathSet;
}
double get_actualRate(vector<VNF> list_colocated, vector<VNF>::iterator it_v, string situation) {
	bool isG = false;
	double actual_rate = 1;
	if (situation == "isG") isG = true;
	for (vector<VNF>::iterator it_co = list_colocated.begin(); it_co != list_colocated.end(); it_co++) {
		if (it_v->get_type() == it_co->get_type()) continue; // same vnf shared, no decline
		else {
			actual_rate *= 1 - get_declineRatio(*it_v, *it_co, isG);
		}
	}
	return actual_rate;
}
double get_declineCalorie(vector<vertex> path, SFC s, vector<int> composition) { // return decline calorie, -1 means invalid placement composition
	double decline_calorie = 0;
	vector<VNF> list_vnf = s.get_VNFList();
	vector<VNF>::iterator it_s = list_vnf.begin();
	vector<vertex>::iterator it_p = path.begin();
	vector<int>::iterator it_c = composition.begin();
	while (it_c != composition.end()) {
		if (*it_c > 0) {
			if (!(*it_p).is_saturated()) { // this vertex is NOT saturated
				double total_capacity = 0, all_capacity = it_p->get_capacity();
				double actual_rate = 1;
				vector<VNF> list_colocated= it_p->get_VNFPlacement(); // union vnfs in vertex
				for (int n = 0; n < *it_c; n++) { // union vnfs in this composition part of sfc
					list_colocated.push_back(*it_s);
					it_s++;
				}
				for (vector<VNF>::iterator it_v = list_colocated.begin(); it_v != list_colocated.end(); it_v++) { // calculate capacity and decline calorie of each co-located vnf
					actual_rate = get_actualRate(list_colocated, it_v, "isG");
					decline_calorie += (1 - actual_rate)*(it_v->get_calorie());
					total_capacity += (it_v->get_QoSCapacity()) / actual_rate;
				}
				if (total_capacity > all_capacity) { // cannot place when QoS cannot achieve
					decline_calorie = -1;
					break;
				}
			}
			else { // cannot place at a saturated vertex
				decline_calorie = -1; // this composition is invalid
				break; // cannot deploy at saturated vertex, pass this combination
			}
		}
		it_c++;
		it_p++;
	}
	return decline_calorie;
}
string deploy_sfc(graph &GlobalMap, vector<vertex> path, SFC s, vector<vertex> &frontSet) {
	int branch_idx = s.get_branchidx(), branch_vertexNum = -1;
	int it_vnf = 0;
	bool is_main = false, frontSet_done = false;
	string deploy_result = "fail";
	vector<VNF> list_vnf = s.get_VNFList();
	vector<vector<int>> set_intcom = integerComposition(path.size(), 2, s.get_VNFList().size());
	vector<pair<double, vector<int>>> set_placement;
	vector<int> best_composition;
	// show_icSet(set_intcom);
	for (vector<vector<int>>::iterator it = set_intcom.begin(); it != set_intcom.end(); it++) { // make pair
		double this_decline = get_declineCalorie(path, s, *it);
		if (this_decline >= 0) set_placement.push_back(make_pair(this_decline, *it)); // if composition valid, push back vector
		else continue;
	}
	if (set_placement.size() > 0) { // can deploy
		sort(set_placement.begin(), set_placement.end(), cmp_placement);
		// show_icSetPair(set_placement);
		best_composition = set_placement.begin()->second;
		cout << "\t";
		show_ic(best_composition);
		vector<int>::iterator it_c = best_composition.begin();
		vector<vertex>::iterator it_p = path.begin();
		if (branch_idx > 0 && !s.is_sub()) is_main = true; // is main branch
		while (it_p != path.end()) { // deployment
			for (int n = 0; n < *it_c; n++) {
				GlobalMap.get_vertexIt(it_p->get_Num())->add_VNF(list_vnf[it_vnf]);
				if (is_main && it_vnf == branch_idx) {
					frontSet_done = true;
					branch_vertexNum = it_p->get_Num(); // record branch vnf position
				}
				else if (!frontSet_done) {
					frontSet.push_back(*(GlobalMap.get_vertexIt(it_p->get_Num())));
				}
				it_vnf++;
			}
			it_c++;
			it_p++;
		}
		if (is_main) deploy_result = toStr(branch_vertexNum); // feedback vertexNum of branching
		else deploy_result = "success";
	}
	else { // mark all vertex of this path as saturated vertex
		for (vector<vertex>::iterator it = path.begin(); it != path.end(); it++) {
			// if (!it->is_pure()) 
			GlobalMap.get_vertexIt(it->get_Num())->update_saturated();
		}
	}
	return deploy_result;
}

// test area ==========
void show_deployment(graph G) {
	int size = G.get_size();
	cout << "Deployment:" << endl;
	for (int n = 0; n < size; n++) {
		double this_offload = 0;
		vector<vertex>::iterator it = G.get_vertexIt(n);
		cout << "#" << it->get_Num() << "\t";
		if (it->is_pure()) cout << "No vnf" << endl;
		else {
			vector<VNF> placement = it->get_VNFPlacement();
			cout.setf(ios::fixed);
			for (vector<VNF>::iterator it_p = placement.begin(); it_p != placement.end(); it_p++) {
				double actual_rate = get_actualRate(placement, it_p, "isG");
				cout << setprecision(2) << "vnf_" << it_p->get_type() << "[" << it_p->get_QoSCapacity() << "(" << actual_rate * 100 << "%)]\t";
				this_offload += (it_p->get_QoSCapacity()) / actual_rate;
			}
			cout << setprecision(2) << endl << " - Offload: " << this_offload << "/" << it->get_capacity() << "(" << 100 * this_offload / (it->get_capacity()) << "%)" << endl;
			cout.unsetf(ios::fixed);
		}
		cout << endl;
	}
}
// test area ==========

int main() {
	graph G;
	vector<SFC> S;
	map<int, VNF> vnfLib;
	map<int, double> vnfSet_capacity;
	string datapath_vnf = "C:/Users/42947/Desktop/Sample-NFV/vnf_data";
	string datapath_sfc = "C:/Users/42947/Desktop/Sample-NFV/sfc_data";
	string datapath_graph = "C:/Users/42947/Desktop/Sample-NFV/graph_data";
	//initializing
	init_getdata(datapath_vnf, vnfLib, datapath_sfc, S, datapath_graph, G);
	inti_getcalorie(S, vnfSet_capacity);
	sort_sfc(S);
	show_input(vnfSet_capacity, vnfLib, S, G);
	int sfc_num = S.size();
	int process = 0;
	int sub_startpoint = -1;
	bool is_this_sub = false;
	string deploy_feedback;
	vector<vertex> path;
	vector<vertex> frontSet;
	while (process < sfc_num) {
		path = find_path(G, S[process], sub_startpoint, frontSet);

		// debug
		cout << endl << "path of #" << S[process].get_num() << ": ";
		show_path(path);

		deploy_feedback = deploy_sfc(G, path, S[process], frontSet);

		cout << "feedback: " << deploy_feedback << endl << endl;

		if (deploy_feedback == "success") { // success, deploy next
			sub_startpoint = -1;
			frontSet.clear();
			process++;
		}
		else if (deploy_feedback == "fail") { // fail, find path again
			continue;
		}
		else { // success and this is main_branch
			sub_startpoint = toNum<int>(deploy_feedback);
			process++;
		}
	}
	cout << "testing done" << endl;
	show_deployment(G);
	system("pause");
	return 0;
}