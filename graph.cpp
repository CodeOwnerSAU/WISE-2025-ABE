#include "graph.h"

int Graph::readBeijingMapDirected(string filename)
{
	ifstream inGraph(filename.c_str());
	if(!inGraph)
	{
		cout << "Cannot open Beijing Map " << filename << endl;
		return -1;
	}

	int nodeNum, i;
	inGraph >> nodeNum >> minX >> maxX >> minY >> maxY;
	cout << "Beijing Node number: " << nodeNum << endl;
	this->nodeNum = nodeNum;

	vector<Edge> vEdges; 
	vCoor.reserve(nodeNum);
	

	double x, y;
	int nodeID, type, j, k;
	int ID2, length;
	for(i = 0; i < nodeNum; i++)
	{
		inGraph >> nodeID >> type >> x >> y >> j >> j; 
		vCoor.push_back(make_pair(x, y));

		for(k = 0; k < j; k++)
		{
			inGraph >> ID2 >> length; 
			struct Edge e, eR;
			e.ID1 = nodeID;
			e.ID2 = ID2;
			e.length = length;
			vEdges.push_back(e);

			eR.ID1 = ID2;
			eR.ID2 = nodeID;
			eR.length = length;
			vEdgeR.push_back(eR);
		}

		inGraph >> j;
		for(k = 0; k < j; k++)
			inGraph >> ID2 >> length;

		inGraph >> j;
		for(k = 0; k < j; k++)
			inGraph >> ID2 >> length;
	}

	cout << "Finish Reading " << filename << endl;
	inGraph.close();

	adjList.resize(nodeNum);
	adjListR.resize(nodeNum);
	adjListEdge.resize(nodeNum);
	adjListEdgeR.resize(nodeNum);
	int edgeCount = 0;
	for(auto &edge: vEdges)
	{
		int ID1 = edge.ID1;
		int ID2 = edge.ID2;
		int length = edge.length;

		bool b = false;
		for(auto ivp = adjList[ID1].begin(); ivp != adjList[ID1].end(); ivp++)
		{
			if((*ivp).first == ID2) 
			{
				b = true;
				break;
			}
		}

		if(!b) 
		{
			adjList[ID1].push_back(make_pair(ID2, length)); 
			adjListR[ID2].push_back(make_pair(ID1, length));

			edge.edgeID = edgeCount; 
			cout << ID2 << "\t" << ID1 << "\t" << edgeCount << endl;
			adjListEdge[ID1].push_back(make_pair(ID2, edgeCount)); 
			adjListEdgeR[ID2].push_back(make_pair(ID1, edgeCount));
			vEdge.push_back(edge);
			edgeCount++;
		}
	}

	cout << "Beijing Road number: " << edgeCount << endl;

	return nodeNum;
}
//readUS
int Graph::readUSMap(const char * edgeFile,const char *nodeFile)
{
    nodeNum=0;
    adjList.clear();
    adjListEdge.clear();
    QueryWord.clear();
    NodesBit.clear();
    RSP.clear();
    KEYS.clear();
    Nodes.clear();
	ifstream inGraph(edgeFile);
	if(!inGraph)
		cout << "Cannot open Map " << edgeFile<< endl;
	cout << "Reading " << edgeFile << endl;

	string line;
    getline(inGraph,line);
    int ID1, ID2, length;
    unsigned long index=line.find(' ');
    nodeNum=stoi(line.substr(0,index));
    adjList.resize(nodeNum);
    adjListR.resize(nodeNum);
   adjListEdge.resize(nodeNum);
    adjListEdgeR.resize(nodeNum);
    int edgeCount = 0;
    char a;
	while(!inGraph.eof())
	{
		inGraph >> a >> ID1 >> ID2 >> length;
        //ID start from 0 because in file ID from start 1 ,thus need -=
		ID1 -= 1;
		ID2 -= 1;
		
		struct Edge e; 
		e.ID1 = ID1;
		e.ID2 = ID2;
		e.length = length;
		e.edgeID = edgeCount; 
		
		bool bExisit = false;
		for(int i = 0; i < (int)adjList[ID1].size(); i++) 
		{
			if(adjList[ID1][i].first == ID2)
			{
				bExisit = true;
				break;
			}
		}

		//cout << ID1 << "\t" << ID2 << "\t" << length << endl;
		if(!bExisit)
		{
			vEdge.push_back(e);
			adjList[ID1].push_back(make_pair(ID2, length));
			adjListR[ID2].push_back(make_pair(ID1, length));
			adjListEdge[ID1].push_back(make_pair(ID2, edgeCount)); 
			adjListEdgeR[ID2].push_back(make_pair(ID1, edgeCount));
			edgeCount++;
		}
		//getline(inGraph,line);
	}

	vbISO.assign(nodeNum, false); 
	inGraph.close();



    //set RSP
    for(int i = 0;i < KIND; i++){
        vector< int >tmp;
        RSP.push_back(tmp);
    }
    //set nodeKEYS

    fstream fp(nodeFile);
    int nid;
    string ss;
    double x,y;
    unordered_map<int,int>ma;
    KEYS.push_back(ma);
    Nodes.resize(nodeNum);
    cout<<"finidshed"<<endl;
    PointList.resize(nodeNum);

    while(fp >> nid >> x >> y >> ss){
        PointList[nid-1].id=nid-1;
        PointList[nid-1].coord[0]=(int)x;
        PointList[nid-1].coord[1]=(int)y;
        Node node={x,y};
        node.keyword.reset();
        //cout<<nid<<" "<<ss<<endl;
        char str[3000]; //according to the keywordstr.length (20*4)
        strcpy(str,ss.c_str());
        const char * split = ",";
        char * p;
        p = strtok (str,split);
        unordered_map<int,int>mp;
        mp.clear();
        bitset<KIND> bit;
        bit.reset();
        int count=0;
        //vector<int> temp;
        while(p!=NULL) {
            if(strcmp(p,"-1") ==0){ break;}
            mp.insert(make_pair(atoi(p),1));
            RSP[atoi(p)].push_back(nid-1);
            bit.set(atoi(p));
            node.keyword.set(atoi(p));
            node.kwList.push_back(atoi(p));
            p = strtok(NULL,split);
            count++;
        }
        if(count!=0){
            POIObjects.push_back(nid-1);
        }
        KEYS.push_back(mp);
        NodesBit.push_back(bit);
        Nodes[nid-1]=node;
    }
//	mCoor["NY"] = make_pair(84000, 69000);
	cout << "--------Graph loading finished-----!!!" << endl;
	return nodeNum;
}

int Graph::readUSMapCost(string filename)
{ 
	ifstream inGraph(filename);
	if(!inGraph)
		cout << "Cannot open Map " << filename << endl; 
	cout << "Reading " << filename << endl;

	string line; 
//	int eNum;
	inGraph >> nodeNum;  
//	inGraph >> nodeNum >> eNum; 

	int ID1, ID2, length, cost;
	adjList.resize(nodeNum);
	adjListR.resize(nodeNum);
	adjListEdge.resize(nodeNum);
	adjListEdgeR.resize(nodeNum); 
	adjListCost.resize(nodeNum);
	adjListCostR.resize(nodeNum);
	int edgeCount = 0;
//	cost = 0;
	while(!inGraph.eof())
	{
		inGraph >> ID1 >> ID2 >> length >> cost;
//		inGraph >> ID1 >> ID2 >> length; 
		ID1 -= 1; 
		ID2 -= 1;  

		struct Edge e; 
		e.ID1 = ID1;
		e.ID2 = ID2;
		e.length = length;
		e.edgeID = edgeCount;  
		e.cost = cost;

		bool bExisit = false;
		for(int i = 0; i < (int)adjListCost[ID1].size(); i++) 
		{
			if(adjListCost[ID1][i].first == ID2)
			{
				bExisit = true;
				if(cost < adjListCost[ID1][i].second)
				{
//					cout << "Cost update:" << adjListCost[ID1][i].second << "\t" << cost << endl; 
					adjListCost[ID1][i].second = cost; 
					for(int j = 0; j < (int)adjListCostR[ID2].size(); j++)
					{
						if(adjListCostR[ID2][j].first == ID1)
						{
							adjListCostR[ID2][j].second = cost;    
							break;
						}
					}

					int eID = adjListEdge[ID1][i].second;
					vEdge[eID].cost = cost;
				}
				break;
			}
		}

		if(!bExisit)
		{
			vEdge.push_back(e);
			adjList[ID1].push_back(make_pair(ID2, length));
			adjListR[ID2].push_back(make_pair(ID1, length));
			adjListEdge[ID1].push_back(make_pair(ID2, edgeCount)); 
			adjListEdgeR[ID2].push_back(make_pair(ID1, edgeCount));
			adjListCost[ID1].push_back(make_pair(ID2, cost));
			adjListCostR[ID2].push_back(make_pair(ID1, cost)); 
			edgeCount++;
		}
	//	getline(inGraph,line);
	}

/*	for(int i = 0; i < (int)adjListCost.size(); i++)
	{
		unordered_map<int, int> us;
		for(auto& c: adjListCost[i])
		{
			if(us.find(c.first) != us.end())
			{
				cout << "Repeated! " << i << "\t" << c.first << "\t" << us[c.first] << "\t" << c.second << endl;
			}
			else
				us[c.first] = c.second;
		}
	}
*/
	vbISO.assign(nodeNum, false); 
	inGraph.close();
/*	
	vector<int> vCorrect;
	ifstream ic("./c"); 
	int icc;
	while(ic >> icc)
		vCorrect.push_back(icc); 
	int c= 0;
	int d = 0;
	int oldv = -1; 
	for(auto&v : vCorrect)
	{
		cout << v << "\t";
		bool b = false;
		if(oldv > -1)
		{
			for(int j= 0 ; j < (int)adjList[oldv].size(); j++)
			{
				if(adjList[oldv][j].first == v)
				{
					d += adjList[oldv][j].second;
					c += adjListCost[oldv][j].second;
					b = true;
					break;
				}
			}
			if(!b)
			{
				cout << endl << "Disconnected!" << endl;
			}
		}
		oldv = v;
		cout << endl << d << "\t" << c << "\t" << endl;
	}
	*/

	return nodeNum;
}

int Graph::readExampleMap(string filename)
{ 
	ifstream inGraph(filename);
	if(!inGraph)
		cout << "Cannot open Map " << filename << endl; 
	cout << "Reading " << filename << endl;

	string line;
	do
	{
		getline(inGraph,line);
		if(line[0]=='p')
		{ 
			vector<string> vs;
			boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
			nodeNum = stoi(vs[2]); //edgenum=stoi(vs[3]);
			cout << "Nodenum " << nodeNum<<endl;
			edgeNum = 0;
		}
	}while(line[0]=='c'|| line[0]=='p');

	int ID1, ID2, length;
	adjList.resize(nodeNum);
	adjListR.resize(nodeNum);
	adjListEdge.resize(nodeNum);
	adjListEdgeR.resize(nodeNum); 
	int edgeCount = 0;
	while(!inGraph.eof())
	{
		vector<string> vs;
		boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
		ID1 = stoi(vs[1]); 
		ID2 = stoi(vs[2]); 
		length = stoi(vs[3]);
		
		struct Edge e; 
		e.ID1 = ID1;
		e.ID2 = ID2;
		e.length = length;
		e.edgeID = edgeCount; 
		vEdge.push_back(e);

		cout << ID1 << "\t" << ID2 << "\t" << length << endl;
		adjList[ID1].push_back(make_pair(ID2, length));
		adjListR[ID2].push_back(make_pair(ID1, length));
		adjListEdge[ID1].push_back(make_pair(ID2, edgeCount)); 
		adjListEdgeR[ID2].push_back(make_pair(ID1, edgeCount));
		edgeCount++;
		getline(inGraph,line);
	}

	vbISO.assign(nodeNum, false); 
	inGraph.close();

//	mCoor["NY"] = make_pair(84000, 69000);

	return nodeNum;
}
	

void Graph::set_nodeKEYS_NodesBit(string filename){
	//set RSP
	for(int i = 0;i < KIND; i++){
		vector< int >tmp;
		RSP.push_back(tmp);
	}

	//set nodeKEYS
	
	fstream fp(filename);
        int nid;
        string ss;
        long x,y;
        unordered_map<int,int>ma;
        KEYS.push_back(ma);
        while(fp >> nid >> x >> y >> ss){
        	    //cout<<nid<<" "<<ss<<endl;
                char str[100]; //according to the keywordstr.length (20*4)
                strcpy(str,ss.c_str());
                const char * split = ",";
                char * p;
                p = strtok (str,split);
                unordered_map<int,int>mp;
                mp.clear();
                bitset<KIND> bit;
                bit.reset();
                while(p!=NULL) {
                                if(strcmp(p,"-1") ==0)
                                        break;
                                mp.insert(make_pair(atoi(p),1));
                                RSP[atoi(p)].push_back(nid-1);
                                bit.set(atoi(p));
                                p = strtok(NULL,split);
                }
                KEYS.push_back(mp);
                NodesBit.push_back(bit);
        }
		cout<<"--------KEYS map and Query init completed-----!!!"<<endl;
        fp.close();
}
bitset<KIND> Graph::getempPathBit(vector<int> temp) {

    bitset<KIND> result;
    for(auto i:temp){
        result|=NodesBit[i];
    }
    return result&QueryBit;
}
int Graph::Clen(int S, int T, bitset<KIND>& query, vector<int> &result){
	// cout<<"-------Clen----------"<<endl;
	set< int > minCandtmp;
	vector< int >minCand;
	for(int i = 0;i < QueryWord.size();i++){
		int MIN_LEN = INT_MAX;
	        int PRE = 40;//the number of node nearest S and T for one keyword
		int min_pos;
		// if(QueryBit.test(QueryWord[i])){
		// 	cout<<QueryWord[i]<<endl;
			int pos = QueryWord[i];
			vector<pair<int,double>>pri_node;
			pri_node.clear();
			for(int k = 0;k < RSP[pos].size();k++){
				int _len = distanceQuery(S+1, RSP[pos][k]+1) + distanceQuery(RSP[pos][k]+1, T+1);
				pri_node.push_back(make_pair(RSP[pos][k], _len));
			}
			sort(pri_node.begin(),pri_node.end(),[](const pair<int, double>& p,const pair<int, double>& q){return p.second <= q.second;});
			if(PRE > RSP[pos].size()){PRE = RSP[pos].size();}
			for(int j = 0;j < PRE ;j++){
				// int len = Dijkstra(S, pri_node[j].first) + Dijkstra(pri_node[j].first, T);
				int len = AStar(S, pri_node[j].first) + AStar(pri_node[j].first, T);
				// int len = distanceQuery(S+1, pri_node[j].first+1) + distanceQuery(pri_node[j].first+1, T+1);
				if(MIN_LEN > len){
					MIN_LEN = len;
					min_pos = pri_node[j].first;
				}
			}	
			minCandtmp.insert(min_pos);
	}
//	cout<<minCandtmp.size()<<endl;
	for(set<int>::iterator it = minCandtmp.begin(); it!=minCandtmp.end(); it++){
//		cout<<*it<<endl;
		minCand.push_back(*it);
	}
//	cout<<endl;
//	cout<<"cand.size() = "<<minCand.size()<<endl;
	for(vector<int>::iterator it = minCand.begin(); it!=minCand.end(); it++){
//		cout<<*it<<endl;
	}
//	cout<<endl;
	int dist = 0;
	int index;
	int Start = S;
	result.push_back(S);
//	cout<<Start<<" -> ";
	while(minCand.size() > 0){
		int min = INT_MAX;
		for(int i=0;i<minCand.size();i++){
			// cout<<"mincand["<<i<<"]="<<minCand[i]<<endl;
			// int len = Dijkstra(Start, minCand[i]);
			int len = AStar(Start, minCand[i]);
			// int len = distanceQuery(Start + 1, minCand[i] + 1);
			if(min > len){
				min = len;
				index = minCand[i];
			}
		}
		result.push_back(index);
		dist += min;
//		cout<<min<<endl;
//		cout<<index<<"("<<minCand.size()<<") -> ";
		for(vector<int>::iterator it = minCand.begin(); it!=minCand.end();){
			if(*it == index){
				it = minCand.erase(it);
				break;
			}else{it++;}
		}
		Start = index; 
	}
//	cout<<T<<endl;;
	dist += AStar(Start, T);
	// dist += Dijkstra(Start, T);
	// dist += distanceQuery(Start + 1, T + 1);
	result.push_back(T);
	return dist;
}

int Graph::readUSCost(string filename)
{
	ifstream inGraph(filename);
	if(!inGraph)
		cout << "Cannot open Map " << filename << endl; 
	cout << "Reading " << filename << endl;

	string line;
	do
	{
		getline(inGraph,line);
	}while(line[0]=='c'||line[0]=='p');

	int ID1, ID2, cost; 
	adjListCost.resize(nodeNum);
	adjListCostR.resize(nodeNum);
	int edgeCount = 0;
	while(!inGraph.eof())
	{
		vector<string> vs;
		boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
		ID1 = stoi(vs[1]) - 1;  
		ID2 = stoi(vs[2]) - 1; 
		cost = stoi(vs[3]);
		
		adjListCost[ID1].push_back(make_pair(ID2, cost));
		adjListCostR[ID2].push_back(make_pair(ID1, cost)); 
		vEdge[edgeCount].cost = cost;
		edgeCount++;
		getline(inGraph,line);
	}
	
	cout << vEdge.size() << endl;
	return nodeNum;
}

int Graph::DijkstraPath(int ID1, int ID2, vector<int>& vPath, vector<int>& vPathEdge)
{
	benchmark::heap<2, int, int> queue(adjList.size());
	queue.update(ID1, 0);

	vector<int> vDistance(adjList.size(), INF);
	vector<int> vPrevious(adjList.size(), -1);
	vector<int> vPreviousEdge(adjList.size(), -1);
	vector<bool> vbVisited(adjList.size(), false);
	int topNodeID, neighborNodeID, neighborLength, neighborRoadID;

	vDistance[ID1] = 0;

	while(!queue.empty())
	{
		int topDistance; 
		queue.extract_min(topNodeID, topDistance); 
		vbVisited[topNodeID] = true;
		if(topNodeID == ID2)
			break;

		for(int i = 0; i < (int)adjList[topNodeID].size(); i++)
		{
			neighborNodeID = adjList[topNodeID][i].first;
			neighborLength = adjList[topNodeID][i].second; 
			neighborRoadID = adjListEdge[topNodeID][i].second;
			int d = vDistance[topNodeID] + neighborLength;
			if(!vbVisited[neighborNodeID])
			{
				if(vDistance[neighborNodeID] > d)
				{
					vDistance[neighborNodeID] = d;
					queue.update(neighborNodeID, d);
					vPrevious[neighborNodeID] = topNodeID;
					vPreviousEdge[neighborNodeID] = neighborRoadID;
				}
			}
		}
	}

	vPath.clear();
	vPathEdge.clear();
	vPath.push_back(ID2);
	int p = vPrevious[ID2];
	int e = vPreviousEdge[ID2];
	while(p != -1)
	{
		vPath.push_back(p);
		vPathEdge.push_back(e);
		e = vPreviousEdge[p];
		p = vPrevious[p];
	}

//	if(vPathEdge.size() > 0)
//		vPathEdge.erase(vPathEdge.end()-1);

	reverse(vPath.begin(), vPath.end());
	reverse(vPathEdge.begin(), vPathEdge.end());

	return vDistance[ID2];
}

int Graph::DijkstraPath2(int ID1, int ID2, unordered_set<int>& sRemovedNode, vector<int>& vPath, vector<int>& vPathEdge)
{
	benchmark::heap<2, int, int> queue(adjList.size());
	queue.update(ID1, 0);

	vector<int> vDistance(adjList.size(), INF);
	vector<int> vPrevious(adjList.size(), -1);
	vector<int> vPreviousEdge(adjList.size(), -1);
	vector<bool> vbVisited(adjList.size(), false);
	int topNodeID, neighborNodeID, neighborLength, neighborRoadID;

	vDistance[ID1] = 0;

	while(!queue.empty())
	{
		int topDistance; 
		queue.extract_min(topNodeID, topDistance); 
		vbVisited[topNodeID] = true;
		if(topNodeID == ID2)
			break;

		for(int i = 0; i < (int)adjList[topNodeID].size(); i++)
		{
			neighborNodeID = adjList[topNodeID][i].first; 
			if(sRemovedNode.find(neighborNodeID) != sRemovedNode.end())
				continue;
			neighborLength = adjList[topNodeID][i].second; 
			neighborRoadID = adjListEdge[topNodeID][i].second;
			int d = vDistance[topNodeID] + neighborLength;
			if(!vbVisited[neighborNodeID])
			{
				if(vDistance[neighborNodeID] > d)
				{
					vDistance[neighborNodeID] = d;
					queue.update(neighborNodeID, d);
					vPrevious[neighborNodeID] = topNodeID;
					vPreviousEdge[neighborNodeID] = neighborRoadID;
				}
			}
		}
	}

	vPath.clear();
	vPathEdge.clear();
	vPath.push_back(ID2);
	int p = vPrevious[ID2];
	int e = vPreviousEdge[ID2];
	while(p != -1)
	{
		vPath.push_back(p);
		vPathEdge.push_back(e);
		e = vPreviousEdge[p];
		p = vPrevious[p];
	}

//	if(vPathEdge.size() > 0)
//		vPathEdge.erase(vPathEdge.end()-1);

	reverse(vPath.begin(), vPath.end());
	reverse(vPathEdge.begin(), vPathEdge.end());

	return vDistance[ID2];
}

int Graph::Dijkstra(int ID1, int ID2)
{
	benchmark::heap<2, int, int> queue(adjList.size());
	queue.update(ID1, 0);

	vector<int> vDistance(nodeNum, INF);
	vector<bool> vbVisited(nodeNum, false);
	int topNodeID, neighborNodeID, neighborLength;
	vector<pair<int, int> >::iterator ivp;

	vDistance[ID1] = 0;

	compareNode cnTop;
	while(!queue.empty())
	{
		int topDistance; 
		queue.extract_min(topNodeID, topDistance); 
		vbVisited[topNodeID] = true; 
//		cout << topNodeID << "\t" << vDistance[topNodeID] << endl;
		if(topNodeID == ID2)
			break;

		for(ivp = adjList[topNodeID].begin(); ivp != adjList[topNodeID].end(); ivp++)
		{
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second; 
			int d = vDistance[topNodeID] + neighborLength;
			if(!vbVisited[neighborNodeID])
			{
				if(vDistance[neighborNodeID] > d)
				{
					vDistance[neighborNodeID] = d;
					queue.update(neighborNodeID, d);
				}
			}
		}
	}

	return vDistance[ID2];
}

int Graph::AStar(int ID1, int ID2)
{
	benchmark::heap<2, int, int> queue(adjList.size());
	vector<int> vDistance(adjList.size(), INF);
	vector<bool> vbVisited(adjList.size(), false);
	int topNodeID, neighborNodeID, neighborLength;
	vector<pair<int, int> >::iterator ivp;

	queue.update(ID1, EuclideanDistance(ID1, ID2)); 
	vDistance[ID1] = 0;

	compareNode cnTop;
	while(!queue.empty())
	{
		int topDistance;
		queue.extract_min(topNodeID, topDistance);
		vbVisited[topNodeID] = true;

		if(topNodeID == ID2)
			break;

		for(ivp = adjList[topNodeID].begin(); ivp != adjList[topNodeID].end(); ivp++)

		{
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second; 
			int d = vDistance[topNodeID] + neighborLength;
			if(!vbVisited[neighborNodeID])
			{
				if(vDistance[neighborNodeID] > d)
				{
					vDistance[neighborNodeID] = d;
					queue.update(neighborNodeID, d+EuclideanDistance(neighborNodeID, ID2));
				}
			}
		}
	}

	return vDistance[ID2];
}

int Graph::AStarPath(int ID1, int ID2, vector<int>& vPath, vector<int>& vPathEdge, string& city)
{
	benchmark::heap<2, int, int> queue(adjList.size());
	vector<int> vDistance(adjList.size(), INF);
	vector<int> vPrevious(adjList.size(), -1);
	vector<int> vPreviousEdge(adjList.size(), -1);
	vector<bool> vbVisited(adjList.size(), false);
	int topNodeID, neighborNodeID, neighborLength, neighborRoadID;

	int latU, lonU;
	if(city == "US")
	{
		lonU = 84000;
		latU = 69000;
	}
	else
	{
		lonU = 83907;
		latU = 111319; 
	}

	queue.update(ID1, EuclideanDistance(ID1, ID2)); 
	vDistance[ID1] = 0;

	compareNode cnTop;
	while(!queue.empty())
	{
		int topDistance;
		queue.extract_min(topNodeID, topDistance);
		vbVisited[topNodeID] = true;

		if(topNodeID == ID2)
			break;

		for(int i = 0; i < (int)adjList[topNodeID].size(); i++)
		{
			neighborNodeID = adjList[topNodeID][i].first;
			neighborLength = adjList[topNodeID][i].second; 
			neighborRoadID = adjListEdge[topNodeID][i].second;
			int d = vDistance[topNodeID] + neighborLength;
			if(!vbVisited[neighborNodeID])
			{
				if(vDistance[neighborNodeID] > d)
				{
					vDistance[neighborNodeID] = d;
					queue.update(neighborNodeID, d+EuclideanDistance(neighborNodeID, ID2));
					vPrevious[neighborNodeID] = topNodeID;
					vPreviousEdge[neighborNodeID] = neighborRoadID;
				}
			}
		}
	}
	
	vPath.clear();
	vPathEdge.clear();
	vPath.push_back(ID2);
	int p = vPrevious[ID2];
	int e = vPreviousEdge[ID2];
	while(p != -1)
	{
		vPath.push_back(p);
		vPathEdge.push_back(e);
		e = vPreviousEdge[p];
		p = vPrevious[p];
	}

	reverse(vPath.begin(), vPath.end());
	reverse(vPathEdge.begin(), vPathEdge.end());

	return vDistance[ID2];
}
double degreesToRadians(double degrees) {
    return degrees * M_PI / 180.0;
}
double Graph::EuclideanDistance(int ID1, int ID2)
{
    double div = 1000000.0;
    double lon1 = ((double) Nodes[ID1].x) / div;
    double lat1 = ((double) Nodes[ID1].y) / div;
    double lon2 = ((double) Nodes[ID2].x) / div;
    double lat2 = ((double) Nodes[ID2].y) / div;
    lat1 = degreesToRadians(lat1);
    lon1 = degreesToRadians(lon1);
    lat2 = degreesToRadians(lat2);
    lon2 = degreesToRadians(lon2);
    double lat_diff = lat2 - lat1;
    double lon_diff = lon2 - lon1;

    double a = std::sin(lat_diff / 2) * std::sin(lat_diff / 2) +
               std::cos(lat1) * std::cos(lat2) *
               std::sin(lon_diff / 2) * std::sin(lon_diff / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
    return 63710000.0* c;
}

inline int Graph::EuclideanDistanceAdaptive(int ID1, int ID2, int latU, int lonU)
{
	int lat=(int)(abs(vCoor[ID1].first - vCoor[ID2].first)*latU);
	int lon=(int)(abs(vCoor[ID1].second - vCoor[ID2].second)*lonU);
	int min, max;
	min = (lat > lon) ? lon : lat;
	max = (lat > lon) ? lat : lon;
	int approx = max*1007 + min*441;
	if(max < (min << 4))
		approx -= max * 40;
	return (approx + 512) >> 10;
}

int Graph::ISONodes()
{
	srand (time(NULL));
	vbISOF.assign((int)adjList.size(), false);	//forward
	vbISOB.assign((int)adjList.size(), false);	//backward
	vbISOU.assign( (int)adjList.size(), false);	//F & B
	vbISO.assign((int)adjList.size(), false);	//F | B

	ifstream ifISOF("./beijingISOF");
	ifstream ifISOB("./beijingISOB");
	ifstream ifISOU("./beijingISOU");
	ifstream ifISO("./beijingISO");

	if(!ifISOF || !ifISOB || !ifISOU || !ifISO)
	{
		//ISO files do not exist
		//Create new ones
		srand(time(NULL));
		cout << "Identifying ISO Nodes" << endl;
		for(int i = 0; i < 10; i++)
		{
			int nodeID = rand() % adjList.size();
			cout <<nodeID <<endl;
			vector<bool> vbVisitedF;
			vector<bool> vbVisitedB;
			int vnF = BFS(nodeID, true, vbVisitedF);
			int vnB = BFS(nodeID, false, vbVisitedB);
			cout <<vnF <<"\t" << vnB <<endl;

			if(vnF < 1000 || vnB < 1000)
				continue;

			for(int j = 0; j < (int)adjList.size(); j++)
			{
				if(!vbVisitedF[j])
					vbISOF[j] = true;

				if(!vbVisitedB[j])
					vbISOB[j] = true;

				if(!vbVisitedF[j] || !vbVisitedB[j])
					vbISOU[j] = true;

				if(!vbVisitedF[j] && !vbVisitedB[j])
					vbISO[j] = true;
			}
		}

		ofstream ofF("./beijingISOF");
		ofstream ofB("./beijingISOB");
		ofstream ofU("./beijingISOU");
		ofstream of("./beijingISO");

		for(int i = 0; i < (int)adjList.size(); i++)
		{
			if(vbISOF[i])
				ofF << i << endl;

			if(vbISOB[i])
				ofB << i << endl;

			if(vbISOU[i])
				ofU << i << endl;

			if(vbISO[i])
				of << i << endl;
		}

		ofF.close();
		ofB.close();
		ofU.close();
		of.close();

		return 0;
	}
	else
	{
		int nodeID;
		cout << "Loading ISO Nodes" << endl;
		while(ifISOF >> nodeID)
			vbISOF[nodeID] = true;

		while(ifISOB >> nodeID)
			vbISOB[nodeID] = true;

		while(ifISOU >> nodeID)
			vbISOU[nodeID] = true;

		while(ifISO >> nodeID)
			vbISO[nodeID] = true;

		return 0;
	}
}

int Graph::BFS(int nodeID, bool bF, vector<bool>& vbVisited)
{
	vbVisited.assign(adjList.size()+1, false);
	queue<int> Q;
	Q.push(nodeID);
	vbVisited[nodeID] = true;
	int count = 0;
	while(!Q.empty())
	{
		int topNodeID = Q.front();
		Q.pop();
		count++;
		if(bF)
		{
			for(auto& it : adjList[topNodeID])
				if(!vbVisited[it.first])
				{
					vbVisited[it.first] = true;
					Q.push(it.first);
				}
		}
		else
		{
			for(auto& it : adjListR[topNodeID])
				if(!vbVisited[it.first])
				{
					vbVisited[it.first] = true;
					Q.push(it.first);
				}
		}
	}

	return count;
}
	
int Graph::readUSCoor(string filename) 
{
	vCoor.reserve(nodeNum);
	ifstream inFile(filename.c_str()); 
	if(!inFile)
	{
		cout << "Cannot open " << filename << endl;
		return -1;
	}

	string s;
	int nodeID;
	double lon, lat;
	while(inFile >> s)
	{
		if(s == "v") 
		{
			inFile >> nodeID >> lon >> lat;
			lon = -lon / 1000000;
			lat = lat / 1000000;
			vCoor.push_back(make_pair(lon, lat));
		}
	}

	return 0;
}
void Graph::SPT(int root, vector<int>& vSPTDistance, vector<int>& vSPTParent, vector<int>& vSPTParentEdge, vector<vector<int> >& vSPT,vector<bitset<KIND> > &vSPTBitS)
{
    benchmark::heap<2, int, int> Queue(nodeNum);
    Queue.update(root, 0);

    vector<bool> vbVisited(nodeNum, false);
    int topNodeID, neighborNodeID, neighborLength, neighborRoadID;
    vector<pair<int, int> >::iterator ivp;

    vSPTDistance[root] = 0;
    compareNode cnTop;
    while(!Queue.empty())
    {
        int topDistance;
        Queue.extract_min(topNodeID, topDistance);
        vbVisited[topNodeID] = true;
        for(int i = 0; i < (int)adjList[topNodeID].size(); i++)
        {
            neighborNodeID = adjList[topNodeID][i].first;
            neighborLength = adjList[topNodeID][i].second;
            neighborRoadID = adjListEdge[topNodeID][i].second;
            int d = vSPTDistance[topNodeID] + neighborLength;
            if(!vbVisited[neighborNodeID])
            {
                if(vSPTDistance[neighborNodeID] > d)
                {
                    vSPTDistance[neighborNodeID] = d;
                    Queue.update(neighborNodeID, d);
                    vSPTParent[neighborNodeID] = topNodeID;
                    vSPTParentEdge[neighborNodeID] = neighborRoadID;
                }
            }
        }
    }
    // cout<<vSPTDistance[1]<<endl;
    // cout<<vSPTDistance[2]<<endl;
    //Construct SPT //father nodeid stores some children nodeid. tree root is ID2

    for(int i = 0; i < nodeNum; i++){
        if(vSPTParent[i] != -1){
            //vSPT[vSPTParent[i]][vSPT[vSPTParent[i]].size()-1]=i;
            vSPT[vSPTParent[i]].push_back(i);
        }
    }
    queue<int>SPTree;
    SPTree.push(root);
    while(!SPTree.empty()){
        int node = SPTree.front();
        if(vSPTParent[node] != -1){
            vSPTBitS[node] |= vSPTBitS[vSPTParent[node]];
        }
        SPTree.pop();
        if(vSPT[node].size()>0){
            for(int i=0;i<vSPT[node].size();i++){
                SPTree.push(vSPT[node][i]);
            }
        }
    }
}
int Graph::buildGreedyPath(int startNode,int endNode,vector<int> POICandidate,vector<int> &Gpath, int &LengthBound){
    Gpath.clear();
    int pathDis=0;
    bitset<KIND> uncover=QueryBit;
    int currentNode=startNode;
    Gpath.push_back(startNode);
    LengthBound=0;
    if(!POICandidate.empty()){
        while (uncover.count()>0){
            int NNdis=INT_MAX;
            int NNPOI;
            for(auto i=POICandidate.begin();i!=POICandidate.end();i++){
                if(std::find(Gpath.begin(), Gpath.end(),*i)==Gpath.end()&&
                   (uncover&NodesBit[*i]).count()>0){
                    int dis= distanceQuery(currentNode+1,*i+1);
                    if(dis<NNdis){
                        NNdis=dis;
                        NNPOI=*i;
                    }
                }
            }
            POICandidate.erase(std::find(POICandidate.begin(), POICandidate.end(),NNPOI));
            int H2HDis= distanceQuery(currentNode+1,NNPOI+1);
            if(LengthBound<H2HDis&&currentNode!=startNode){
                LengthBound=H2HDis;
            }
            pathDis+= H2HDis;
            Gpath.push_back(NNPOI);
            uncover^=(QueryBit&(NodesBit[NNPOI]&uncover));
            currentNode=NNPOI;
        }
    }
    pathDis+= distanceQuery(currentNode+1,endNode+1);
    Gpath.push_back(endNode);
    return pathDis;
}
int Graph::generateMultPath(int s, int t, int pathNum, vector<pair<int, vector<int>>> &allPath, vector<int> greedyPath,
                            map<int, vector<int>> keyNodeMap, vector<vector<int>> &vPaths, vector<int> &vPathsDis) {
    int count=0;
    if(((pathNum)%(greedyPath.size() - 3))!=0){
        count=1;
    }
    int eachNodeKum=((pathNum)/(greedyPath.size() - 3));
    if(eachNodeKum==0){
        eachNodeKum=1;
    }
    bitset<KIND> currentBit;
    for(int i=1;i<=greedyPath.size()-3;i++) {
        map<int,vector<int>> keyMap=keyNodeMap;
        int forNum=eachNodeKum;
        if(count==1&&i==greedyPath.size()-3){
            forNum=pathNum-(greedyPath.size()-4)*eachNodeKum;
        }
        for (int j = 0; j <forNum; j++) {
            int feasible=1;
            currentBit |= NodesBit[greedyPath[i]] & QueryBit;
            bitset<KIND> uncover = (QueryBit ^ currentBit);
            vector<int> currentPath(greedyPath.begin(), greedyPath.begin() + i + 1);
            vector<int> newPOI(greedyPath.begin() + 1, greedyPath.begin() + i + 1);
            int num = uncover.count();
            int flag = 0;

            while (uncover.count() != 0) {
                if(num==0){
                    feasible=0;
                    break;
                }
                vector<int> POIs;
                vector<int> bestPOIList;
                int bestPOI;
                int POIDis = INT_MAX;
                for (auto node: keyMap[num]) {
                    if ((NodesBit[node] & uncover).count() != 0) {
                        int dis = distanceQuery(s + 1, node + 1) + distanceQuery(node + 1, t + 1);
                        if (dis < POIDis) {
                            flag = 1;
                            POIDis = dis;
                            bestPOI = node;
                        }
                    }
                }
                if (flag == 1) {
                    keyMap[num].erase(std::find(keyMap[num].begin(), keyMap[num].end(), bestPOI));
                }
                if (flag == 0) {
                    num--;
                } else {
                    newPOI.push_back(bestPOI);
                    uncover ^= NodesBit[bestPOI] & QueryBit & uncover;
                    num = uncover.count();
                    flag=0;
                }
            }
            if(feasible==0){
                break;
            }
            vector<int> path;
            int A;
            int pathDis = buildGreedyPath(s, t, newPOI, path,A);
            vPaths.push_back(path);
            vPathsDis.push_back(PruneRepeatedPoiPath(path));
            allPath.emplace_back(pathDis,path);
        }
    }
}
void Graph::getKPathByIGTree(int s, int t, int K, vector<int> POICandidate, vector<vector<int>> &vkPath,
                             vector<int> &vkPathDis) {
    vector<vector<int>> kPath;
    vector<int> kPathDis;
    vector<int> PrunedPOIList;
    vector<int> GPath;
    vector<int> GpathS;
    vector<int> GpathT;
    int A,B;
    //first Greedy path
    buildGreedyPath(s,t,POICandidate,GpathS,A);
    int disA=PruneRepeatedPoiPath(GpathS);
    //second Greey path
    buildGreedyPath(t,s,POICandidate,GpathT,B);
    int disB=PruneRepeatedPoiPath(GpathT);
    int pathDis=0;
    if(disA<=disB){
        GPath=GpathS;
        pathDis=disA;
    } else{
        GPath=GpathT;
        pathDis=disB;
    }
    if(K==1){
        vkPath.push_back(GPath);
        vkPathDis.push_back(pathDis);
        //return Top-1 path
        return;
    }
    //build MAP
    map<int,vector<int>> keyNodeMap;
    for(int i=0;i<QueryBit.count();i++){
        vector<int> l;
        keyNodeMap.insert({i+1,l});
    }
    for(auto node:POICandidate){
        int count=(NodesBit[node]&QueryBit).count();
        keyNodeMap[count].push_back(node);
    }
    //save all path first is pathDis ,second is path
    vector<pair<int,vector<int>>> allPath;
    if(K<GPath.size()-3){
        //only one
        generateMultPath(s,t,K, allPath,GPath,keyNodeMap,kPath,kPathDis);
    }else{
        generateMultPath(s,t,K,allPath,GPath,keyNodeMap,kPath,kPathDis);
        kPath.clear();
        kPathDis.clear();
        std::reverse(GPath.begin(), GPath.end());
        generateMultPath(t,s,K,allPath,GPath,keyNodeMap,kPath,kPathDis);
    }
    allPath.emplace_back(disA,GpathS);
    allPath.emplace_back(disB,GpathT);
    //sort select top-K path
    sort(allPath.begin(),allPath.end());
    int count=0;
    for(int i=0;i<allPath.size();i++){
        if(vkPath.size()==K)
            break;
        if(std::find(vkPathDis.begin(), vkPathDis.end(),allPath[i].first)==vkPathDis.end()){
            vkPath.push_back(allPath[i].second);
            vkPathDis.push_back(allPath[i].first);
        }
    }
}

int Graph::PruneRepeatedPoiPath(vector<int> &bestpath) {
    reverse(bestpath.begin(),bestpath.end());
    vector<int> bestpoi;
    bestpoi.push_back(bestpath.front());
    bitset<KIND> qu(Qu);
    for(auto &p : bestpath){
        bitset<KIND>testpoi(NodesBit[p]);
        testpoi &= qu;
        qu ^= testpoi;
        if(testpoi.count() > 0&& std::find(bestpoi.begin(), bestpoi.end(),p)==bestpoi.end()){
            bestpoi.push_back(p);
        };//keyword unique
    }
    bestpoi.push_back(bestpath.back());
    reverse(bestpath.begin(),bestpath.end());
    int da = 0;
    for(int i = 0;i < bestpoi.size() - 1; i++){
        da += distanceQuery(bestpoi[i] + 1, bestpoi[i+1] + 1);
    }
    return da;
}
void Graph::InitQueryInfo(int s,int t,vector<int> Query){
    QueryBit.reset();
    Qu.reset();
    QueryWord.clear();
    for(auto keyword:Query){
        QueryBit.set(keyword);
        Qu.set(keyword);
        if(!NodesBit[s].test(keyword)&&!NodesBit[t].test(keyword)){
            QueryWord.push_back(keyword);
        }
    }
    QueryBit^=(QueryBit&(NodesBit[s]|NodesBit[t]));
}
bool Graph::testDominated(vector<int> &nodePathTable, vector<vector<int>> pathList, vector<bitset<KIND>> pathListBit,
                          vector<int> pathListLength, vector<int> newPath, bitset<KIND> newPathBit, int newPathLength) {
    if(nodePathTable.empty()){
        return false;
    }
    for(auto pathID:nodePathTable){
        if(pathList[pathID][0]==newPath[0]&&pathList[pathID].back()==newPath.back()&&
        pathListLength[pathID]<newPathLength&&
           (pathListBit[pathID]&QueryBit).count()>=((newPathBit&QueryBit).count())){
            return true;
        }
    }
    return false;
}



