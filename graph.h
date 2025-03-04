#ifndef GRAPH_H
#define	GRAPH_H
#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <chrono>
#include "tools.h"
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include "heap.h"
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set> 
#include <algorithm> 
#include <stack>
#include <bitset>
#include <sys/time.h>
#include <xmmintrin.h>
#include <cmath>
#include <queue>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/bind/bind.hpp>
#include <metis.h>
#include "Point.h"
#include "MBR.h"
#include <thread>
using namespace std;
using namespace benchmark; 
#define KIND 500
#define para 1.0
//grid info 8000000

#define  SIZE1 8000000
#define  SIZE2 35000
#define  SIZE3 500000

static const int CELLSIZE=SIZE1;
static const int MINPOINT=4;
//static const int KIND = 4000;	//the number of keywords kind
typedef struct COMPARENODE
{
	pair<int, int> pif;//ID, minCost
	bool operator() (const struct COMPARENODE& a, const struct COMPARENODE& b) const
	{
		return a.pif.second > b.pif.second;
	}
}compareNode;

struct Edge
{
	int ID1, ID2, length, edgeID, cost;
};
typedef struct{
    double x,y;
    vector<int> adjnodes;
    vector<int> adjweight;
    bitset<KIND> keyword;
    vector<int> kwList;
    //bitset<20> keyword;
}Node;
struct PTNode
{
	int nodeID;
	unordered_map<int, int>	PTRoad;	//childrenID, roadID
	unordered_map<int, int> PTChildren; //childrenID, node Pos
};
class Graph
{
public:
	Graph(){}
	Graph(const char * edgeFile,const char * nodeFile)
	{
			readUSMap(edgeFile,nodeFile);
	}

	//True: Only adjList and adjListMap are used
	//False: Also use R versions
	int nodeNum;
	int edgeNum;
    int avgDisSum;
	vector<vector<pair<int, int> > > adjList;		//neighborID, Distance
	vector<vector<pair<int, int> > > adjListR;
	vector<vector<pair<int, int> > > adjListEdge;	//neighbor,edgeID
	vector<vector<pair<int, int> > > adjListEdgeR;

	vector<vector<pair<int, int> > > adjListCost;		//neighborID, cost
	vector<vector<pair<int, int> > > adjListCostR;


	vector<Edge> vEdge;
	vector<Edge> vEdgeR; 

	int readCost(string filename);

    vector<Node> Nodes;
	//Identify the ISO nodes
	vector<bool> vbISOF;	//forward
	vector<bool> vbISOB;	//backward
	vector<bool> vbISOU;	//F & B
	vector<bool> vbISO;		//F | B
	int ISONodes();
	int BFS(int nodeID, bool bF, vector<bool>& vbVisited);

	vector<pair<double, double> > vCoor; //存储经纬度
	unordered_map<string, pair<int, int> > mCoor;
	double minX, minY, maxX, maxY;

	//lon, lat
	int readBeijingMapDirected(string filename); 
	int readUSMap(const char * edgeFile,const char *nodeFile);
    void InitQueryInfo(int s,int t,vector<int> Query);
	int readUSCost(string filename);
	int readUSMapCost(string filename);
	int readExampleMap(string filename);  
	int readUSCoor(string filename);

	//test.cpp
	void testCSP(string filename);
	void SCPT(int root, vector<int>& vSPTDistance, vector<int>& vSPTCost, vector<int>& vSPTParent, vector<int>& vSPTParentEdge, vector<vector<int> >& vSPT, int C);
	void rCPT(int root, int ID1, int C, vector<int>& vrSPTCost, vector<int>& vrSPTDistance);

	// void contractNode(int threshold);
	int Dijkstra(int ID1, int ID2);
	int DijkstraPath(int ID1, int ID2, vector<int>& vPath, vector<int>& vPathEdge);
	int DijkstraPath2(int ID1, int ID2, unordered_set<int>& sRemovedNode, vector<int>& vPath, vector<int>& vPathEdge);
	int AStar(int ID1, int ID2);
	int AStarPath(int ID1, int ID2, vector<int>& vPath, vector<int>& vPathEdge, string& city); 
	double EuclideanDistance(int ID1, int ID2);
	int EuclideanDistanceAdaptive(int ID1, int ID2, int latU, int lonU);

	//Cache (not use)
	unordered_map<string, unordered_map<int, int> > Cache;
	unordered_map<string, unordered_map<int, int> > Cache_node;
	//set KEYS
	bitset<KIND>Qu;		// all query keywords bit
	bitset<KIND>QueryBit; 	// keywords bit not in ID1 and ID2
	vector<int>QueryWord;	// keywords not in ID1 and ID2
	vector<vector<int>> RSP;	// each keyword contain pois (0 ~ KIND-1) RSP[i] id is from 1
	vector<unordered_map<int, int>>KEYS;	//store each poi's keywords using map, and node start from zero
	vector<bitset<KIND>>NodesBit;	// store each node's keywords bit
	void set_nodeKEYS_NodesBit(string filename);	// init Query and each node keywors using map
	int Clen(int S,int T,bitset<KIND> &query, vector<int> &result);	// constraint len by Optimal shortest source-keyword-destination


	// H2H
    const int infinity = INT_MAX;
    const int SIZEOFINT = 4;
    int *toRMQ, *height, **RMQIndex;
    int *belong;
    int root, TreeSize;
    int **rootToRoot, *rootSite;
    int **dis,**dis2, **pos, **pos2;
    int *posSize, *pos2Size;
    int *chSize;
    int ** ch;
    int *LOG2, *LOGD;
    int rootSize;
    int *DFSList, *toDFS;
    int ***BS;
    int *Degree;
    int **Neighbor, **Weight;
    FILE *fin_Index;
    int *EulerSeq;
    long long queryCnt;
    int LCAQuery(int _p, int _q);
    /**
     * when use distanceQuery  p q
     * because in H2H.index is from id 1
     * but in road network id is form 0 ,so id must be id+1
     * @param p
     * @param q
     * @return
     */
    int distanceQuery(int p, int q);
    void readGraph(const char * filename);
    int shortestPathQuery(int p, int q);
    void scanIntArray(int *a, int n);
    int* scanIntVector(int *a);
    void readIndex(const char *file);
    void Init_H2H_Index(const char *index, const char* graph);
    int H2HPathPlusBit(int ID1, int ID2, vector<int>& vPath,vector<bitset<KIND>> &H2HPathBit,int& pBegin,int &pEnd,int pathCount);
    int H2HPath(int ID1, int ID2, vector<int>& vPath,vector<bitset<KIND>> &H2HPathBit);
	// DAmrrp
    int H2HPath(int ID1, int ID2, vector<int>& vPath,bitset<KIND> &H2HPathBit);
	int eKSPNew(int ID1, int ID2, int k, vector<int>& query,vector<int>& kResults, vector<vector<int> >& vkPath);
    //void SPT();
	void SPT(int root, vector<int>& vSPTDistance, vector<int>& vSPTParent, vector<int>& vSPTParentEdge, vector<vector<int> >& vSPT,vector<bitset<KIND> > &vSPTBitS);
    int FindNN(int deviation,int t,bitset<KIND> &vPathBit);
	void FindRepeatedPath(vector<vector<int> >& vvPath);
	int PruneRepeatedPoiPath(vector<int>& bestpath, vector<int>& bestpoi);
	int PruneRepeatedPoiPath(vector<int> &bestpath);
    bitset<KIND> getempPathBit(vector<int> vector1);

    int getTempPathDistance(vector<int> vector1);


    /**
     * Grid info
     */
    vector<vector<int>> nodeMap;
    struct gridNode{
        bitset<KIND> gridNodeBits;
        vector<int> POI;
    };
    vector<int> POIGridMap;//存储每一个点所处的网格id
    vector<gridNode> gridNodeMap;
    int colSize;	// # of columns
    int rowSize;	// # of rows
    int colLength;
    int rowLength;
    MBR *extend;
    vector<int> cellContain;
    void createGridIndex();
    void buildIndex();
    void calculateGridSize(const int i);
    void get_POI_Straight_From_StoT(int s, int t, vector<int> &cell);
    void getPOICoordinateByGrid(int s,int t,vector<int> Query,vector<int> &POICandidate);
    int getCell(Point &point);
    int rightOf(int cid); // return the right cell of given cell
    int upperOf(int cid); // return the upper cell of given cell
    int belowOf(int cid); // return the below cell of given cell
    int leftOf(int cid);

    /**
     * IG-Tree info
     */
    typedef struct{
        bitset<KIND> Gkeyword;
        //bitset<20> Gkeyword;
        vector<int> borders;
        map<int,vector<int>> borderIGNodeMap;
        vector<int> children;
        bool isleaf;
        vector<int> leafnodes;
        int father;
        int level;
// ----- min dis -----
        vector<int> union_borders; // for non leaf node
        vector<int> mind; // min dis, row by row of union_borders
// ----- for pre query init, OCCURENCE LIST in paper -----
        vector<int> nonleafinvlist;
        vector<int> leafinvlist;
        vector<int> up_pos;
        vector<int> current_pos;
        map<int,vector<int>> keywordsNodeMap;
        vector<int> POINode;
        vector<double> kwDisMatrix;
    }TreeNode;

    vector<vector<int>> gtreepath;
    vector<bool> isborder;

    vector<TreeNode> GTree;
    idx_t nvtxs; // |vertices|
    idx_t ncon; // number of weight per vertex
    idx_t* xadj; // array of adjacency of indices
    idx_t* adjncy; // array of adjacency nodes
    idx_t* vwgt; // array of weight of nodes
    idx_t* adjwgt; // array of weight of edges in adjncy
    idx_t nparts; // number of parts to partition
    idx_t objval; // edge cut for partitioning solution
    idx_t* part; // array of partition vector
    idx_t options[METIS_NOPTIONS]; // option array
    bitset<KIND> setIGTreeNodeKeyword(int k);
    vector<int> setIGTreePOINode(int k);
    map<int, vector<int>> setIGTreeBorderNodeMap(int k);
    vector<double> setIGTreeNodeKWDisMatrix(int k);
    void Init_IGTree_Index(const char *IGTreeIndex, const char *IgTreeMIND, const char *IGTreePath);
    void Load_IGTree_MIND(const char *IGTreeMIND);

    int getLowestCommonAncestor(int s, int t);
    void  getSAT(int LCA,vector<int> &satTemp);
    void getSATPOI(int LCA,vector<int> & POIList, bitset<KIND> &uncover);
    void getPOICoordinate(int s, int t, vector<int> queryKey, vector<int> &POICoordinate);
    void getKPathByIGTree(int s,int t,int K,vector<int> POICandidate,vector<vector<int>> &vkPath,vector<int> &vkPathDis);
    int buildGreedyPath(int startNode,int endNode,vector<int> POICandidate,vector<int> &Gpath,int &LengthBound);
    int generateMultPath(int s,int t,int pathNum,vector<pair<int,vector<int>>> &allPath,vector<int> greedyPath,map<int,vector<int>> keyNodeMap,vector<vector<int>> &vPaths,vector<int> &vPathsDis);


    vector<int> POIObjects;
    vector<Point> PointList;
    int extendGridNum;//num of extend grid num to find all keywords
    void GreedyFiler(int s, int t,  vector<int> &POICandidate, vector<int> &POICandidate_Greedy_Filter);
    void TriangleFiltering(int s,int t,vector<int> &POICandidate, vector<int> &POICandidate_Triangle_Filter);
    void getQueryNearOD(int s,int t,vector<int> POICandidate, vector<int> &closeS, vector<int> &closeT );
    void getPPathSet(int s,int t,benchmark::heap<2, int, int> qPath,int lengthThreshold,int GLength,
                     vector<vector<int>> &vvNodePath,vector<int> &vvPPathLength,vector<bitset<KIND>> &vvPPathBit);
    void RangeSearch(bitset<KIND> uncover,int rangeLength,int POINode,vector<int> &inRangeKeywords);
    void expendGrid(bitset<KIND> uncover,int rangeLength,vector<int> &currentGridList,vector<int> &newExpandGridList,int flag,vector<int> &newGridList);
    void RangeSearchByGrid(bitset<KIND> uncover, int rangeLength,int POINode, vector<int> &inRangeKeywords);
    int getPartialPathSet(int s,int t,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength);
    int getPartialPathSet_NoSPT(int s,int t,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength);
    int getMaxLB(vector<int> &vpath, bitset<KIND> &vpathBit,map<int,vector<int>> keyNodeMap, int pos,int ID2);
    int test;
    //verison 1:queue save path Length
    void SY_SKORP_V1(int s,int t,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                  int &resultPathLength,vector<int> &resultPath);

    //version 2:queue save path Length/cover kw num
    void SY_SKORP_V2(int s,int t,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                     int &resultPathLength,vector<int> &resultPath);
    //version 3:get new node consider OD distance
    void SY_SKORP_V3(int s,int t,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                     int &resultPathLength,vector<int> &resultPath);
    void Top_k_SKORP_V1(int s,int t,int k,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                     vector<int> &resultPathLength,vector<vector<int>> &resultPath);
//    void expendGrid(Grid& grid,bitset<KIND> uncover,int rangeLength,vector<int> &currentGridList,vector<int> &newExpandGridList,int flag,vector<int> &newGridList);


    void Top_k_SKORP_Plus(int s,int t,int K,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                     vector<int> &resultPathLength,vector<vector<int>> &resultPath);

    //Range Prune
    vector<int> RPOI;
    long long initPairSum;
    int RangePOI(int s, int t, int dist,int RangeParameter);
    bool testDominated(vector<int> &nodePathTable,
                       vector<vector<int>> pathList,
                       vector<bitset<KIND>> pathListBit,
                              vector<int> pathListLength,
                              vector<int> newPath,
                              bitset<KIND> newPathBit,
                              int newPathLength);
    void Top_K_SKORP_LOAD(int s,int t,int K,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                                vector<int> &resultPathLength,vector<vector<int>> &resultPath);
    void CBIE_NO_improve(int s,int t,int K,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                         vector<int> &resultPathLength,vector<vector<int>> &resultPath);
    int determineInsertPosition(vector<int> path,int node);
    int H2HPath_SY(int ID1, int ID2,vector<int> &plusNode,bitset<KIND> uncoverKW,vector<int> &spPath);
    void TWE(int ID1,int ID2,vector<int> QueryList,vector<vector<int>> &kPaths, vector<int> &kPathLength,ofstream &ofst);
    void TWE_Base(int ID1,int ID2,vector<int> QueryList,vector<vector<int>> &kPaths, vector<int> &kPathLength,ofstream &ofst);
    void getRangePairPath(map<int,vector<int>> &nodeStartMap,
                          map<int,vector<int>> &nodeEndMap,
                          vector<int> &pairPathPosition,
                          benchmark::heap<2, int, int> &qPath,
                          vector<pair<int,int>> &nodeST,
                          int subgraphID,int poi_node,
                          int aveDis,
                          bitset<SIZE1> &pairPathLengthBit,
                          bitset<SIZE1> &subgraphBit,
                          vector<vector<int>> &vvNodePath,
                          vector<int> &vvNodePathLength,
                          vector<bitset<KIND>> &vvNodePathBit);
    vector<int> CPOI_list;
    int s,t;
    int getMaxLB(vector<int> currentPath,int location,vector<int> IG_poiList,bitset<KIND> current_pathBit,int ID1,int ID2);
    void KTWE(int ID1,int ID2,int K,vector<int> QueryList,vector<vector<int>> &kPaths, vector<int> &kPathLength,ofstream &ofst);
};
#endif
