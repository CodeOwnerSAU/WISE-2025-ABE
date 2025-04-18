//
// Created by xiongxing on 3/8/24.
//
#include "graph.h"
int Graph::LCAQuery(int _p, int _q){
    int p = toRMQ[_p], q = toRMQ[_q];

    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;

    int i = LOGD[len], k = LOG2[len];

    q = q - i + 1;
    if (height[RMQIndex[k][p]] < height[RMQIndex[k][q]])
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}
int Graph::distanceQuery(int p, int q){
    if (p == q) return 0;
    int x = belong[p], y = belong[q];
    int lca = LCAQuery(x, y);
    if (lca == x || lca == y){
        queryCnt++;
        //		aCnt++;
        if (lca == y){//x是y的祖先节点
            int v = y;
            y = x;
            x = v;
            v = p;
            p = q;
            q = v;
        }
        return dis[y][pos[x][posSize[x] - 1]];
    }
    else {
        int res = infinity;
        int *dx = dis[x], *dy = dis[y],*p2 = pos2[lca];
        _mm_prefetch(dx, _MM_HINT_T0);//从地址P处预取尺寸为cache line大小的数据缓存
        _mm_prefetch(dy, _MM_HINT_T0);
        _mm_prefetch(p2, _MM_HINT_T0);
        int ps = pos2Size[lca];
        for (int i = 0; i < ps; i++){
            queryCnt ++;
            int tmp = dx[p2[i]] + dy[p2[i]];//dis+dis
            if (res > tmp)
                res = tmp;
        }
        return res;
    }
}
void Graph::readGraph(const char *filename){
    ifstream file(filename);
    //FILE * file = fopen(filename, "r");
    int n, m;
    file>>n>>m;
    //fscanf(file, "%d %d", &n, &m);
    Degree = (int*)malloc(sizeof(int) * (n + 1));
    vector< vector<pair<int, int> > > nb;
    vector<pair<int, int> > v;
    v.clear();
    for (int i = 0; i <= n; i++){
        //	Degree[i] = 0;
        nb.push_back(v);
    }
    //	cout << n << " " << m << endl;
    char c;int x, y, z;
    while (file>>c>>x>>y>>z){
        nb[x].push_back(make_pair(y, z));
    }
    Neighbor = (int**)malloc(sizeof(int*) * (n + 1));
    Weight = (int**)malloc(sizeof(int*) * (n + 1));
    for (int i = 1; i <= n; i++){
        Degree[i] = nb[i].size();
        Neighbor[i] = (int*)malloc(sizeof(int) * nb[i].size());
        Weight[i] = (int*)malloc(sizeof(int) * nb[i].size());
        for (int j = 0; j < nb[i].size(); j++){
            Neighbor[i][j] = nb[i][j].first;
            Weight[i][j] = nb[i][j].second;
        }
    }
}
int Graph::shortestPathQuery(int p, int q){
    vector<int> path;
    path.push_back(p);
    int res = 0;
    while (p != q){
        res++;
        int pq = distanceQuery(p, q);
        for (int i = 0; i < Degree[p]; i++){
            int x = Neighbor[p][i];
            path.push_back(x);
            //	int y = Weight[p][i];
            int xq = distanceQuery(x, q);
            if (xq + Weight[p][i] == pq){
                p = x;
                break;
            }
        }
    }
    return res;
}
void Graph::scanIntArray(int *a, int n){
    fread(a, SIZEOFINT, n, fin_Index);
}
int* Graph::scanIntVector(int *a){
    int _n;
    fread(&_n, SIZEOFINT, 1, fin_Index);
    a = (int*)malloc(sizeof(int) * _n);
    scanIntArray(a, _n);
    return a;
}
void Graph::readIndex(const char *file) {
    int tree_height = 0, tree_width = 0, most_sp = 0;
    fin_Index = fopen(file, "rb");
    fread(&nodeNum, SIZEOFINT, 1, fin_Index);
    int ts;
    fread(&ts, SIZEOFINT, 1, fin_Index);
    TreeSize = ts;
    height = (int *) malloc(sizeof(int) * (ts + 1));
    for (int i = 0; i < ts; i++) {
        fread(&height[i], SIZEOFINT, 1, fin_Index);
    }
    belong = (int *) malloc(sizeof(int) * (nodeNum + 1));
    fread(belong, SIZEOFINT, nodeNum + 1, fin_Index);
    toRMQ = (int *) malloc(sizeof(int) * (nodeNum + 1));
    fread(toRMQ, SIZEOFINT, nodeNum + 1, fin_Index);
    int ris;
    fread(&ris, SIZEOFINT, 1, fin_Index);
    fread(&ts, SIZEOFINT, 1, fin_Index);
    EulerSeq = (int *) malloc(sizeof(int) * (ts + 1));
    RMQIndex = (int **) malloc(sizeof(int *) * (ris + 1));
    for (int i = 0; i < ris; i++) {
        RMQIndex[i] = scanIntVector(RMQIndex[i]);
    }
    fread(&root, SIZEOFINT, 1, fin_Index);
    cout << "root: " << root << endl;

    posSize = (int *) malloc(sizeof(int) * (nodeNum + 1));
    pos2Size = (int *) malloc(sizeof(int) * (nodeNum + 1));
    pos = (int **) malloc(sizeof(int *) * (TreeSize));
    pos2 = (int **) malloc(sizeof(int *) * (TreeSize));
    dis = (int **) malloc(sizeof(int *) * (TreeSize));
    chSize = (int *) malloc(sizeof(int) * (TreeSize));
    ch = (int **) malloc(sizeof(int *) * (TreeSize));

    for (int i = 0; i < TreeSize; i++) {
        fread(&chSize[i], SIZEOFINT, 1, fin_Index);
        ch[i] = (int *) malloc(sizeof(int) * chSize[i]);
        for (int j = 0; j < chSize[i]; j++) {
            int x;
            fread(&x, SIZEOFINT, 1, fin_Index);
            ch[i][j] = x;
        }
    }
    for (int i = 0; i < TreeSize; i++) {
        int x;
        fread(&x, SIZEOFINT, 1, fin_Index);
        fread(&posSize[x], SIZEOFINT, 1, fin_Index);
        pos[x] = (int *) malloc(sizeof(int) * (posSize[x] + 1));
        fread(pos[x], SIZEOFINT, posSize[x], fin_Index);
        if (posSize[x] > tree_width)
            tree_width = posSize[x];
        int _n;
        fread(&_n, SIZEOFINT, 1, fin_Index);
        dis[x] = (int *) malloc(sizeof(int) * _n);
        fread(dis[x], SIZEOFINT, _n, fin_Index);
        if (_n > tree_height)
            tree_height = _n;
    }
    printf("dis read finished!\n");
    for (int i = 0; i < TreeSize; i++) {
        int x;
        fread(&x, SIZEOFINT, 1, fin_Index);
        fread(&pos2Size[x], SIZEOFINT, 1, fin_Index);
        pos2[x] = (int *) malloc(sizeof(int) * (pos2Size[x] + 1));
        fread(pos2[x], SIZEOFINT, pos2Size[x], fin_Index);
        if (pos2Size[x] > most_sp)
            most_sp = pos2Size[x];
    }

    fclose(fin_Index);
    //
    printf("tree height: %d\n", tree_height);
    printf("tree width: %d\n", tree_width);
    printf("most search space: %d\n", most_sp);
}
void Graph::Init_H2H_Index(const char *index, const char* graph){
    readIndex(index);
    readGraph(graph);
    cout<<nodeNum<<endl;
    LOG2 = (int*)malloc(sizeof(int) * (nodeNum * 2 + 10));
    LOGD = (int*)malloc(sizeof(int) * (nodeNum * 2 + 10));
    int k = 0, j = 1;
    for (int i = 0; i < nodeNum * 2 + 10; i++){
        if (i > j * 2){
            j *= 2;
            k++;
        }
        LOG2[i] = k;
        LOGD[i] = j;
    }
    cout << "--------load index finished---------!!!" << endl;
}
int Graph::H2HPath(int ID1, int ID2, vector<int>& vPath,vector<bitset<KIND>> &H2HPathBit){
    vPath.clear();
    H2HPathBit.clear();
    int p = ID1+1 ,q = ID2+1;
    int res = distanceQuery(p, q);
    vPath.push_back(p-1);
    H2HPathBit.push_back(NodesBit[p-1]&QueryBit);
    int tn = p;
    while (p != q){
        int pq = distanceQuery(p, q);
        for (int i = 0; i < Degree[p]; i++){
            int x = Neighbor[p][i];
            int xq = distanceQuery(x, q);
            if (xq + Weight[p][i] == pq){
                p = x;
                vPath.push_back(p-1);
                H2HPathBit.push_back(H2HPathBit.back()|(NodesBit[p-1]&QueryBit));
                //vPathEdge.push_back(adjListEdge[tn-1][i].second);
                tn = x;
                break;
            }
        }
    }
    return res;
}
int Graph::H2HPath(int ID1, int ID2, vector<int>& vPath,bitset<KIND> &H2HPathBit){
    vPath.clear();
    H2HPathBit.reset();
    int p = ID1+1 ,q = ID2+1;
    int res = distanceQuery(p, q);
    vPath.push_back(p-1);
    int tn = p;
    while (p != q){
        int pq = distanceQuery(p, q);
        for (int i = 0; i < Degree[p]; i++){
            int x = Neighbor[p][i];
            int xq = distanceQuery(x, q);
            if (xq + Weight[p][i] == pq){
                p = x;
                vPath.push_back(p-1);
                H2HPathBit|=QueryBit&NodesBit[p-1];
                //vPathEdge.push_back(adjListEdge[tn-1][i].second);
                tn = x;
                break;
            }
        }
    }
    return res;
}
int Graph::H2HPathPlusBit(int ID1, int ID2, vector<int>& vPath,vector<bitset<KIND>> &H2HPathBit,int& pBegin,int &pEnd,int pathCount){
    vPath.clear();
    H2HPathBit.clear();
    pBegin=-1;
    pEnd=-1;
    int p = ID1+1 ,q = ID2+1;
    int res = distanceQuery(p, q);
    vPath.push_back(p-1);
    H2HPathBit.push_back(NodesBit[p-1]&QueryBit);
    int tn = p;
    int index=0;
    int flag=0;
    while (p != q){
        int pq = distanceQuery(p, q);
        for (int i = 0; i < Degree[p]; i++){
            int x = Neighbor[p][i];
            int xq = distanceQuery(x, q);
            if (xq + Weight[p][i] == pq){
                p = x;
                vPath.push_back(p-1);
                index++;
                if(H2HPathBit.back().count()==0&&(H2HPathBit.back()|(NodesBit[p-1]&QueryBit)).count()!=0){
                    pBegin=index;
                }
                H2HPathBit.push_back(H2HPathBit.back()|(NodesBit[p-1]&QueryBit));
                if(H2HPathBit.back().count()==pathCount&&flag==0){
                    pEnd=index;
                    flag=1;
                }
                //vPathEdge.push_back(adjListEdge[tn-1][i].second);
                tn = x;
                break;
            }
        }
    }
    return res;
}

int Graph::getTempPathDistance(vector<int> tempPath) {
    if(tempPath.empty())
        return 0;
    int result=0;
    for(int i=0;i<tempPath.size()-1;i++){
        result+= distanceQuery(tempPath[i]+1,tempPath[i+1]+1);
    }
    return result;
}

int Graph::H2HPath_SY(int ID1, int ID2,vector<int> &plusNode,bitset<KIND> uncoverKW,vector<int> &spPath){
    spPath.clear();
    int p = ID1+1 ,q = ID2+1;
    int res = distanceQuery(p, q);
    int tn = p;
    while (p != q){
        int pq = distanceQuery(p, q);
        for (int i = 0; i < Degree[p]; i++){
            int x = Neighbor[p][i];
            int xq = distanceQuery(x, q);
            if (xq + Weight[p][i] == pq){
                p = x;
                spPath.push_back(p-1);
                if((NodesBit[p-1]&uncoverKW).count()>0&&p-1!=ID1){
                    plusNode.push_back(p-1);
                }
                //vPathEdge.push_back(adjListEdge[tn-1][i].second);
                tn = x;
                break;
            }
        }
    }
    return res;
}


