#include "graph.h"
void Graph::Init_IGTree_Index(const char *IGTreeIndex, const char *IgTreeMIND, const char *IGTreePath){
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // _RB
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // _VOL
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // _RM
    options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM; // _GROW _EDGE _NODE
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // _GREEDY _SEP2SIDED _SEP1SIDED
    // options[METIS_OPTION_NCUTS] = 1;
    // options[METIS_OPTION_NITER] = 10;
    /* balance factor, used to be 500 */
    options[METIS_OPTION_UFACTOR] = 500;
    // options[METIS_OPTION_MINCONN];
    options[METIS_OPTION_CONTIG] = 1;
    // options[METIS_OPTION_SEED];
    options[METIS_OPTION_NUMBERING] = 0;
    // options[METIS_OPTION_DBGLVL] = 0;
    FILE *fin= fopen(IGTreeIndex,"rb");
    int *buf=new int[this->nodeNum];
    int count_borders, count_children, count_leafnodes;
    bool isleaf;
    int father;
    GTree.clear();
    while( fread( &count_borders, sizeof(int), 1, fin ) ){
        TreeNode tn;
        // borders
        tn.borders.clear();
        fread( buf, sizeof(int), count_borders, fin );
        for ( int i = 0; i < count_borders; i++ ){
            tn.borders.push_back(buf[i]);
        }
        // children
        fread( &count_children, sizeof(int), 1, fin );
        fread( buf, sizeof(int), count_children, fin );
        for ( int i = 0; i < count_children; i++ ){
            tn.children.push_back(buf[i]);
        }
        // isleaf
        fread( &isleaf, sizeof(bool), 1, fin );
        tn.isleaf = isleaf;
        // leafnodes
        fread( &count_leafnodes, sizeof(int), 1, fin );
        fread( buf, sizeof(int), count_leafnodes, fin );
        for ( int i = 0; i < count_leafnodes; i++ ){
            tn.leafnodes.push_back(buf[i]);
        }
        // father
        fread( &father, sizeof(int), 1, fin );
        tn.father = father;
        GTree.push_back(tn);
    }
    fclose(fin);

    int count;
    fin = fopen( IGTreePath, "rb" );
    int pos = 0;
    gtreepath.resize(nodeNum);
    while( fread( &count, sizeof(int), 1, fin ) ){
        fread( buf, sizeof(int), count, fin );
        // clear gtreepath
        gtreepath[pos].clear();
        //Nodes[pos].gtreepath.clear();
        for ( int i = 0; i < count; i++ ){
            gtreepath[pos].push_back( buf[i] );
        }
        // pos increase
        pos ++;
    }
    fclose(fin);
    delete[] buf;
    setIGTreeNodeKeyword(0);
    setIGTreePOINode(0);
    for(int i=0;i<GTree.size();i++){
        if(GTree[i].isleaf){
            for(auto border_node:GTree[i].borders){
                for(auto edge:adjList[border_node]){
                    if(gtreepath[edge.first].back()!=i){
                        GTree[i].borderIGNodeMap[gtreepath[edge.first].back()].push_back(border_node);
                    }
                }
            }
        }
    }
    //setIGTreeBorderNodeMap(0);
    Load_IGTree_MIND(IgTreeMIND);
    levelMap=collectLeafByLevel();
    cout<<"IG-Tree Index load finished"<<endl;
}
map<int, vector<int>> Graph:: setIGTreeBorderNodeMap(int k){
    vector<pair<int,int>> temp;
    if(GTree[k].isleaf){
        temp.clear();
        for(auto border_node:GTree[k].borders){
            for(auto edge:adjList[border_node]){
                if(gtreepath[edge.first].back()==k){
                    continue;
                }
                if(std::find(GTree[k].borderIGNodeMap[gtreepath[edge.first].back()].begin(),
                             GTree[k].borderIGNodeMap[gtreepath[edge.first].back()].end(),border_node)==
                        GTree[k].borderIGNodeMap[gtreepath[edge.first].back()].end()){
                    temp.emplace_back(gtreepath[edge.first].back(),border_node);
                    GTree[k].borderIGNodeMap[gtreepath[edge.first].back()].push_back(border_node);
                }
            }
        }
        return GTree[k].borderIGNodeMap;
    } else{
        map<int, vector<int>> pairNode;
        for(auto children:GTree[k].children){
            pairNode=setIGTreeBorderNodeMap(children);
            for(auto t:pairNode){
                   for(auto b_node:t.second){
                       if(std::find(GTree[k].borderIGNodeMap[t.first].begin(),
                                    GTree[k].borderIGNodeMap[t.first].end(),b_node)
                                    ==GTree[k].borderIGNodeMap[t.first].end()){
                       GTree[k].borderIGNodeMap[t.first].push_back(b_node);
                   }
                }
            }
        }
        return GTree[k].borderIGNodeMap;
    }
}
void Graph::Load_IGTree_MIND(const char *IGTreeMIND){
    FILE* fin = fopen( IGTreeMIND, "rb" );
    int* buf;
    int count, pos = 0;
    while( fread( &count, sizeof(int), 1, fin ) ){
        // union borders
        buf = new int[count];
        fread( buf, sizeof(int), count, fin );
        GTree[pos].union_borders.clear();
        for ( int i = 0; i < count; i++ ){
            GTree[pos].union_borders.push_back(buf[i]);
        }
        delete[] buf;
        // mind
        fread( &count, sizeof(int), 1, fin );
        buf = new int[count];
        fread( buf, sizeof(int), count, fin );
        GTree[pos].mind.clear();
        for ( int i = 0; i < count; i++ ){
            GTree[pos].mind.push_back(buf[i]);
        }
        pos++;
        delete[] buf;
    }
    fclose(fin);
}
vector<double> Graph::setIGTreeNodeKWDisMatrix() {
    vector<int> level_order;
    for(auto level:levelMap){
        level_order.push_back(level.first);
    }
    level_order.erase(level_order.begin(),level_order.end());
    for(auto level:level_order){
        for(auto g_node:levelMap[level]){
            GTree[g_node].kwDisMatrix.resize(KIND*(KIND+1)/2);
            if(GTree[g_node].isleaf){
                //GTree[g_node].kwDisMatrix.resize(KIND*(KIND+1)/2);
                for(auto first_kw:GTree[g_node].keywordsNodeMap){
                    vector<int> first_poi=first_kw.second;
                    for(auto second_kw:GTree[g_node].keywordsNodeMap){
                       if(first_kw.first!=second_kw.first){
                           vector<int> second_poi=second_kw.second;
                           int totalDis=0;
                           int pairNum=0;
                           for(auto node1:first_poi){
                               for(auto node2:second_poi){
                                   if(node1!=node2){
                                       totalDis+= distanceQuery(node1+1,node2+1);
                                       pairNum++;
                                   }
                               }
                           }
                           double aveDis=totalDis/(pairNum*1.0);
                           GTree[g_node].kwDisMatrix[first_kw.first*(2*KIND-first_kw.first-1)/2+abs(first_kw.first-second_kw.first)]=aveDis;
                       }
                    }
                }
            } else{
                int leftNodeID=GTree[g_node].children[0];
                int rightNodeID=GTree[g_node].children[1];
                for(int i=0;i<GTree[g_node].kwDisMatrix.size();i++){
                    GTree[g_node].kwDisMatrix[i]=(GTree[leftNodeID].kwDisMatrix[i]+GTree[rightNodeID].kwDisMatrix[i])/2;
                }
                for(auto first_kw:GTree[leftNodeID].keywordsNodeMap){
                    vector<int> first_poi=first_kw.second;
                    for(auto second_kw:GTree[rightNodeID].keywordsNodeMap){
                        vector<int> second_poi=second_kw.second;
                        int totalDis=0;
                        int pairNum=0;
                        for(auto node1:first_poi) {
                            for (auto node2: second_poi) {
                                if (node1 != node2) {
                                    totalDis += distanceQuery(node1 + 1, node2 + 1);
                                    pairNum++;
                                }
                            }
                        }
                        double aveDis=totalDis/(pairNum*1.0);
                        GTree[g_node].kwDisMatrix[first_kw.first*(2*KIND-first_kw.first-1)/2+abs(first_kw.first-second_kw.first)]+=aveDis;
                    }
                }
            }
        }
    }
}
vector<int> Graph::setIGTreePOINode(int k){
    if(GTree[k].isleaf){
        for(int i=0;i<GTree[k].leafnodes.size();i++){
            if(NodesBit[GTree[k].leafnodes[i]].count()>0){
                GTree[k].POINode.push_back(GTree[k].leafnodes[i]);
            }
            for(auto node:GTree[k].POINode){
                for(auto kw:Nodes[node].kwList){
                    if(std::find(GTree[k].keywordsNodeMap[kw].begin(), GTree[k].keywordsNodeMap[kw].end(),node)==GTree[k].keywordsNodeMap[kw].end()){
                        GTree[k].keywordsNodeMap[kw].push_back(node);
                    }
                }
            }
        }
        return GTree[k].POINode;
    }else{
        //bitset<20>temp;
        vector<int> temp;
        temp.clear();
        for(int j=0;j<GTree[k].children.size();j++){
            temp=setIGTreePOINode(GTree[k].children[j]);
            GTree[k].POINode.insert(GTree[k].POINode.end(), temp.begin(),temp.end());
            for(auto node:GTree[k].POINode){
                for(auto kw:Nodes[node].kwList){
                    if(std::find(GTree[k].keywordsNodeMap[kw].begin(), GTree[k].keywordsNodeMap[kw].end(),node)==GTree[k].keywordsNodeMap[kw].end()){
                        GTree[k].keywordsNodeMap[kw].push_back(node);
                    }
                }
            }
        }
        temp=GTree[k].POINode;
        return temp;
    }
}
bitset<KIND> Graph::setIGTreeNodeKeyword(int k) {
    if(GTree[k].isleaf){
        for(int i=0;i<GTree[k].leafnodes.size();i++){
            GTree[k].Gkeyword|=Nodes[GTree[k].leafnodes[i]].keyword;
        }
        return GTree[k].Gkeyword;
    }else{
        bitset<KIND>temp;
        temp.reset();
        for(int j=0;j<GTree[k].children.size();j++){
            GTree[k].Gkeyword|= setIGTreeNodeKeyword(GTree[k].children[j]);
        }
        temp|=GTree[k].Gkeyword;
        return temp;
    }
}

int Graph::getLowestCommonAncestor(int s, int t) {
    int pnode=-1;
    for(int i=0,j=0;i<gtreepath[s].size(),j<gtreepath[t].size();++i,++j){
        if(gtreepath[s][i]==gtreepath[t][j])
        {
            pnode = gtreepath[s][i];
        }
    }
    return pnode;
}

/**
 * get POI list use IGTree save in POICoordinate
 * @param s
 * @param t
 * @param queryKey
 * @param POICoordinate
 */
void Graph::getPOICoordinate(int s, int t, vector<int> queryKey, vector<int> &POICoordinate) {
    //get LCA that can cover all keywords
    POICoordinate.clear();
    int LCA=getLowestCommonAncestor(s,t);
    //cout<<LCA<<endl;
    while((QueryBit & GTree[LCA].Gkeyword).count() != QueryBit.count())
    {
        LCA = GTree[LCA].father;
    }
    vector<int> satTemp;
    getSAT(LCA,satTemp);
    vector<int> neset(satTemp.begin(),satTemp.end());
    POICoordinate=neset;
}
/**
 * get all POI point in leafNode in IGTree
 * @param LCA
 * @param satTemp
 */
void Graph::getSAT(int LCA,vector<int> &satTemp) {
//    vector<int>node;
//    node.clear();
    for(auto kw:QueryWord){
        for(auto node:GTree[LCA].keywordsNodeMap[kw]){
            if(std::find(satTemp.begin(), satTemp.end(),node)==satTemp.end()){
                satTemp.push_back(node);
            }
        }
    }
}
void Graph::getSATPOI(int LCA,vector<int> & POIList, bitset<KIND> & uncover) {
    vector<int>node;
    node.clear();
    if(GTree[LCA].isleaf){
        for(int j=0;j<GTree[LCA].POINode.size();j++){
            if((NodesBit[GTree[LCA].POINode[j]]&uncover).count()>0){
                POIList.push_back(GTree[LCA].POINode[j]);
            }
        }
    }
    for(int i=0;i<GTree[LCA].children.size();i++){
        getSATPOI(GTree[LCA].children[i],POIList,uncover);
    }
}

map<int, std::vector<int>> Graph::collectLeafByLevel() {
    std::map<int, std::vector<int>> levelMap;
    std::queue<std::pair<int, int>> q;
    q.push({0, 1});
    while (!q.empty()) {
        auto [currentId, currentLevel] = q.front();
        q.pop();
        const TreeNode &currentNode = GTree[currentId];
        if (!currentNode.children.empty()) {
            levelMap[currentLevel].push_back(currentId);
            for (int childId: currentNode.children) {
                q.push({childId, currentLevel + 1});
            }
        }
    }
    return levelMap;
}

