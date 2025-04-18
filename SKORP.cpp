#include "graph.h"
#include "thread"
void Graph::GreedyFiler(int s, int t,  int GreedyPath,vector<int> &POICandidate, vector<int> &POICandidate_Greedy_Filter) {
    //1 get Greedy path from s and from t
    bitset<KIND> greedyBit;
    POICandidate_Greedy_Filter.clear();
    for(auto node:POICandidate){
        int dis1=distanceQuery(s+1,node+1);
        int dis2=distanceQuery(node+1,t+1);
        if(dis1+dis2<GreedyPath){
            POICandidate_Greedy_Filter.push_back(node);
            greedyBit|=NodesBit[node]&QueryBit;
        }
    }
    if(greedyBit.count()!=QueryBit.count()){
        POICandidate_Greedy_Filter.clear();
        POICandidate_Greedy_Filter=POICandidate;
    }
}
void Graph::TriangleFiltering(int s,int t,vector<int> &POICandidate, vector<int> &POICandidate_Triangle_Filter){
    bitset<KIND> POISetBit;
    for(auto node:POICandidate){
        //EuclideanDistance(s,t);
        double disST= distanceQuery(s+1,t+1);
        double disSNode= distanceQuery(s+1,node+1);
        double disTNode= distanceQuery(node+1,t+1);
//        double disST= EuclideanDistance(s,t);
//        double disSNode= EuclideanDistance(s,node);
//        double disTNode= EuclideanDistance(node,t);
        if(!((disTNode>disST&&disTNode>disSNode)||(disSNode>disST&&disSNode>disTNode))){
            POICandidate_Triangle_Filter.push_back(node);
            POISetBit|=NodesBit[node]&QueryBit;
        }
    }
    if(POICandidate_Triangle_Filter.empty()||POISetBit.count()!=QueryBit.count()){
        POICandidate_Triangle_Filter.clear();
        POICandidate_Triangle_Filter=POICandidate;
    }
}
void Graph::getQueryNearOD(int s,int t,vector<int> POICandidate, vector<int> &closeS, vector<int> &closeT) {
    map<int,vector<int>> keyNodeMap;
    for(int i=0;i<QueryWord.size();i++){
        vector<int> l;
        keyNodeMap.insert({QueryWord[i],l});
    }
    for(auto node:POICandidate){
        for(auto key:keyNodeMap){
            if(NodesBit[node].test(key.first)){
                keyNodeMap[key.first].push_back(node);
            }
        }
    }
    double disST= EuclideanDistance(s,t);
    for(auto keyNode:keyNodeMap){
        int CSDis=0;
        int CTDis=0;
        for(auto node:keyNode.second){
            double disSNode= EuclideanDistance(s,node);
            double disTNode= EuclideanDistance(node,t);
            CSDis+=disSNode;
            CTDis+=disTNode;
//            if((double)(disSNode/disTNode)<0.5){
//                //current kw is closer s
//
//                cout<<node<<" is closer s"<<endl;
//            }
//            if((double)(disTNode/disSNode)<0.5){
//                //current kw is closer t
//                CT++;
//                cout<<node<<" is closer t"<<endl;
//            }
        }
        cout<<CSDis<<endl;
        cout<<CTDis<<endl;
        if((CSDis<CTDis&&(double)(CTDis-CSDis)/CTDis>(double)(1.0/4.0))){
            closeS.push_back(keyNode.first);
        }
        if((CTDis<CSDis&&(double)(CSDis-CTDis)/CSDis>(double)(1.0/4.0))){
            closeT.push_back(keyNode.first);
        }
    }
}
void Graph::getPPathSet(int s,int t,benchmark::heap<2, int, int> qPath,int lengthThreshold,int GLength,
                        vector<vector<int>> &vvNodePath,vector<int> &vvPPathLength,vector<bitset<KIND>> &vvPPathBit){

}
void Graph::RangeSearch(bitset<KIND> uncover,int rangeLength,int POINode,vector<int> &inRangeKeywords){
    vector<int> POINodeInGrid;
    vector<int> haveKWGrid;
    //POI Node id
    int startGridID=POIGridMap[POINode];
    if((gridNodeMap[startGridID].gridNodeBits&uncover).count()==uncover.count()){
        haveKWGrid.push_back(startGridID);
        for(auto node:gridNodeMap[startGridID].POI){
            if((NodesBit[node]&uncover).count()>0){
                inRangeKeywords.push_back(node);
            }
        }
    }
    int upGrid,downGrid,leftGrid,rightGrid;
    long long upDis,downDis,leftDis,rightDis;
    upDis=0;
    downDis=0;
    leftDis=0;
    rightDis=0;
    vector<int> extendGridList;
    vector<int> GridList;
    vector<int> temp;
    temp.push_back(startGridID);
    int flag=0;
    while(1){
        for(int i=0;i<temp.size();i++){
            int currentGrid=temp[i];
            upGrid=upperOf(currentGrid);
            if(std::find(extendGridList.begin(), extendGridList.end(),upGrid)==extendGridList.end()&&upGrid!=startGridID){
                if((gridNodeMap[upGrid].gridNodeBits&uncover).count()>0){
                    if(std::find(haveKWGrid.begin(), haveKWGrid.end(),upGrid)==haveKWGrid.end()){
                        haveKWGrid.push_back(upGrid);
                        for(auto node:gridNodeMap[upGrid].POI){
                            if((NodesBit[node]&uncover).count()>0){
                                inRangeKeywords.push_back(node);
                            }
                        }
                    }
                }
                extendGridList.push_back(upGrid);
                upDis+=colLength;
            }
            downGrid=belowOf(currentGrid);
            if(std::find(extendGridList.begin(), extendGridList.end(),downGrid)==extendGridList.end()&&downGrid!=startGridID){
                if((gridNodeMap[downGrid].gridNodeBits&uncover).count()>0){
                    if(std::find(haveKWGrid.begin(), haveKWGrid.end(),downGrid)==haveKWGrid.end()){
                        haveKWGrid.push_back(downGrid);
                        for(auto node:gridNodeMap[downGrid].POI){
                            if((NodesBit[node]&uncover).count()>0){
                                inRangeKeywords.push_back(node);
                            }
                        }
                    }
                }
                extendGridList.push_back(downGrid);
                downDis+=colLength;
            }
            leftGrid=leftOf(currentGrid);
            if(std::find(extendGridList.begin(), extendGridList.end(),leftGrid)==extendGridList.end()&&leftGrid!=startGridID){
                if((gridNodeMap[leftGrid].gridNodeBits&uncover).count()>0){
                    if(std::find(haveKWGrid.begin(), haveKWGrid.end(),leftGrid)==haveKWGrid.end()){
                        haveKWGrid.push_back(leftGrid);
                        for(auto node:gridNodeMap[leftGrid].POI){
                            if((NodesBit[node]&uncover).count()>0){
                                inRangeKeywords.push_back(node);
                            }
                        }
                    }
                }
                extendGridList.push_back(leftGrid);
                leftDis+=rowLength;
            }
            rightGrid=rightOf(currentGrid);
            if(std::find(extendGridList.begin(), extendGridList.end(),rightGrid)==extendGridList.end()&&rightGrid!=startGridID){
                if((gridNodeMap[rightGrid].gridNodeBits&uncover).count()>0){
                    if(std::find(haveKWGrid.begin(), haveKWGrid.end(),rightGrid)==haveKWGrid.end()){
                        haveKWGrid.push_back(rightGrid);
                        for(auto node:gridNodeMap[rightGrid].POI){
                            if((NodesBit[node]&uncover).count()>0){
                                inRangeKeywords.push_back(node);
                            }
                        }
                    }
                }
                extendGridList.push_back(rightGrid);
                rightDis+=rowLength;
            }
            if((upDis>rangeLength*1.5)&&(downDis>rangeLength*1.5)&&(leftDis>rangeLength*1.5)&&(rightDis>rangeLength*1.5)){
                flag=1;
                break;
            }
        }
        if(flag==1){
            break;
        }
        sort(extendGridList.begin(),extendGridList.end());
        //delete same cell
        extendGridList.erase(unique(extendGridList.begin(),extendGridList.end()),extendGridList.end());
        temp.clear();
        temp.insert(temp.begin(),extendGridList.begin(),extendGridList.end());
    }
}
void Graph::expendGrid(bitset<KIND> uncover,int rangeLength,vector<int> &currentGridList,vector<int> &newExpandGridList,int flag,vector<int> &newGridList){
    newGridList.clear();
    if(flag==1){
        for(auto cell:currentGridList){
            newExpandGridList.push_back(upperOf(cell));
            newGridList.push_back(upperOf(cell));
        }
    } else if(flag==2){

        for(auto cell:currentGridList){
            newExpandGridList.push_back(belowOf(cell));
            newGridList.push_back(belowOf(cell));
        }
    } else if(flag==3){
        for(auto cell:currentGridList){
            newExpandGridList.push_back(cell-1);
            newGridList.push_back(cell-1);
        }
    }else if(flag==4){
        for(auto cell:currentGridList){
            newExpandGridList.push_back(cell+1);
            newGridList.push_back(cell+1);
        }
    }
}
void Graph::RangeSearchByGrid(bitset<KIND> uncover, int rangeLength,int POINode, vector<int> &inRangeKeywords){
    vector<int> POINodeInGrid;
    vector<int> haveKWGrid;
    int startGridID=POIGridMap[POINode];
    if((gridNodeMap[startGridID].gridNodeBits&uncover).count()>0){
        haveKWGrid.push_back(startGridID);
        for(auto node:gridNodeMap[startGridID].POI){
            if((NodesBit[node]&uncover).count()>0){
                inRangeKeywords.push_back(node);
            }
        }
    }
    int extendCount=1;
    vector<int> currentExpendGridList;
    vector<int> upGrid;
    vector<int> downGrid;
    vector<int> leftGrid;
    vector<int> rightGrid;
    int initUpGrid=upperOf(startGridID);
    upGrid.push_back(initUpGrid);
    upGrid.push_back(initUpGrid-1);
    upGrid.push_back(initUpGrid+1);
    currentExpendGridList.insert(currentExpendGridList.end(),upGrid.begin(),upGrid.end());
    int initDownGrid=belowOf(startGridID);
    downGrid.push_back(initDownGrid);
    downGrid.push_back(initDownGrid-1);
    downGrid.push_back(initDownGrid+1);
    currentExpendGridList.insert(currentExpendGridList.end(),downGrid.begin(),downGrid.end());
    int initLeftGrid=leftOf(startGridID);
    leftGrid.push_back(initLeftGrid);
    leftGrid.push_back(upperOf(initLeftGrid));
    leftGrid.push_back(belowOf(initLeftGrid));
    currentExpendGridList.insert(currentExpendGridList.end(),leftGrid.begin(),leftGrid.end());
    int initRightGrid=rightOf(startGridID);
    rightGrid.push_back(initRightGrid);
    rightGrid.push_back(upperOf(initRightGrid));
    rightGrid.push_back(belowOf(initDownGrid));
    int upLeftGrid=0,upRightGrid=0,downLeftGrid=0,downRightGrid=0;
    upLeftGrid=upGrid[0]-1;
    upRightGrid=upGrid[0]+1;
    downLeftGrid=downGrid[0]-1;
    downRightGrid=downGrid[0]+1;
    int couLeft=upLeftGrid-startGridID;
    int couRight=upRightGrid-startGridID;
    currentExpendGridList.insert(currentExpendGridList.end(),rightGrid.begin(),rightGrid.end());
    std::sort(currentExpendGridList.begin(), currentExpendGridList.end());
    currentExpendGridList.erase(unique(currentExpendGridList.begin(),currentExpendGridList.end()),currentExpendGridList.end());
    long long upDis=0,downDis=0,leftDis=0,rightDis=0;
    //get all POI hat can serve kw
    for(auto cell:currentExpendGridList){
        if((gridNodeMap[cell].gridNodeBits&uncover).count()>0){
            for(auto node:gridNodeMap[cell].POI){
                if((NodesBit[node]&uncover).count()>0){
                    inRangeKeywords.push_back(node);
                }
            }
        }
    }
    upDis+=colLength;
    leftDis+=rowLength;
    int count=0;
    while(upDis<(int)rangeLength*1.5||leftDis<(int)rangeLength*1.5){
        vector<int> newGrid;
        currentExpendGridList.clear();
        expendGrid(uncover,rangeLength,upGrid,currentExpendGridList,1,newGrid);
        upGrid=newGrid;
        expendGrid(uncover,rangeLength,downGrid,currentExpendGridList,2,newGrid);
        downGrid=newGrid;
        expendGrid(uncover,rangeLength,leftGrid,currentExpendGridList,3,newGrid);
        leftGrid=newGrid;
        expendGrid(uncover,rangeLength,rightGrid,currentExpendGridList,4,newGrid);
        rightGrid=newGrid;
        extendCount++;
        currentExpendGridList.push_back(startGridID+couLeft*extendCount);
        currentExpendGridList.push_back(startGridID-couLeft*extendCount);
        upGrid.push_back(startGridID+couLeft*extendCount);
        upGrid.push_back(startGridID-couLeft*extendCount);
        currentExpendGridList.push_back(startGridID+couRight*extendCount);
        currentExpendGridList.push_back(startGridID-couRight*extendCount);
        downGrid.push_back(startGridID+couRight*extendCount);
        downGrid.push_back(startGridID-couRight*extendCount);
        leftGrid.push_back(startGridID-couLeft*extendCount);
        leftGrid.push_back(startGridID-couRight*extendCount);
        rightGrid.push_back(startGridID+couLeft*extendCount);
        rightGrid.push_back(startGridID+couRight*extendCount);
        upDis+=colLength;
        leftDis+=rowLength;
        for(auto cell:currentExpendGridList){
            if((gridNodeMap[cell].gridNodeBits&uncover).count()>0){
                for(auto node:gridNodeMap[cell].POI){
                    if((NodesBit[node]&uncover).count()>0){
                        inRangeKeywords.push_back(node);
                    }
                }
            }
        }
    }
}
int Graph::getPartialPathSet(int s,int t,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength) {
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span1(0.0);
    std::chrono::duration<double> time_span2(0.0);
    t1 = std::chrono::high_resolution_clock::now();
    int OD= distanceQuery(s+1,t+1);
    vector<int> vSPTDistanceS(nodeNum, INF); //record dist from node to root ID2
    vector<int> vSPTDistanceT(nodeNum, INF); //record dist from node to root ID2
    vector<int> vSPTParentS(nodeNum, -1);  //record node's prenode in shorest path
    vector<int> vSPTParentT(nodeNum, -1);  //record node's prenode in shorest path
    vector<int> vSPTParentEdgeS(nodeNum, -1);
    vector<int> vSPTParentEdgeT(nodeNum, -1);
    vector<int> vTmpS;
    vector<int> vTmpT;
    vector<vector<int> > vSPTS(nodeNum, vTmpS); //Tree from root  ID2
    vector<vector<int> > vSPTT(nodeNum, vTmpT); //Tree from root  ID2
    vector<bitset<KIND> > vSPTBitS(NodesBit);
    vector<bitset<KIND> > vSPTBitT(NodesBit);
//    SPT(s, vSPTDistanceS, vSPTParentS, vSPTParentEdgeS, vSPTS,vSPTBitS);
//    SPT(t, vSPTDistanceT, vSPTParentT, vSPTParentEdgeT, vSPTT,vSPTBitT);
    //并发构建SPT

//    SPT(s, vSPTDistanceS, vSPTParentS, vSPTParentEdgeS, vSPTS,vSPTBitS);
//    SPT(t, vSPTDistanceT, vSPTParentT, vSPTParentEdgeT, vSPTT,vSPTBitT);
    t1 = std::chrono::high_resolution_clock::now();
    std::thread thread1([&]() { SPT(s, vSPTDistanceS, vSPTParentS, vSPTParentEdgeS, vSPTS, vSPTBitS); });
    std::thread thread2([&]() { SPT(t, vSPTDistanceT, vSPTParentT, vSPTParentEdgeT, vSPTT, vSPTBitT); });
    thread1.join();
    thread2.join();
    t2 = std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    cout<<"SPT time is "<<time_span1.count()<<endl;
    t1 = std::chrono::high_resolution_clock::now();
    //cout<<"SPT Time is"<<time_span.count()<<endl;
    if((vSPTBitS[t]&QueryBit).count()==QueryBit.count()||(vSPTBitT[s]&QueryBit).count()==QueryBit.count()){
        cout<<"SKORP path is "<<vSPTDistanceS[t]<<" "<<endl;
        return vSPTDistanceS[t];
    }
    vector<vector<int>> vvNodePath;
    vector<int> vvPPathLength;
    vector<bitset<KIND>> vvPPathBit;
    map<int,vector<int>> keyNodeMap;
    keyNodeMap.clear();
    t1 = std::chrono::high_resolution_clock::now();
    for(int i=0;i<QueryWord.size();i++){
        vector<int> l;
        keyNodeMap.insert({QueryWord[i],l});
    }
    for(auto node:POICandidate){
        for(auto key:keyNodeMap){
            if(NodesBit[node].test(key.first)){
                keyNodeMap[key.first].push_back(node);
            }
        }
    }
    vector<pair<int,int>> keyNum;
    for(auto map:keyNodeMap){
        keyNum.emplace_back(map.second.size(),map.first);
    }
    std::sort(keyNum.begin(), keyNum.end());
    vector<bool> isCP(nodeNum,false);
    benchmark::heap<2, int, int> qPath(nodeNum*50);
    t2 = std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    cout<<"Map time is "<<time_span1.count()<<endl;
    //benchmark::heap<2,int,int> qPathGrid(nodeNum);
    int count=0;
//    for(int i=0;i<QueryWord.size();i++){
//        if(uncover.test(QueryWord[i])){
//            for(auto node1:keyNodeMap[QueryWord[i]]){
//                bitset<KIND> temp(NodesBit[node1]&QueryBit);
//                vector<int> inRangeKW;
//                RangeSearchByGrid(temp^QueryBit,lengthThreshold,node1,inRangeKW);
//                //RangeSearchByGrid(this,grid,temp^G.QueryBit,lengthThreshold,node1,inRangeKW);
//                if(inRangeKW.size()>0){
//                    for(auto node2:inRangeKW){
//                        vector<int> PPath;
//                        PPath.push_back(node1);
//                        PPath.push_back(node2);
//                        bitset<KIND> pathBit = (NodesBit[node1] & uncover | NodesBit[node2] & uncover);
//                        int dis= distanceQuery(node1+1,node2+1);
//                        if(std::find(vvPPathLength.begin(), vvPPathLength.end(),dis)==vvPPathLength.end()&&dis<lengthThreshold){
//                            vvNodePath.push_back(PPath);
//                            vvPPathBit.push_back(pathBit);
//                            vvPPathLength.push_back(dis);
//                            int dis1=dis+vSPTDistanceS[PPath[0]]+vSPTDistanceT[PPath[PPath.size()-1]]-OD;
//                            int dis2=dis+vSPTDistanceT[PPath[0]]+vSPTDistanceS[PPath[PPath.size()-1]]-OD;
//                            int smallDis=(dis1<dis2)?dis1:dis2;
//                            qPathGrid.update(vvNodePath.size()-1,smallDis);
//                            count++;
//                        }
////                            cout<<node1<<endl;
////                            cout<<"have kw"<<endl;
//
//                    }
//                }
//            }
//        }
//    }
//    cout<<count<<endl;
//    t2 = std::chrono::high_resolution_clock::now();
//    time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    t1 = std::chrono::high_resolution_clock::now();
    int topPathID=0;
    int topPathDis=0;
    //qPathGrid.extract_min(topPathID ,topPathDis);
    int pNum=0;
    vector<int> H2Hpath;
    vector<int> H2HEdge;
    vector<bitset<KIND>> H2HPathBit;
    //H2HPath(s,t,H2Hpath,H2HEdge,H2HPathBit);
    for(int i=0;i<QueryWord.size();i++) {
        for (int j = i + 1; j < QueryWord.size(); j++) {
            if (uncover.test(QueryWord[i]) && uncover.test(QueryWord[j])) {
                for (auto node1: keyNodeMap[QueryWord[i]]) {
                    int c=1;
                    int POIDis=0;
                    int lastDis=0;
                    for (auto node2: keyNodeMap[QueryWord[j]]) {
                        if(c==1){
                            POIDis = distanceQuery(node1 + 1, node2 + 1);
                            lastDis=POIDis;
                            c++;
                        } else{
                            int Edis= (int)EuclideanDistance(node1,node2);
                            if(Edis>lengthThreshold){
                                continue;
                            } else{
                                if(Edis>lastDis){
                                    continue;
                                } else{
                                    POIDis= distanceQuery(node1+1,node2+1);
                                    lastDis= MAX(lastDis,POIDis);
                                    //lastDis=POIDis;
                                }
                            }
                        }
                        if (POIDis < lengthThreshold && POIDis!= 0) {
                            int id1=POIGridMap[node1];
                            int id2=POIGridMap[node2];
                            vector<int> PPath;
                            PPath.push_back(node1);
                            PPath.push_back(node2);
                            isCP[node2]=true;
                            bitset<KIND> pathBit = (NodesBit[node1] & uncover | NodesBit[node2] & uncover);
                            vvNodePath.push_back(PPath);
                            vvPPathBit.push_back(pathBit);
                            vvPPathLength.push_back(POIDis);
                            int dis1=POIDis+ distanceQuery(PPath[0]+1,s+1)+ distanceQuery(PPath[PPath.size()-1]+1,t+1)-OD;
                            int dis2=POIDis+distanceQuery(PPath[0]+1,t+1)+ distanceQuery(PPath[PPath.size()-1]+1,s+1)-OD;
                            //int dis1=POIDis+vSPTDistanceS[PPath[0]]+vSPTDistanceT[PPath[PPath.size()-1]]-OD;
                            //int dis2=POIDis+vSPTDistanceT[PPath[0]]+vSPTDistanceS[PPath[PPath.size()-1]]-OD;
                            int smallDis=(dis1<dis2)?dis1:dis2;
                            qPath.update(vvNodePath.size()-1,smallDis);
                            //qPath.update(vvNodePath.size()-1,dis);
                            count++;
                            pNum++;
                            //cout<<count<<endl;
                        }
                    }
                }
            }
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span2 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    cout<<"Init time is "<<time_span2.count()<<endl;
//    for(int i=0;i<QueryWord.size();i++) {
//        for (int j = i + 1; j < QueryWord.size(); j++) {
//            if (uncover.test(QueryWord[i]) && uncover.test(QueryWord[j])) {
//                for (auto node1: keyNodeMap[QueryWord[i]]) {
//                    int POIDis=0;
//                    for (auto node2: keyNodeMap[QueryWord[j]]) {
//                            POIDis = distanceQuery(node1 + 1, node2 + 1);
//                        if (POIDis < lengthThreshold && POIDis != 0) {
//                            int id1=POIGridMap[node1];
//                            int id2=POIGridMap[node2];
//                            vector<int> PPath;
//                            PPath.push_back(node1);
//                            PPath.push_back(node2);
//                            isCP[node2]=true;
//                            bitset<KIND> pathBit = (NodesBit[node1] & uncover | NodesBit[node2] & uncover);
//                            vvNodePath.push_back(PPath);
//                            vvPPathBit.push_back(pathBit);
//                            vvPPathLength.push_back(POIDis);
//                            int dis1=POIDis+vSPTDistanceS[PPath[0]]+vSPTDistanceT[PPath[PPath.size()-1]]-OD;
//                            int dis2=POIDis+vSPTDistanceT[PPath[0]]+vSPTDistanceS[PPath[PPath.size()-1]]-OD;
//                            int smallDis=(dis1<dis2)?dis1:dis2;
//                            qPath.update(vvNodePath.size()-1,smallDis);
//                            //qPath.update(vvNodePath.size()-1,dis);
//                            count++;
//                            pNum++;
//                            //cout<<count<<endl;
//                        }
//                    }
//                }
//            }
//        }
//    }
    cout<<pNum<<endl;
    //qPath.extract_min(topPathID,topPathDis);
//
////   bitset<KIND> pathBit=spBit|vvPPathBit[topPathID];
////   cout<<pathBit.count()<<endl;
    int NUM=0;
    int num=0;
    vector<int> resultPath;

   while(!qPath.empty()){
       //cout<<qPath.size()<<endl;
       NUM++;
       //cout<<NUM<<endl;
       qPath.extract_min(topPathID,topPathDis);
       vector<int> currentPath(vvNodePath[topPathID]);
       bitset<KIND> currentPathBit(vvPPathBit[topPathID]);
       bitset<KIND> spBit=QueryBit
               &(vSPTBitS[vSPTParentS[vvNodePath[topPathID][0]]]|vSPTBitT[vSPTParentT[vvNodePath[topPathID][vvNodePath[topPathID].size()-1]]]);
       if((spBit|currentPathBit).count()==QueryBit.count()){
           cout<<"current path is a feasible"<<endl;
           int  startNode=vvNodePath[topPathID][0];
           int endNode=vvNodePath[topPathID][vvNodePath[topPathID].size()-1];
           int newSPDis1=vSPTDistanceS[startNode]+vvPPathLength[topPathID]+
                         vSPTDistanceT[endNode];
           int newSPDis2=vSPTDistanceT[startNode]+vvPPathLength[topPathID]+
                         vSPTDistanceS[endNode];
           resultPath.clear();
           if(newSPDis1<newSPDis2){
               resultPath.push_back(s);
               resultPath.insert(resultPath.end(),vvNodePath[topPathID].begin(),vvNodePath[topPathID].end());
               resultPath.push_back(t);
           } else{
               resultPath.push_back(s);
               std::reverse(currentPath.begin(), currentPath.end());
               resultPath.insert(resultPath.end(),currentPath.begin(),currentPath.end());
               resultPath.push_back(t);
           }
           break;
       }
       bitset<KIND> uncoverKW=currentPathBit^QueryBit;
       int PPathLength1=vSPTParentS[vvNodePath[topPathID][0]]+
               vvPPathLength[topPathID]+vSPTParentT[vvNodePath[topPathID][vvNodePath[topPathID].size()-1]];
       int PPathLength2=vSPTParentT[vvNodePath[topPathID][0]]+
                        vvPPathLength[topPathID]+vSPTParentS[vvNodePath[topPathID][vvNodePath[topPathID].size()-1]];
       if(PPathLength1>GLength||PPathLength2>GLength){
           //cout<<"GLength Pruned"<<endl;
           continue;
       }
       for(auto kw:QueryWord){
           if(uncoverKW.test(kw)){
               for(auto node:keyNodeMap[kw]){
                   num=vvNodePath[topPathID].size();
                   int dis1=distanceQuery(node+1,vvNodePath[topPathID][0]+1)+vvPPathLength[topPathID];
                   int dis2=distanceQuery(vvNodePath[topPathID][vvNodePath[topPathID].size()-1]+1,node+1)+vvPPathLength[topPathID];
                   if(dis1<(lengthThreshold*num)||dis2<(lengthThreshold*num)){
                       int insertStart=0;
                       if(dis1<dis2){
                           insertStart=1;
                       }
                       vector<int> newPPath;
                       bitset<KIND> newPPahBit(vvPPathBit[topPathID]);
                       if(insertStart==1){
                           newPPath.push_back(node);
                           newPPath.insert(newPPath.end(),currentPath.begin(),currentPath.end());
                           vvPPathLength.push_back(dis1);
                       } else{
                           newPPath.insert(newPPath.end(),currentPath.begin(),currentPath.end());
                           newPPath.push_back(node);
                           vvPPathLength.push_back(dis2);
                       }
                       newPPahBit|=NodesBit[node]&uncoverKW;
                       vvNodePath.push_back(newPPath);
                       vvPPathBit.push_back(newPPahBit);
                       int startNode=vvNodePath[vvNodePath.size()-1][0];
                       int endNode=vvNodePath[vvNodePath.size()-1][vvNodePath.back().size()-1];
                       int newSPDis1=vSPTDistanceS[startNode]+vvPPathLength.back()+
                               vSPTParentT[endNode];
                       int newSPDis2=vSPTDistanceT[startNode]+vvPPathLength.back()+
                                     vSPTParentS[endNode];
                       int newSPDis=(newSPDis1<newSPDis2)?newSPDis1:newSPDis2;
                       //cout<<newSPDis<<endl;
                       qPath.update(vvNodePath.size()-1,newSPDis-OD);
                   }
               }
           }
       }
   }
   int Sum=0;

   if(!resultPath.empty()){
       PruneRepeatedPoiPath(resultPath);
       for(int i=0;i<resultPath.size()-1;i++){
           //cout<<resultPath[i]<<" "<<resultPath[i+1]<<endl;
           Sum+= distanceQuery(resultPath[i]+1,resultPath[i+1]+1);
       }
   } else{
       Sum=GLength;
       cout<<"error"<<endl;
   }
   cout<<"SKORP path is "<<Sum<<" "<<endl;
    return Sum;
    //cout<<topPathID<<" "<<topPathDis<<endl;
}
int Graph::getPartialPathSet_NoSPT(int s,int t,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength){
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span1(0.0);
    std::chrono::duration<double> time_span2(0.0);
    t1 = std::chrono::high_resolution_clock::now();
    int OD= distanceQuery(s+1,t+1);
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    vector<int> H2Hpath;
    vector<int> H2HEdge;
    vector<bitset<KIND>> H2HPathBit;
    H2HPath(s,t,H2Hpath,H2HPathBit);
    //cout<<"SPT Time is"<<time_span.count()<<endl;
    if((H2HPathBit.back()&QueryBit).count()==QueryBit.count()){
        //cout<<"shortest path is feasible path"<<endl;
        return OD;
    }
    vector<vector<int>> vvNodePath;
    vector<int> vvPPathLength;
    vector<bitset<KIND>> vvPPathBit;
    map<int,vector<int>> keyNodeMap;
    keyNodeMap.clear();
    t1 = std::chrono::high_resolution_clock::now();
    for(int i=0;i<QueryWord.size();i++){
        vector<int> l;
        keyNodeMap.insert({QueryWord[i],l});
    }
    for(auto node:POICandidate){
        for(auto key:keyNodeMap){
            if(NodesBit[node].test(key.first)){
                keyNodeMap[key.first].push_back(node);
            }
        }
    }
    vector<bool> isCP(nodeNum,false);
    benchmark::heap<2, int, int> qPath(nodeNum*50);
    t2 = std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    //cout<<"Map time is "<<time_span1.count()<<endl;
    //benchmark::heap<2,int,int> qPathGrid(nodeNum);
    int count=0;
    t1 = std::chrono::high_resolution_clock::now();
    int topPathID=0;
    int topPathDis=0;
    //qPathGrid.extract_min(topPathID ,topPathDis);
    int pNum=0;
    //H2HPath(s,t,H2Hpath,H2HEdge,H2HPathBit);

    /**
     * Init path by IG-Tree
     */
//    vector<bool> isDo(POICandidate.size(),false);
//    for(int i=0;i<POICandidate.size();i++){
//        if(!isDo[i]){
//            int GTreeNode=gtreepath[POICandidate[i]].back();
//            auto nodePoint = std::find(GTree[GTreeNode].leafnodes.begin(),
//
//                                       GTree[GTreeNode].leafnodes.end(), POICandidate[i]);
//            int minDis=INT_MAX;
//            for(auto bnode:GTree[GTreeNode].union_borders){
//                minDis= min(distanceQuery(bnode+1,POICandidate[i]+1),minDis);
//            }
//
//            if(minDis>lengthThreshold){
//                cout<<"1"<<endl;
//                //get all POI node in current Node;
//            }else{
//                GTreeNode=GTree[GTreeNode].father;
//            }
////            int index= std::distance(GTree[GTreeNode].leafnodes.begin(),nodePoint);
////            int cindex=GTree[GTreeNode].leafnodes.size()*(index-1);
////            vector<int> minds(GTree[GTreeNode].mind.begin()+cindex,GTree[GTreeNode].mind.begin()+cindex+GTree[GTreeNode].union_borders.size());
//
//        }
//    }


    for(int i=0;i<QueryWord.size();i++) {
        for (int j = i + 1; j < QueryWord.size(); j++) {
            if (uncover.test(QueryWord[i]) && uncover.test(QueryWord[j])) {
                for (auto node1: keyNodeMap[QueryWord[i]]) {
                    int c=1;
                    int POIDis=0;
                    int lastDis=0;
                    for (auto node2: keyNodeMap[QueryWord[j]]) {
                        if(c==1){
                            POIDis = distanceQuery(node1 + 1, node2 + 1);
                            lastDis=POIDis;
                            c++;
                        } else{
                            int Edis= (int)EuclideanDistance(node1,node2);
                            Edis=Edis/10;
                            if(Edis>lengthThreshold){
                                continue;
                            } else{
                                if(Edis>lastDis){
                                    continue;
                                } else{
                                    POIDis= distanceQuery(node1+1,node2+1);
                                    lastDis= MAX(lastDis,POIDis);
                                    //lastDis=POIDis;
                                }
                            }
                        }
                        if (POIDis < lengthThreshold && POIDis!= 0) {
                            int id1=POIGridMap[node1];
                            int id2=POIGridMap[node2];
                            vector<int> PPath;
                            PPath.push_back(node1);
                            PPath.push_back(node2);
                            isCP[node2]=true;
                            bitset<KIND> pathBit = (NodesBit[node1] & uncover | NodesBit[node2] & uncover);
                            vvNodePath.push_back(PPath);
                            vvPPathBit.push_back(pathBit);
                            vvPPathLength.push_back(POIDis);
                            int dis1=POIDis+ distanceQuery(PPath[0]+1,s+1)+ distanceQuery(PPath[PPath.size()-1]+1,t+1)-OD;
                            int dis2=POIDis+distanceQuery(PPath[0]+1,t+1)+ distanceQuery(PPath[PPath.size()-1]+1,s+1)-OD;
                            //int dis1=POIDis+vSPTDistanceS[PPath[0]]+vSPTDistanceT[PPath[PPath.size()-1]]-OD;
                            //int dis2=POIDis+vSPTDistanceT[PPath[0]]+vSPTDistanceS[PPath[PPath.size()-1]]-OD;
                            int smallDis=(dis1<dis2)?dis1:dis2;
                            qPath.update(vvNodePath.size()-1,smallDis);
                            //qPath.update(vvNodePath.size()-1,dis);
                            count++;
                            pNum++;
                            //cout<<count<<endl;
                        }
                    }
                }
            }
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span2 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    //cout<<"Init time is "<<time_span2.count()<<endl;

    //cout<<pNum<<endl;

//
////   bitset<KIND> pathBit=spBit|vvPPathBit[topPathID];
////   cout<<pathBit.count()<<endl;
    int NUM=0;
    int num=0;
    vector<int> resultPath;

   while(!qPath.empty()){
       vector<int> H2HPA;
       vector<int> H2HPB;
       vector<bitset<KIND>> H2HPABit;
       vector<bitset<KIND>> H2HPBBit;
       //cout<<qPath.size()<<endl;
       NUM++;
       //cout<<NUM<<endl;
       qPath.extract_min(topPathID,topPathDis);
       vector<int> currentPath(vvNodePath[topPathID]);
       bitset<KIND> currentPathBit(vvPPathBit[topPathID]);
       int pathBegin=vvNodePath[topPathID][0];
       int pathEnd=vvNodePath[topPathID].back();
       int H2HD1;
       int H2HD2;
       int H2HD3;
       int H2HD4;
       t1 = std::chrono::high_resolution_clock::now();
       vector<int> H2HSBegin;vector<bitset<KIND>> H2HSBeginBit;
       vector<int> H2HSEnd;vector<bitset<KIND>> H2HSEndBit;
       vector<int> H2HTBegin;vector<bitset<KIND>> H2HTBeginBit;
       vector<int> H2HTEnd;vector<bitset<KIND>> H2HTEndBit;

       H2HD1=H2HPath(pathBegin,s,H2HSBegin,H2HSBeginBit);
       H2HD2=H2HPath(pathEnd,s,H2HSEnd,H2HSEndBit);
       H2HD3=H2HPath(pathBegin,t,H2HTBegin,H2HTBeginBit);
       H2HD4=H2HPath(pathEnd,t,H2HTEnd,H2HTEndBit);

//       std::thread thread1([&]() {  H2HD1=H2HPath(pathBegin,s,H2HSBegin,H2HSBeginBit);
//           H2HD2=H2HPath(pathEnd,s,H2HSEnd,H2HSEndBit); });
//       //std::thread thread2([&]() {   });
//       std::thread thread3([&]() {  H2HD3=H2HPath(pathBegin,t,H2HTBegin,H2HTBeginBit);H2HD4=H2HPath(pathEnd,t,H2HTEnd,H2HTEndBit); });
//       //std::thread thread4([&]() {   });
//       thread1.join();
//       //thread2.join();
//       thread3.join();
//       //thread4.join();
//       t2 = std::chrono::high_resolution_clock::now();
//       time_span2 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
//       cout<<"H2H Time is "<<time_span2.count()<<endl;
       bitset<KIND> spBit1=QueryBit&(H2HSBeginBit.back()|H2HTEndBit.back());
       bitset<KIND> spBit2=QueryBit&(H2HSEndBit.back()|H2HTBeginBit.back());

       if(((spBit1|currentPathBit).count()==QueryBit.count())||((spBit2|currentPathBit).count()==QueryBit.count())){
           cout<<"current path is a feasible"<<endl;
           int  startNode=vvNodePath[topPathID][0];
           int endNode=vvNodePath[topPathID][vvNodePath[topPathID].size()-1];
           int newSPDis1=H2HD1+vvPPathLength[topPathID]+
                         H2HD4;
           int newSPDis2=H2HD2+vvPPathLength[topPathID]+
                         H2HD3;
           resultPath.clear();
           if(newSPDis1<newSPDis2){
               resultPath.push_back(s);
               resultPath.insert(resultPath.end(),vvNodePath[topPathID].begin(),vvNodePath[topPathID].end());
               resultPath.push_back(t);
           } else{
               resultPath.push_back(s);
               std::reverse(currentPath.begin(), currentPath.end());
               resultPath.insert(resultPath.end(),currentPath.begin(),currentPath.end());
               resultPath.push_back(t);
           }
           break;
       }
       bitset<KIND> uncoverKW=currentPathBit^QueryBit;
       int PPathLength1=H2HD1+vvPPathLength[topPathID]+
                        H2HD4;
       int PPathLength2=H2HD2+vvPPathLength[topPathID]+
                        H2HD3;
       if(min(PPathLength1,PPathLength2)>GLength){
           //cout<<"GLength Pruned"<<endl;
           continue;
       }
//       if(PPathLength1<PPathLength2){
//           if((H2HSBeginBit.back()&vvPPathBit[topPathID]).count()==vvPPathBit[topPathID].count()||
//                   (H2HTEndBit.back()&vvPPathBit[topPathID]).count()==vvPPathBit[topPathID].count()){
//               //cout<<"SP Pruned"<<endl;
//               continue;
//           }
//       } else{
//           if((H2HSEndBit.back()&vvPPathBit[topPathID]).count()==vvPPathBit[topPathID].count()||
//              (H2HTBeginBit.back()&vvPPathBit[topPathID]).count()==vvPPathBit[topPathID].count()){
//               //cout<<"SP Pruned"<<endl;
//               continue;
//           }
//       }
       for(auto kw:QueryWord){
           if(uncoverKW.test(kw)){
               for(auto node:keyNodeMap[kw]){
                   num=vvNodePath[topPathID].size();
                   int dis1=distanceQuery(node+1,vvNodePath[topPathID][0]+1)+vvPPathLength[topPathID];
                   int dis2=distanceQuery(vvNodePath[topPathID].back()+1,node+1)+vvPPathLength[topPathID];
                   if(dis1<(lengthThreshold*num)||dis2<(lengthThreshold*num)){
                       int insertStart=0;
                       if(dis1<dis2){
                           insertStart=1;
                       }
                       vector<int> newPPath;
                       bitset<KIND> newPPahBit(vvPPathBit[topPathID]);
                       if(insertStart==1){
                           newPPath.push_back(node);
                           newPPath.insert(newPPath.end(),currentPath.begin(),currentPath.end());
                           vvPPathLength.push_back(dis1);
                       } else{
                           newPPath.insert(newPPath.end(),currentPath.begin(),currentPath.end());
                           newPPath.push_back(node);
                           vvPPathLength.push_back(dis2);
                       }
                       newPPahBit|=NodesBit[node]&uncoverKW;
                       vvNodePath.push_back(newPPath);
                       vvPPathBit.push_back(newPPahBit);
                       int startNode=vvNodePath[vvNodePath.size()-1][0];
                       int endNode=vvNodePath[vvNodePath.size()-1].back();
                       int newSPDis1=distanceQuery(s+1,startNode+1)+vvPPathLength.back()+
                               distanceQuery(endNode+1,t+1);
                       int newSPDis2= distanceQuery(startNode+1,t+1)+vvPPathLength.back()+
                               distanceQuery(endNode+1,t+1);
                       int newSPDis=(newSPDis1<newSPDis2)?newSPDis1:newSPDis2;
                       //cout<<newSPDis<<endl;
                       qPath.update(vvNodePath.size()-1,newSPDis-OD);
                   }
               }
           }
       }
   }

    int Sum=0;
    if(!resultPath.empty()){
       Sum=PruneRepeatedPoiPath(resultPath);
//       for(int i=0;i<resultPath.size()-1;i++){
//           //cout<<resultPath[i]<<" "<<resultPath[i+1]<<endl;
//           Sum+=distanceQuery(resultPath[i] + 1, resultPath[i + 1] + 1);
//       }
   } else{
        cout<<"error"<<endl;
   }
    cout<<"SKORP path is "<<Sum<<" "<<endl;
    return Sum;
    //cout<<topPathID<<" "<<topPathDis<<endl;
}
void Graph::SY_SKORP_V1(int s,int t,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
              int& resultPathLength,vector<int> &resultPath){
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span1(0.0);
    std::chrono::duration<double> time_span2(0.0);
    t1 = std::chrono::high_resolution_clock::now();
    int OD= distanceQuery(s+1,t+1);
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    vector<int> H2Hpath;
    vector<int> H2HEdge;
    vector<bitset<KIND>> H2HPathBit;

    H2HPath(s,t,H2Hpath,H2HPathBit);
    //cout<<"SPT Time is"<<time_span.count()<<endl;
    if((H2HPathBit.back()&QueryBit).count()==QueryBit.count()){
        //cout<<"shortest path is feasible path"<<endl;
       resultPathLength=OD;
       resultPath=H2Hpath;
        return;
    }
    vector<vector<int>> vvNodePath;
    vector<int> vvPPathLength;
    vector<bitset<KIND>> vvPPathBit;
    map<int,vector<int>> keyNodeMap;
    keyNodeMap.clear();
    t1 = std::chrono::high_resolution_clock::now();
    for(int i=0;i<QueryWord.size();i++){
        vector<int> l;
        keyNodeMap.insert({QueryWord[i],l});
    }
    for(auto node:POICandidate){
        for(auto key:keyNodeMap){
            if(NodesBit[node].test(key.first)){
                keyNodeMap[key.first].push_back(node);
            }
        }
    }
    vector<bool> isCP(nodeNum,false);
    benchmark::heap<2, int, int> qPath(nodeNum*50);
    t2 = std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    //cout<<"Map time is "<<time_span1.count()<<endl;
    //benchmark::heap<2,int,int> qPathGrid(nodeNum);
    int count=0;
    t1 = std::chrono::high_resolution_clock::now();
    int topPathID=0;
    int topPathDis=0;
    //qPathGrid.extract_min(topPathID ,topPathDis);
    int pNum=0;
    //H2HPath(s,t,H2Hpath,H2HEdge,H2HPathBit);

    /**
     * Init path by IG-Tree
     */
//    vector<bool> isDo(POICandidate.size(),false);
//    for(int i=0;i<POICandidate.size();i++){
//        if(!isDo[i]){
//            int GTreeNode=gtreepath[POICandidate[i]].back();
//            auto nodePoint = std::find(GTree[GTreeNode].leafnodes.begin(),
//
//                                       GTree[GTreeNode].leafnodes.end(), POICandidate[i]);
//            int minDis=INT_MAX;
//            for(auto bnode:GTree[GTreeNode].union_borders){
//                minDis= min(distanceQuery(bnode+1,POICandidate[i]+1),minDis);
//            }
//
//            if(minDis>lengthThreshold){
//                cout<<"1"<<endl;
//                //get all POI node in current Node;
//            }else{
//                GTreeNode=GTree[GTreeNode].father;
//            }
////            int index= std::distance(GTree[GTreeNode].leafnodes.begin(),nodePoint);
////            int cindex=GTree[GTreeNode].leafnodes.size()*(index-1);
////            vector<int> minds(GTree[GTreeNode].mind.begin()+cindex,GTree[GTreeNode].mind.begin()+cindex+GTree[GTreeNode].union_borders.size());
//
//        }
//    }

    for(int i=0;i<QueryWord.size();i++) {
        for (int j = i + 1; j < QueryWord.size(); j++) {
            if (uncover.test(QueryWord[i]) && uncover.test(QueryWord[j])) {
                for (auto node1: keyNodeMap[QueryWord[i]]) {
                    int c=1;
                    int POIDis=0;
                    int lastDis=0;
                    for (auto node2: keyNodeMap[QueryWord[j]]) {
                        if(c==1){
                            POIDis = distanceQuery(node1 + 1, node2 + 1);
                            lastDis=POIDis;
                            c++;
                        } else{
                            int Edis= (int)EuclideanDistance(node1,node2);
                            Edis=Edis/10;
                            if(Edis>lengthThreshold){
                                continue;
                            } else{
                                if(Edis>lastDis){
                                    continue;
                                } else{
                                    POIDis= distanceQuery(node1+1,node2+1);
                                    lastDis= MAX(lastDis,POIDis);
                                    //lastDis=POIDis;
                                }
                            }
                        }
                        if (POIDis < lengthThreshold && POIDis!= 0) {
                            int id1=POIGridMap[node1];
                            int id2=POIGridMap[node2];
                            vector<int> PPath;
                            PPath.push_back(node1);
                            PPath.push_back(node2);
                            isCP[node2]=true;
                            bitset<KIND> pathBit = (NodesBit[node1] & uncover | NodesBit[node2] & uncover);
                            vvNodePath.push_back(PPath);
                            vvPPathBit.push_back(pathBit);
                            vvPPathLength.push_back(POIDis);
                            int dis1=POIDis+distanceQuery(PPath[0]+1,s+1)+ distanceQuery(PPath[PPath.size()-1]+1,t+1);
                            int dis2=POIDis+distanceQuery(PPath[0]+1,t+1)+ distanceQuery(PPath[PPath.size()-1]+1,s+1);
                            //int dis1=POIDis+vSPTDistanceS[PPath[0]]+vSPTDistanceT[PPath[PPath.size()-1]]-OD;
                            //int dis2=POIDis+vSPTDistanceT[PPath[0]]+vSPTDistanceS[PPath[PPath.size()-1]]-OD;
                            int smallDis=(dis1<dis2)?dis1:dis2;
                            double ita=smallDis/2.0;
                            //method1
                            //
                            qPath.update(vvNodePath.size()-1,smallDis);

                            //method2
                            //qPath.update(vvNodePath.size()-1,ita*1000);

                            //qPath.update(vvNodePath.size()-1,dis);
                            count++;
                            pNum++;
                            //cout<<count<<endl;
                        }
                    }
                }
            }
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span2 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    //cout<<"Init time is "<<time_span2.count()<<endl;

    //cout<<pNum<<endl;

//
////   bitset<KIND> pathBit=spBit|vvPPathBit[topPathID];
////   cout<<pathBit.count()<<endl;
    int NUM=0;
    int num=0;
    while(!qPath.empty()){
        vector<int> H2HPA;
        vector<int> H2HPB;
        vector<bitset<KIND>> H2HPABit;
        vector<bitset<KIND>> H2HPBBit;
        //cout<<qPath.size()<<endl;
        NUM++;
        //cout<<NUM<<endl;
        qPath.extract_min(topPathID,topPathDis);
        vector<int> currentPath(vvNodePath[topPathID]);
        bitset<KIND> currentPathBit(vvPPathBit[topPathID]);
        int pathBegin=vvNodePath[topPathID][0];
        int pathEnd=vvNodePath[topPathID].back();
        int H2HD1;
        int H2HD2;
        int H2HD3;
        int H2HD4;
        int SBBeginIndex,SBEndIndex;
        int SEBeginIndex,SEEndIndex;
        int TBBeginIndex,TBEndIndex;
        int TEBeginIndex,TEEndIndex;
        t1 = std::chrono::high_resolution_clock::now();
        vector<int> H2HSBegin;vector<bitset<KIND>> H2HSBeginBit;
        vector<int> H2HSEnd;vector<bitset<KIND>> H2HSEndBit;
        vector<int> H2HTBegin;vector<bitset<KIND>> H2HTBeginBit;
        vector<int> H2HTEnd;vector<bitset<KIND>> H2HTEndBit;

        H2HD1=H2HPathPlusBit(s,pathBegin,H2HSBegin,H2HSBeginBit,SBBeginIndex,SBEndIndex,currentPathBit.count());
        H2HD2=H2HPathPlusBit(s,pathEnd,H2HSEnd,H2HSEndBit,SEBeginIndex,SEEndIndex,currentPathBit.count());
        H2HD3=H2HPathPlusBit(pathBegin,t,H2HTBegin,H2HTBeginBit,TBBeginIndex,TBEndIndex,currentPathBit.count());
        H2HD4=H2HPathPlusBit(pathEnd,t,H2HTEnd,H2HTEndBit,TEBeginIndex,TEEndIndex,currentPathBit.count());

//       std::thread thread1([&]() {  H2HD1=H2HPath(pathBegin,s,H2HSBegin,H2HSBeginBit);
//           H2HD2=H2HPath(pathEnd,s,H2HSEnd,H2HSEndBit); });
//       //std::thread thread2([&]() {   });
//       std::thread thread3([&]() {  H2HD3=H2HPath(pathBegin,t,H2HTBegin,H2HTBeginBit);H2HD4=H2HPath(pathEnd,t,H2HTEnd,H2HTEndBit); });
//       //std::thread thread4([&]() {   });
//       thread1.join();
//       //thread2.join();
//       thread3.join();
//       //thread4.join();
//       t2 = std::chrono::high_resolution_clock::now();
//       time_span2 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
//       cout<<"H2H Time is "<<time_span2.count()<<endl;
        bitset<KIND> spBit1=QueryBit&(H2HSBeginBit.back()|H2HTEndBit.back());
        bitset<KIND> spBit2=QueryBit&(H2HSEndBit.back()|H2HTBeginBit.back());
//        cout<<"SP1 "<<currentPathBit.count()<<endl;
//        cout<<"SP2 "<<currentPathBit.count()<<endl;








        int PPathLength1=H2HD1+vvPPathLength[topPathID]+
                         H2HD4;
        int PPathLength2=H2HD2+vvPPathLength[topPathID]+
                         H2HD3;
//        if(min(PPathLength1,PPathLength2)>GLength){
//            //cout<<"GLength Pruned"<<endl;
//            continue;
//        }

        if(PPathLength1<PPathLength2){
            //current path is a feasible path
            if(((spBit1|currentPathBit).count()==QueryBit.count())){
                int  startNode=vvNodePath[topPathID][0];
                int endNode=vvNodePath[topPathID][vvNodePath[topPathID].size()-1];
                if(spBit1.count()==currentPathBit.count()){
                    cout<<"sp feasible"<<endl;
                    resultPath.clear();
                    resultPath.push_back(s);
                    resultPath.insert(resultPath.end(),H2HSBegin.begin()+1,H2HSBegin.end()-1);
                    resultPath.push_back(currentPath[0]);
                    resultPath.push_back(currentPath.back());
                    resultPath.insert(resultPath.end(),H2HTEnd.begin()+1,H2HTEnd.end()-1);
                    resultPath.push_back(t);
                } else{
                    cout<<"current path feasible"<<endl;
                    resultPath.clear();
                    resultPath.push_back(s);
                    resultPath.insert(resultPath.end(),H2HSBegin.begin()+1,H2HSBegin.end()-1);
                    resultPath.insert(resultPath.end(),currentPath.begin(),currentPath.end());
                    resultPath.insert(resultPath.end(),H2HTEnd.begin()+1,H2HTEnd.end()-1);
                    resultPath.push_back(t);
                }
                break;
            }
        }
        if(PPathLength2<PPathLength1){
            if(((spBit2|currentPathBit).count()==QueryBit.count())){
                int  startNode=vvNodePath[topPathID][0];
                int endNode=vvNodePath[topPathID][vvNodePath[topPathID].size()-1];
                std::reverse(currentPath.begin(), currentPath.end());
                if(spBit2.count()==currentPathBit.count()){
                    cout<<"sp feasible"<<endl;
                    resultPath.clear();
                    resultPath.push_back(s);
                    resultPath.insert(resultPath.end(),H2HSEnd.begin()+1,H2HSEnd.end()-1);
                    resultPath.push_back(currentPath[0]);
                    resultPath.push_back(currentPath.back());
                    resultPath.insert(resultPath.end(),H2HTBegin.begin()+1,H2HTBegin.end()-1);
                    resultPath.push_back(t);
                } else{
                    cout<<"current path feasible"<<endl;
                    resultPath.clear();
                    resultPath.push_back(s);
                    resultPath.insert(resultPath.end(),H2HSEnd.begin()+1,H2HSEnd.end()-1);
                    resultPath.insert(resultPath.end(),currentPath.begin(),currentPath.end());
                    resultPath.insert(resultPath.end(),H2HTBegin.begin()+1,H2HTBegin.end()-1);
                    resultPath.push_back(t);
                }
                break;
            }
        }
//        int minBegin=INT_MAX;
//        int minBeginNode;
//        int minEnd=INT_MAX;
//        int minEndNode;
//        if(currentPathBit.count()==(QueryBit.count()-1)){
//            cout<<s<<" "<<t<<" "<<"1111111111"<<endl;
//        }
//        if(currentPathBit.count()==(QueryBit.count()-1)){
//            bitset<KIND> uncover(currentPathBit^QueryBit);
//            for(auto kw:QueryWord){
//                if(uncover.test(kw)){
//
//                    for(auto node:RSP[kw]){
//                        int Len1= distanceQuery(s+1,node+1)+ distanceQuery(node+1,currentPath[0]+1);
//                        int Len2= distanceQuery(currentPath.back()+1,node+1)+ distanceQuery(node+1,t+1);
//                        if(Len1<minBegin){
//                            minBegin=Len1;
//                            minBeginNode=node;
//                        }
//                        if(Len2<minEnd){
//                            minEnd=Len2;
//                            minEndNode=node;
//                        }
//                    }
//                }
//            }
//            if(vvPPathLength[topPathID]+ distanceQuery(minBeginNode+1,currentPath[0]+1)<vvPPathLength[topPathID]+
//                                                                                        distanceQuery(currentPath.back()+1,minEndNode+1)){
//                currentPath.insert(currentPath.begin(),minBeginNode);
//                resultPath.push_back(s);
//                resultPath.insert(resultPath.end(),currentPath.begin(),currentPath.end());
//                resultPath.push_back(t);
//                break;
//            } else{
//                currentPath.insert(currentPath.end(),minEndNode);
//                resultPath.push_back(s);
//                resultPath.insert(resultPath.end(),currentPath.begin(),currentPath.end());
//                resultPath.push_back(t);
//                break;
//            }
//        }
        if(min(PPathLength1,PPathLength2)>GLength){
              //cout<<"GLength Pruned"<<endl;
              continue;
        }
//        if(((spBit1|currentPathBit).count()==QueryBit.count())||((spBit2|currentPathBit).count()==QueryBit.count())){
//            cout<<"current path is a feasible"<<endl;
//            int  startNode=vvNodePath[topPathID][0];
//            int endNode=vvNodePath[topPathID][vvNodePath[topPathID].size()-1];
//            int newSPDis1=H2HD1+vvPPathLength[topPathID]+
//                          H2HD4;
//            int newSPDis2=H2HD2+vvPPathLength[topPathID]+
//                          H2HD3;
//            resultPath.clear();
//            if(newSPDis1<newSPDis2){
//                resultPath.push_back(s);
//                resultPath.insert(resultPath.end(),vvNodePath[topPathID].begin(),vvNodePath[topPathID].end());
//                resultPath.push_back(t);
//            } else{
//                resultPath.push_back(s);
//                std::reverse(currentPath.begin(), currentPath.end());
//                resultPath.insert(resultPath.end(),currentPath.begin(),currentPath.end());
//                resultPath.push_back(t);
//            }
//            break;
//        }
        bitset<KIND> uncoverKW=currentPathBit^QueryBit;
//       if(PPathLength1<PPathLength2){
//           if((H2HSBeginBit.back()&vvPPathBit[topPathID]).count()==vvPPathBit[topPathID].count()||
//                   (H2HTEndBit.back()&vvPPathBit[topPathID]).count()==vvPPathBit[topPathID].count()){
//               //cout<<"SP Pruned"<<endl;
//               continue;
//           }
//       } else{
//           if((H2HSEndBit.back()&vvPPathBit[topPathID]).count()==vvPPathBit[topPathID].count()||
//              (H2HTBeginBit.back()&vvPPathBit[topPathID]).count()==vvPPathBit[topPathID].count()){
//               //cout<<"SP Pruned"<<endl;
//               continue;
//           }
//       }
        for(auto kw:QueryWord){
            if(uncoverKW.test(kw)){
                for(auto node:keyNodeMap[kw]){
                    num=vvNodePath[topPathID].size();
                    int dis1=distanceQuery(node+1,vvNodePath[topPathID][0]+1)+vvPPathLength[topPathID];
                    int dis2=distanceQuery(vvNodePath[topPathID].back()+1,node+1)+vvPPathLength[topPathID];
                    if(dis1<(lengthThreshold*num)||dis2<(lengthThreshold*num)){
                        int insertStart=0;
                        if(dis1<dis2){
                            insertStart=1;
                        }
                        vector<int> newPPath;
                        bitset<KIND> newPPahBit(vvPPathBit[topPathID]);
                        if(insertStart==1){
                            newPPath.push_back(node);
                            newPPath.insert(newPPath.end(),currentPath.begin(),currentPath.end());
                            vvPPathLength.push_back(dis1);
                        } else{
                            newPPath.insert(newPPath.end(),currentPath.begin(),currentPath.end());
                            newPPath.push_back(node);
                            vvPPathLength.push_back(dis2);
                        }
                        newPPahBit|=NodesBit[node]&uncoverKW;
                        vvNodePath.push_back(newPPath);
                        vvPPathBit.push_back(newPPahBit);
                        int startNode=vvNodePath[vvNodePath.size()-1][0];
                        int endNode=vvNodePath[vvNodePath.size()-1].back();
                        int newSPDis1=distanceQuery(s+1,startNode+1)+vvPPathLength.back()+
                                      distanceQuery(endNode+1,t+1);
                        int newSPDis2= distanceQuery(startNode+1,t+1)+vvPPathLength.back()+
                                       distanceQuery(endNode+1,t+1);
                        double ita=0;
                        if(newSPDis1<newSPDis2){
                             ita=vvPPathLength.back()/((newPPahBit.count())*1.0);
                        } else{
                             ita=vvPPathLength.back()/((newPPahBit.count())*1.0);
                        }
                        int newSPDis=(newSPDis1<newSPDis2)?newSPDis1:newSPDis2;
                        //cout<<newSPDis<<endl;

                        //method1
                        //qPath.update(vvNodePath.size()-1,ita*1000);

                        //method2
                        qPath.update(vvNodePath.size()-1,newSPDis);
                    }
                }
            }
        }
    }

    int Sum=0;
    if(!resultPath.empty()){
        resultPathLength=PruneRepeatedPoiPath(resultPath);
//       for(int i=0;i<resultPath.size()-1;i++){
//           //cout<<resultPath[i]<<" "<<resultPath[i+1]<<endl;
//           Sum+=distanceQuery(resultPath[i] + 1, resultPath[i + 1] + 1);
//       }
    } else{
        resultPathLength=GLength;
        cout<<"error"<<endl;
    }
    //cout<<topPathID<<" "<<topPathDis<<endl;

}

int Graph::getMaxLB(vector<int> &vpath, bitset<KIND> &vpathBit, map<int,vector<int>> keyNodeMap,int DANode,int ID2) {
    //get keyword that no covered
    bitset<KIND> uncoveredKeyWords(vpathBit^QueryBit);
    int MAX=INT64_MIN;
    int k=0;
    int maxLBNode=0;
    if(uncoveredKeyWords.count()==0)
        return 0;
    for(auto i:keyNodeMap){
        if(uncoveredKeyWords.test(i.first)){
            int MIN=INT_MAX;
            for(auto keyword:i.second){
                int a= distanceQuery(DANode+1,keyword+1);
                int b= distanceQuery(keyword+1,ID2+1);
                if(a+b<=MIN){
                    MIN=a+b;
                }
                if(MIN<MAX&&k!=0){
                    break;
                }
            }
            MAX=max(MAX,MIN);
        }
    }
//    for(auto i:QueryWord){
//        if(uncoveredKeyWords.test(i)){
//            // if test ==true meas that keyword is uncovered
//            int MIN=INT_MAX;
//            for(auto keyword:RSP[i]){
//                //get each POI that has this keyword
//                int a= distanceQuery(DANode+1,keyword+1);
//                int b= distanceQuery(keyword+1,ID2+1);
//                if(a+b<=MIN){
//                    MIN=a+b;
//                }
//                /**
//                 * a method that can prune some POI
//                 */
//                if(MIN<MAX&&k!=0){
//                    break;
//                }
//            }
//            MAX=max(MAX,MIN);
//        }
//        k++;
//    }
    return MAX;
}
int Graph::RangePOI(int s, int t, int dist,int RangeParameter){
    RPOI.clear();
    for(auto kw:QueryWord){
        //vector<int> temp;
        for(auto node:RSP[kw]){
            int ds= distanceQuery(s+1,node+1);
            int dt= distanceQuery(node+1,t+1);
            if((ds+dt)<RangeParameter*dist){
                //temp.push_back(node);
                RPOI.push_back(node);
            }
        }
        while(RPOI.empty()){
            //increase Range
          RangeParameter+=1;
            for(auto node:RSP[kw]){
                int ds= distanceQuery(s+1,node+1);
                int dt= distanceQuery(node+1,t+1);
                if((ds+dt)<RangeParameter*dist){
                    RPOI.push_back(node);
                }
            }
        }
    }
    return RPOI.size();
}
struct Point3D{
    double x,y,z;
    Point3D(double x_,double y_,double z_):x(x_),y(y_),z(z_){}
};
Point3D vectorSubtraction(double x1,double y1,double z1,double x2,double y2,double z2){
    return  Point3D(x1-x2,y1-y2,z1-z2);
}
double dotProduct(Point3D a,Point3D b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}
int Graph::determineInsertPosition(vector<int> path,int node){
    double lat1=Nodes[path[0]].x*M_PI*180.0;
    double lng1=Nodes[path[0]].y*M_PI*180.0;
    double x1=6371.0*cos(lat1)*cos(lng1);
    double y1=6371.0*cos(lat1)*sin(lng1);
    double z1=6371.0*sin(lat1);

    double lat2=Nodes[path.back()].x*M_PI*180.0;
    double lng2=Nodes[path.back()].y*M_PI*180.0;
    double x2=6371.0*cos(lat2)*cos(lng2);
    double y2=6371.0*cos(lat2)*sin(lng2);
    double z2=6371.0*sin(lat2);


    double lat3=Nodes[node].x*M_PI*180.0;
    double lng3=Nodes[node].y*M_PI*180.0;
    double x3=6371.0*cos(lat3)*cos(lng3);
    double y3=6371.0*cos(lat3)*sin(lng3);
    double z3=6371.0*sin(lat3);
    Point3D AB= vectorSubtraction(x2,y2,z2,x1,y1,z1);
    Point3D AP= vectorSubtraction(x3,y3,z3,x1,y1,z1);
    Point3D BP= vectorSubtraction(x3,y3,z3,x2,y2,z2);
    double dotAP= dotProduct(AB,AP);
    double dotBP= dotProduct(AB,BP);
    if(dotAP>0){
        return 1;
    } else if(dotBP<0){
        return 0;
    } else{
        return 1;
    }
}
void Graph::getRangePairPath(map<int,vector<int>> &nodeStartMap, map<int,vector<int>> &nodeEndMap,vector<int> &pairPathPosition,benchmark::heap<2, int, int> &qPath, vector<pair<int,int>> &nodeST,int subgraph,int poi_node,int aveDis,bitset<SIZE1> &pairPathLengthBit,bitset<SIZE1> &subgraphBit,vector<vector<int>> &vvNodePath,vector<int> &vvNodePathLength, vector<bitset<KIND>> &vvNodePathBit) {
    //returnPath.clear();
    //get poi_node's all adj subgraph
    vector<int> newPairPath;
    bitset<KIND> newPairPathBit;
    bitset<SIZE1> border_subgraphBit;
    bitset<KIND> uncoverKW = (QueryBit ^ NodesBit[poi_node]) & QueryBit;
    int pair_dis=0;
    //if current subgraph has uncover kw
    int s_newPath_t_1,s_newPath_t_2;

    //test current subgraph
    if((GTree[subgraph].Gkeyword&uncoverKW).count()>0){
        for(auto poi:GTree[subgraph].POINode){
            if((NodesBit[poi]&uncoverKW).count()>0&&subgraphBit.test(poi)){
                pair_dis= distanceQuery(poi_node+1,poi+1);
                if(pair_dis<aveDis*para) {
                    newPairPath.clear();
                    newPairPath.push_back(poi_node);
                    newPairPath.push_back(poi);
                    //std::sort(newPairPath.begin(), newPairPath.end());
                    newPairPathBit.reset();
                    newPairPathBit |= (NodesBit[poi_node] | NodesBit[poi]) & QueryBit;
                    if((newPairPathBit&NodesBit[poi_node])==newPairPathBit||(newPairPathBit&NodesBit[poi])==newPairPathBit){
                        //keywords dominated
                        continue;
                    }
                    if(!pairPathLengthBit.test(pair_dis)){
                        vvNodePath.push_back(newPairPath);
                        vvNodePathLength.push_back(pair_dis);
                        vvNodePathBit.push_back(newPairPathBit);
                        s_newPath_t_1 = pair_dis + nodeST[poi_node].first + nodeST[poi].second;
                        s_newPath_t_2 = pair_dis + nodeST[poi_node].second + nodeST[poi].first;
                        if(s_newPath_t_1<=s_newPath_t_2){
                            pairPathPosition.push_back(1);// path is s->poi_node->poi->t;
                        } else{
                            pairPathPosition.push_back(0);//path is t->poi->poi_node->s;
                        }
                        nodeStartMap[poi_node].push_back(vvNodePath.size()-1);
                        nodeEndMap[poi].push_back(vvNodePath.size()-1);
                        pairPathLengthBit.set(pair_dis);
                        double ita=0;
                        ita=vvNodePathLength.back()/(newPairPathBit.count()*1.0);
                        //qPath.update(vvNodePath.size() - 1, ita*1000);
                        // A*distance+(1/keywords)*B
                        int newPathKWCount=newPairPathBit.count();
                        double disCount=(min(s_newPath_t_1, s_newPath_t_2)-sp)/(sp*1.0);
                        double kwCount=(QueryBit.count()-newPairPathBit.count())/(QueryBit.count()*1.0);
                        //qPath.update(vvNodePath.size() - 1, min(s_newPath_t_1, s_newPath_t_2));
                        //cout<<(disCount*K1+kwCount*K2)*100<<endl;
                        //qPath.update(vvNodePath.size() - 1,(disCount*K1+kwCount*K2)*10000);
                        qPath.update(vvNodePath.size() - 1, min(s_newPath_t_1, s_newPath_t_2));
                    }
                }
            }
        }
    }
    int node_to_subgraph_dis=0;
    vector<pair<int,vector<int>>> subgraphList;
    for(auto& adj_subgraph:GTree[subgraph].borderIGNodeMap){
        subgraphList.emplace_back(adj_subgraph);
    }
    int levelCount=0;
    while (node_to_subgraph_dis<=aveDis*para){
        //levelCount++;
        int last_min_dis=INT_MAX;
        //cout<<node_to_subgraph_dis<<endl;
        vector<pair<int,vector<int>>> temp_subgraphList;
        //new_subgraphbit.reset();
       for(auto current_sub_graph:subgraphList){
           //current sub graph not enter
           if(!border_subgraphBit.test(current_sub_graph.first)){
               // set current sub graph has enter
               border_subgraphBit.set(current_sub_graph.first);
               // if current graph as uncover kw
               if((GTree[current_sub_graph.first].Gkeyword&uncoverKW).count()>0) {
                   // first is get its first level
                   for (auto poi: GTree[current_sub_graph.first].POINode) {
                       if ((NodesBit[poi] & uncoverKW).count() > 0&&subgraphBit.test(poi)) {
                           pair_dis = distanceQuery(poi_node + 1, poi + 1);
                           if (pair_dis <aveDis*para) {
                               newPairPath.clear();
                               newPairPath.push_back(poi_node);
                               newPairPath.push_back(poi);
                               //std::sort(newPairPath.begin(), newPairPath.end());
                               newPairPathBit.reset();
                               newPairPathBit |= (NodesBit[poi_node] | NodesBit[poi]) & QueryBit;
                               if((newPairPathBit&NodesBit[poi_node])==newPairPathBit||(newPairPathBit&NodesBit[poi])==newPairPathBit){
                                   continue;
                               }
                               if(!pairPathLengthBit.test(pair_dis)){
                                   pairPathLengthBit.set(pair_dis);
                                   vvNodePath.push_back(newPairPath);
                                   vvNodePathLength.push_back(pair_dis);
                                   vvNodePathBit.push_back(newPairPathBit);
                                   s_newPath_t_1 = pair_dis + nodeST[poi_node].first + nodeST[poi].second;
                                   s_newPath_t_2 = pair_dis + nodeST[poi_node].second + nodeST[poi].first;
                                   if(s_newPath_t_1<=s_newPath_t_2){
                                       pairPathPosition.push_back(1);// path is s->poi_node->poi->t;
                                   } else{
                                       pairPathPosition.push_back(0);//path is t->poi->poi_node->s;
                                   }
                                   nodeStartMap[poi_node].push_back(vvNodePath.size()-1);
                                   nodeEndMap[poi].push_back(vvNodePath.size()-1);
                                   double ita=0;
                                   ita=vvNodePathLength.back()/(newPairPathBit.count()*1.0);
                                   //qPath.update(vvNodePath.size() - 1,ita*1000);
                                   int newPathKWCount=newPairPathBit.count();
                                   double disCount=(min(s_newPath_t_1, s_newPath_t_2)-sp)/(sp*1.0);
                                   double kwCount=(QueryBit.count()-newPairPathBit.count())/(QueryBit.count()*1.0);
                                   //qPath.update(vvNodePath.size() - 1,(disCount*K1+kwCount*K2)*10000);
                                   qPath.update(vvNodePath.size() - 1, min(s_newPath_t_1, s_newPath_t_2));
                                   //qPath.update(vvNodePath.size() - 1,(disCount*K1+kwCount*K2)*10000);
                               }
                           }
                       }
                   }
                   for (auto next_level_subgraph: GTree[current_sub_graph.first].borderIGNodeMap) {
                       if (!border_subgraphBit.test(next_level_subgraph.first)) {
                           temp_subgraphList.emplace_back(next_level_subgraph);
                       }
                   }
               }
               else{
                   node_to_subgraph_dis=INT_MAX;
                   for(auto border_node:current_sub_graph.second){
                       if(border_node!=poi_node){
                           int poi_to_border= distanceQuery(poi_node+1,border_node+1);
                           node_to_subgraph_dis= min(node_to_subgraph_dis,poi_to_border);
                       }
                   }
                   if(node_to_subgraph_dis<aveDis*para){
                       // its next level sub_graph maybe have kw
                       for(auto next_level_subgraph:GTree[current_sub_graph.first].borderIGNodeMap){
                           if(!border_subgraphBit.test(next_level_subgraph.first)){
                               temp_subgraphList.emplace_back(next_level_subgraph);
                           }
                       }
                   }
               }
               last_min_dis= min(last_min_dis,node_to_subgraph_dis);
           }
       }
       node_to_subgraph_dis=last_min_dis;
       subgraphList.clear();
       subgraphList=temp_subgraphList;
    }
    //cout<<"extend level "<<levelCount<<endl;
}
void Graph::TWE(int ID1,int ID2,vector<int> QueryList,vector<vector<int>> &kPaths, vector<int> &kPathLength,ofstream &ofst){
    //sp path has cover all kw
    vector<int> H2Hpath;
    vector<bitset<KIND>> H2HPathBit;
    int ODDis=H2HPath(ID1, ID2, H2Hpath, H2HPathBit);
    s=ID1;
    t=ID2;
    vector<vector<int>> CPOI_RSP(KIND);
    if ((H2HPathBit.back() & QueryBit).count() == QueryBit.count()) {
        kPaths.push_back(H2Hpath);
        kPathLength.push_back(ODDis);
//        return;
    }
    CPOI_list.clear();
    vector<int> CPOIByIGTree;
    vector<int> CPOIByTriangle;
    vector<int> CPOIByGreedy;
    map<int,vector<int>> nodeIGNodeMap;
    bitset<SIZE1> subGraphBit;
    vector<vector<int>> vvNodePath;
    vector<int> vvNodePathLength;
    vector<int> vvPathLCAMap;
    vector<int> vvPathLocation;// 1 is left ; 0 is right
    vector<bitset<KIND>> vvNodePathBit;
    int globalUB;
    int aveDis=0;
    getPOICoordinate(ID1,ID2,QueryWord,CPOIByIGTree);
    int LCANode= getLowestCommonAncestor(ID1,ID2);
    TriangleFiltering(ID1,ID2,CPOIByIGTree,CPOIByTriangle);
    //CPOIByTriangle=CPOIByIGTree;
    vector<int> GPath1,GPath2;
    int maxDisA=0;
    int maxDisB=0;
    int GreedyPath1= buildGreedyPath(ID1,ID2,CPOIByTriangle,GPath1);
    int GreedyPath2= buildGreedyPath(ID2,ID1,CPOIByTriangle,GPath2);
    int GreedyPath=min(GreedyPath1,GreedyPath2);
    //set IGpath as UB
    if(GreedyPath1<GreedyPath2){
        globalUB=GreedyPath1- distanceQuery(ID1+1,GPath1[1]+1)- distanceQuery(GPath1[GPath1.size()-2]+1,ID2+1);
        aveDis=globalUB/(GPath1.size()-3);
    } else{
        globalUB=GreedyPath2- distanceQuery(ID2+1,GPath2[1]+1)- distanceQuery(GPath2[GPath2.size()-2]+1,ID1+1);
        aveDis=globalUB/(GPath2.size()-3);
    }
    bitset<KIND> greedyBit;
    //pair to save s -> node -> t dis
    vector<pair<int,int>> nodeST(nodeNum);
    int s_node,node_t,lbDis;
    for(auto node:CPOIByTriangle){
        s_node=distanceQuery(ID1+1,node+1);
        node_t=distanceQuery(node+1,ID2+1);
        lbDis=s_node+node_t;
        if(lbDis<GreedyPath){
            nodeST[node]= make_pair(s_node,node_t);
            CPOIByGreedy.push_back(node);
            CPOI_list.push_back(node);
            greedyBit|=NodesBit[node]&QueryBit;
            nodeIGNodeMap[gtreepath[node].back()].push_back(node);
            subGraphBit.set(node);
            for(auto kw: Nodes[node].kwList){
                CPOI_RSP[kw].push_back(node);
            }
        }
    }
    if(greedyBit.count()!=QueryBit.count()){
        CPOIByGreedy.clear();
        CPOI_list.clear();
        CPOI_list=CPOIByTriangle;
        CPOIByGreedy=CPOIByTriangle;
    }
    vector<int> newPath;
    bitset<KIND> newPathBit;
    int newPathLength;
    int s_newPath_t_1,s_newPath_t_2;
    benchmark::heap<2, int, int> qPath;
    vector<vector<int>> returnPath;
    int c=0;
    bitset<SIZE1> pairPathLengthBit;
    map<int,vector<int>> nodeStartMap;
    map<int,vector<int>> nodeEndMap;
    for(auto IGNode:nodeIGNodeMap){
        //IGNode first is subgraph ID, next is its poi node has query kw
        int subgraphID=IGNode.first;
        vector<int> subgraph_node=IGNode.second;
        for(auto poi_node:subgraph_node){
            getRangePairPath(nodeStartMap,nodeEndMap,vvPathLocation,qPath,nodeST,subgraphID,poi_node,aveDis,pairPathLengthBit,subGraphBit,vvNodePath,vvNodePathLength,vvNodePathBit);
        }
    }
    //cout<<"Init path size is "<<vvNodePath.size()<<endl;
    int topPathID;
    int topPathAdis;
    //generate init two POI path
    vector<int> temp_path;
    bitset<KIND> temp_pathBit;
    int temp_path_Length;
    int temp_path_location;
    int poi_dis;
    initPairSum+=vvNodePath.size();
    vector<int> H2HPathS;
    vector<int> H2HPathT;
    int HS,HT;
    bitset<KIND> H2HPathSBit;
    bitset<KIND> H2HPathTBit;
    int resultPathLength;
    bitset<KIND> uncoverKW;
    int edgeNum=0;
    int start_poi;
    int end_poi;
    int count=0;
    int flag=1;
    while (!qPath.empty()){
        //cout<<++count<<endl;
        qPath.extract_min(topPathID,topPathAdis);
        temp_path=vvNodePath[topPathID];
        temp_pathBit=vvNodePathBit[topPathID];
        //cout<<temp_pathBit.count()<<endl;
        temp_path_location=vvPathLocation[topPathID];
        temp_path_Length=vvNodePathLength[topPathID];
        uncoverKW = (QueryBit ^ temp_pathBit) & QueryBit;
        // 1: dequeue path --------------------------------------------------
        if((temp_pathBit&QueryBit).count()==QueryBit.count()){
            //current path is a feasible path
            if(temp_path_location==1){
                resultPathLength=nodeST[temp_path[0]].first+nodeST[temp_path.back()].second+vvNodePathLength[topPathID];
            } else{
                resultPathLength=nodeST[temp_path[0]].second+nodeST[temp_path.back()].first+vvNodePathLength[topPathID];
            }
            //cout<<"total success"<<endl;
            newPath.push_back(ID1);
            newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
            newPath.push_back(ID2);
            kPaths.push_back(newPath);
            kPathLength.push_back(resultPathLength);
            break;
        } else{
            // if s->current->t path also can cover all kw
            if(temp_path_location==1){
                flag=1;
                HS=H2HPath(ID1,temp_path[0],H2HPathS,H2HPathSBit);
                HT=H2HPath(temp_path.back(),ID2,H2HPathT,H2HPathTBit);
            } else{
                flag=0;
                // if t->current->s path also can cover all kw
                HS=H2HPath(ID1,temp_path.back(),H2HPathS,H2HPathSBit);
                HT=H2HPath(ID2,temp_path[0],H2HPathT,H2HPathTBit);
            }
            if(((H2HPathSBit|H2HPathTBit|temp_pathBit)&QueryBit).count()==QueryBit.count()){
                newPath.clear();
                //cout<<"shortest path success "<<topPathID<<endl;
                resultPathLength=HS+HT+vvNodePathLength[topPathID];
                kPathLength.push_back(resultPathLength);
                break;
            }
        }
        //---------------------------------------------------

        // 2: extend current temp_path-----------------------
        // 2.1 find other pair path in vvNodePath
         //first is to find pair_path that end node is temp path start node
         start_poi=temp_path[0];
         end_poi=temp_path.back();
         if(!nodeEndMap[start_poi].empty()){
             vector<int> pair_pathIDList = nodeEndMap[start_poi];
             for(auto pair_pathID:pair_pathIDList){
                if((vvNodePathBit[pair_pathID]&uncoverKW).count()>0){
                    newPath.clear();
                    newPath.insert(newPath.end(),vvNodePath[pair_pathID].begin(),vvNodePath[pair_pathID].end());
                    newPath.insert(newPath.end(),temp_path.begin()+1,temp_path.end());
                    newPathBit=(vvNodePathBit[pair_pathID]|temp_pathBit)&QueryBit;
                    newPathLength=vvNodePathLength[pair_pathID]+temp_path_Length;
                    vvNodePath.push_back(newPath);
                    vvNodePathLength.push_back(newPathLength);
                    vvNodePathBit.push_back(newPathBit);
                    s_newPath_t_1 = newPathLength + nodeST[newPath[0]].first + nodeST[newPath.back()].second;
                    s_newPath_t_2 = newPathLength + nodeST[newPath.back()].first + nodeST[newPath[0]].second;
                    if(s_newPath_t_1<=s_newPath_t_2){
                        vvPathLocation.push_back(1);// path is s->poi_node->poi->t;
                    } else{
                        vvPathLocation.push_back(0);//path is s->poi->poi_node->t;
                    }
                    nodeStartMap[newPath[0]].push_back(vvNodePath.size()-1);
                    nodeEndMap[newPath.back()].push_back(vvNodePath.size()-1);
                    double ita=0;
                    ita=vvNodePathLength.back()/(newPathBit.count()*1.0);
                    //int maxLB=getMaxLB(newPairPath,pairPathPosition.back(),CPOI_list,newPairPathBit,s,t);
                    //qPath.update(vvNodePath.size() - 1,ita*1000);
                    qPath.update(vvNodePath.size()-1, min(s_newPath_t_1,s_newPath_t_2));
                }
             }
         }else if(!nodeStartMap[end_poi].empty()){
             vector<int> pair_pathIDList = nodeStartMap[end_poi];
             for(auto pair_pathID:pair_pathIDList){
                 if((vvNodePathBit[pair_pathID]&uncoverKW).count()>0){
                     newPath.clear();
                     newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                     newPath.insert(newPath.end(),vvNodePath[pair_pathID].begin()+1,vvNodePath[pair_pathID].end());
                     newPathBit=(vvNodePathBit[pair_pathID]|temp_pathBit)&QueryBit;
                     newPathLength=vvNodePathLength[pair_pathID]+temp_path_Length;
                     vvNodePath.push_back(newPath);
                     vvNodePathLength.push_back(newPathLength);
                     vvNodePathBit.push_back(newPathBit);
                     s_newPath_t_1 = newPathLength + nodeST[newPath[0]].first + nodeST[newPath.back()].second;
                     s_newPath_t_2 = newPathLength + nodeST[newPath.back()].first + nodeST[newPath[0]].second;
                     if(s_newPath_t_1<=s_newPath_t_2){
                         vvPathLocation.push_back(1);// path is s->poi_node->poi->t;
                     } else {
                         vvPathLocation.push_back(0);//path is s->poi->poi_node->t;
                     }
                     nodeStartMap[newPath[0]].push_back(vvNodePath.size()-1);
                     nodeEndMap[newPath.back()].push_back(vvNodePath.size()-1);
                     double ita=0;
                     ita=vvNodePathLength.back()/(newPathBit.count()*1.0);
                     //int maxLB=getMaxLB(newPairPath,pairPathPosition.back(),CPOI_list,newPairPathBit,s,t);
                     //qPath.update(vvNodePath.size() - 1,ita*1000);
                     qPath.update(vvNodePath.size()-1, min(s_newPath_t_1,s_newPath_t_2));
                     //qPath.update(vvNodePath.size() - 1, newPathLength+maxLB);
                 }
             }
         }
         else{
             bool front=false;
             bool end=false;
             int front_dis;
             int end_dis;
             for(auto kw:QueryWord){
                 if(uncoverKW.test(kw)){
                     for(auto poi:CPOI_RSP[kw]){
                         front_dis= distanceQuery(poi+1,temp_path[0]+1);
                         end_dis= distanceQuery(temp_path.back()+1,poi+1);
                         if(front_dis+temp_path_Length<=(temp_path.size()*aveDis*para)&&front_dis>aveDis*para){
                             front= true;
                         }
                         if(end_dis+temp_path_Length<=(temp_path.size()*aveDis*para)&&front_dis>aveDis*para){
                             end= true;
                         }
                         if(!front&&!end){
                             continue;
                         }
                         if(front&&end){
                             if(front_dis<end_dis){
                                 newPath.clear();
                                 newPath.push_back(poi);
                                 newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                                 newPathLength=temp_path_Length+front_dis;
                             } else{
                                 newPath.clear();
                                 newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                                 newPath.push_back(poi);
                                 newPathLength=temp_path_Length+end_dis;
                             }
                         } else if(front && !end){
                             newPath.clear();
                             newPath.push_back(poi);
                             newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                             newPathLength=temp_path_Length+front_dis;
                         }else if(!front&& end){
                             newPath.clear();
                             newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                             newPath.push_back(poi);
                             newPathLength=temp_path_Length+end_dis;
                         }
                         newPathBit=(NodesBit[poi]|temp_pathBit)&QueryBit;
                         vvNodePath.push_back(newPath);
                         vvNodePathLength.push_back(newPathLength);
                         vvNodePathBit.push_back(newPathBit);
                         s_newPath_t_1 = newPathLength + nodeST[newPath[0]].first + nodeST[newPath.back()].second;
                         s_newPath_t_2 = newPathLength + nodeST[newPath.back()].first + nodeST[newPath[0]].second;
                         if(s_newPath_t_1<=s_newPath_t_2){
                             vvPathLocation.push_back(1);// path is s->poi_node->poi->t;
                         } else {
                             vvPathLocation.push_back(0);//path is s->poi->poi_node->t;
                         }
                         nodeStartMap[newPath[0]].push_back(vvNodePath.size()-1);
                         nodeEndMap[newPath.back()].push_back(vvNodePath.size()-1);
                         double ita=0;
                         ita=vvNodePathLength.back()/(newPathBit.count()*1.0);
                         //int maxLB=getMaxLB(newPairPath,pairPathPosition.back(),CPOI_list,newPairPathBit,s,t);
                         //qPath.update(vvNodePath.size() - 1,ita*1000);
                         qPath.update(vvNodePath.size()-1, min(s_newPath_t_1,s_newPath_t_2));
                         //int maxLB=getMaxLB(newPath,vvPathLocation.back(),CPOI_list,newPathBit,s,t);
                     }
                 }
             }
         }
    }
    if(kPathLength.empty()){
        cout<<"error"<<endl;
        kPathLength.push_back(GreedyPath);
    }
}
void Graph::KTWE(int ID1,int ID2,int K,vector<int> QueryList,vector<vector<int>> &kPaths, vector<int> &kPathLength,ofstream &ofst){

    //collectLeafByLevel();
    //sp path has cover all kw
    vector<int> H2Hpath;
    vector<bitset<KIND>> H2HPathBit;
    int ODDis=H2HPath(ID1, ID2, H2Hpath, H2HPathBit);
    s=ID1;
    t=ID2;
    vector<vector<int>> CPOI_RSP(KIND);
    if ((H2HPathBit.back() & QueryBit).count() == QueryBit.count()) {
        kPaths.push_back(H2Hpath);
        kPathLength.push_back(ODDis);
        return;
    }
    CPOI_list.clear();
    vector<int> CPOIByIGTree;
    vector<int> CPOIByTriangle;
    vector<int> CPOIByGreedy;
    map<int,vector<int>> nodeIGNodeMap;
    bitset<SIZE1> subGraphBit;
    vector<vector<int>> vvNodePath;
    vector<int> vvNodePathLength;
    vector<int> vvPathLCAMap;
    vector<int> vvPathLocation;// 1 is left ; 0 is right
    vector<bitset<KIND>> vvNodePathBit;
    int globalUB;
    int aveDis=0;
    getPOICoordinate(ID1,ID2,QueryWord,CPOIByIGTree);
    vector<int> GPath1,GPath2;
    TriangleFiltering(ID1,ID2,CPOIByIGTree,CPOIByTriangle);
    //CPOIByTriangle=CPOIByIGTree;
    int maxDisA=0;
    int maxDisB=0;
    //CPOIByIGTree=POIObjects;
    int GreedyPath1= buildGreedyPath(ID1,ID2,CPOIByTriangle,GPath1);
    int GreedyPath2= buildGreedyPath(ID2,ID1,CPOIByTriangle,GPath2);
    int GreedyPath=min(GreedyPath1,GreedyPath2);
    //kPathLength.push_back(GreedyPath);
    //set IGpath as UB
    if(GreedyPath1<GreedyPath2){
        globalUB=GreedyPath1- distanceQuery(ID1+1,GPath1[1]+1)- distanceQuery(GPath1[GPath1.size()-2]+1,ID2+1);
        aveDis=globalUB/(GPath1.size()-3);
    } else{
        globalUB=GreedyPath2- distanceQuery(ID2+1,GPath2[1]+1)- distanceQuery(GPath2[GPath2.size()-2]+1,ID1+1);
        aveDis=globalUB/(GPath2.size()-3);
    }

    bitset<KIND> greedyBit;
    //pair to save s -> node -> t dis
    vector<pair<int,int>> nodeST(nodeNum);
    int s_node,node_t,lbDis;
    for(auto node:CPOIByTriangle){
        s_node=distanceQuery(ID1+1,node+1);
        node_t=distanceQuery(node+1,ID2+1);
        if((s_node+node_t)<GreedyPath){
            nodeST[node]= make_pair(s_node,node_t);
            CPOIByGreedy.push_back(node);
            CPOI_list.push_back(node);
            greedyBit|=NodesBit[node]&QueryBit;
            nodeIGNodeMap[gtreepath[node].back()].push_back(node);
            subGraphBit.set(node);
            for(auto kw: Nodes[node].kwList){
                CPOI_RSP[kw].push_back(node);
            }
        }
    }
    if(greedyBit.count()!=QueryBit.count()){
        CPOIByGreedy.clear();
        CPOI_list.clear();
        CPOI_list=CPOIByTriangle;
        CPOIByGreedy=CPOIByTriangle;
        //cout<<"-------------------------------"<<endl;
    }
    vector<int> newPath;
    bitset<KIND> newPathBit;
    int newPathLength;
    int s_newPath_t_1,s_newPath_t_2;
    benchmark::heap<2, int, int> qPath;
    vector<vector<int>> returnPath;
    int c=0;
    bitset<SIZE1> pairPathLengthBit;
    map<int,vector<int>> nodeStartMap;
    map<int,vector<int>> nodeEndMap;
    for(auto IGNode:nodeIGNodeMap){
        //IGNode first is subgraph ID, next is its poi node has query kw
        int subgraphID=IGNode.first;
        vector<int> subgraph_node=IGNode.second;
        for(auto poi_node:subgraph_node){
            getRangePairPath(nodeStartMap,nodeEndMap,vvPathLocation,qPath,nodeST,subgraphID,poi_node,aveDis,pairPathLengthBit,subGraphBit,vvNodePath,vvNodePathLength,vvNodePathBit);
        }
    }
    //cout<<"Init path size is "<<vvNodePath.size()<<endl;
    int topPathID;
    int topPathAdis;
    //generate init two POI path
    vector<int> temp_path;
    bitset<KIND> temp_pathBit;
    int temp_path_Length;
    int temp_path_location;
    int poi_dis;
    initPairSum+=vvNodePath.size();
    vector<int> H2HPathS;
    vector<int> H2HPathT;
    int HS,HT;
    bitset<KIND> H2HPathSBit;
    bitset<KIND> H2HPathTBit;
    int resultPathLength;
    bitset<KIND> uncoverKW;
    int edgeNum=0;
    int start_poi;
    int end_poi;
    int count=0;
    int flag=1;
    int maxPathLegnth=0;
    while (!qPath.empty()){
        //cout<<++count<<endl;
        qPath.extract_min(topPathID,topPathAdis);
        temp_path=vvNodePath[topPathID];
        temp_pathBit=vvNodePathBit[topPathID];
        //cout<<temp_pathBit.count()<<endl;
        temp_path_location=vvPathLocation[topPathID];
        temp_path_Length=vvNodePathLength[topPathID];
        uncoverKW = (QueryBit ^ temp_pathBit) & QueryBit;
        // 1: dequeue path --------------------------------------------------
        if((temp_pathBit&QueryBit).count()==QueryBit.count()){
            //current path is a feasible path
            if(temp_path_location==1){
                resultPathLength=nodeST[temp_path[0]].first+nodeST[temp_path.back()].second+vvNodePathLength[topPathID];
            } else{
                resultPathLength=nodeST[temp_path[0]].second+nodeST[temp_path.back()].first+vvNodePathLength[topPathID];
            }
            if(!kPathLength.empty()&&resultPathLength==kPathLength.back()){
                continue;
            }
            //cout<<"total success"<<endl;
            newPath.clear();
            newPath.push_back(ID1);
            newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
            newPath.push_back(ID2);
            kPathLength.push_back(resultPathLength);
            kPaths.push_back(newPath);
            //kPathLength.push_back(resultPathLength);
            if(resultPathLength<maxPathLegnth){
                resultPathLength=maxPathLegnth;
            }
            if(kPathLength.size()==K){
                break;
            }
        } else{
            // if s->current->t path also can cover all kw
            if(temp_path_location==1){
                flag=1;
                HS=H2HPath(ID1,temp_path[0],H2HPathS,H2HPathSBit);
                HT=H2HPath(temp_path.back(),ID2,H2HPathT,H2HPathTBit);
            } else{
                flag=0;
                // if t->current->s path also can cover all kw
                HS=H2HPath(ID1,temp_path.back(),H2HPathS,H2HPathSBit);
                HT=H2HPath(ID2,temp_path[0],H2HPathT,H2HPathTBit);
            }
            if(((H2HPathSBit|H2HPathTBit|temp_pathBit)&QueryBit).count()==QueryBit.count()){
                newPath.clear();
                resultPathLength=HS+HT+vvNodePathLength[topPathID];
                if(!kPathLength.empty()&&resultPathLength==kPathLength.back()){
                    continue;
                }
                newPath.clear();
                kPathLength.push_back(resultPathLength);
                if(resultPathLength<maxPathLegnth){
                    resultPathLength=maxPathLegnth;
                }
                if(kPathLength.size()==K){
                    break;
                }
            }
        }
        //---------------------------------------------------

        // 2: extend current temp_path-----------------------
        // 2.1 find other pair path in vvNodePath
        //first is to find pair_path that end node is temp path start node
        start_poi=temp_path[0];
        end_poi=temp_path.back();
        if(!nodeEndMap[start_poi].empty()){
            vector<int> pair_pathIDList = nodeEndMap[start_poi];
            for(auto pair_pathID:pair_pathIDList){
                if((vvNodePathBit[pair_pathID]&uncoverKW).count()>0){
                    newPath.clear();
                    newPath.insert(newPath.end(),vvNodePath[pair_pathID].begin(),vvNodePath[pair_pathID].end());
                    newPath.insert(newPath.end(),temp_path.begin()+1,temp_path.end());
                    newPathBit=(vvNodePathBit[pair_pathID]|temp_pathBit)&QueryBit;
                    newPathLength=vvNodePathLength[pair_pathID]+temp_path_Length;
                    s_newPath_t_1 = newPathLength + nodeST[newPath[0]].first + nodeST[newPath.back()].second;
                    s_newPath_t_2 = newPathLength + nodeST[newPath.back()].first + nodeST[newPath[0]].second;
                    if(min(s_newPath_t_1,s_newPath_t_2)>maxPathLegnth&&newPathBit.count()<QueryBit.count()){
                        //cout<<"UB prune"<<endl;
                        continue;
                    }
                    vvNodePath.push_back(newPath);
                    vvNodePathLength.push_back(newPathLength);
                    vvNodePathBit.push_back(newPathBit);
                    if(s_newPath_t_1<=s_newPath_t_2){
                        vvPathLocation.push_back(1);// path is s->poi_node->poi->t;
                    } else{
                        vvPathLocation.push_back(0);//path is s->poi->poi_node->t;
                    }
                    nodeStartMap[newPath[0]].push_back(vvNodePath.size()-1);
                    nodeEndMap[newPath.back()].push_back(vvNodePath.size()-1);
                    double ita=0;
                    ita=vvNodePathLength.back()/(newPathBit.count()*1.0);
                    //int maxLB=getMaxLB(newPairPath,pairPathPosition.back(),CPOI_list,newPairPathBit,s,t);
                    //qPath.update(vvNodePath.size() - 1,ita*1000);
                    //min(s_newPath_t_1, s_newPath_t_2)*A+(1.0/(vvNodePathBit.back().count()))*B

                    double disCount=(min(s_newPath_t_1, s_newPath_t_2)-sp)/(sp*1.0);
                    double kwCount=(QueryBit.count()-newPathBit.count())/(QueryBit.count()*1.0);
                    qPath.update(vvNodePath.size() - 1, min(s_newPath_t_1, s_newPath_t_2));
                    //qPath.update(vvNodePath.size() - 1,(disCount*K1+kwCount*K2)*10000);
                }
            }
        }else if(!nodeStartMap[end_poi].empty()){
            vector<int> pair_pathIDList = nodeStartMap[end_poi];
            for(auto pair_pathID:pair_pathIDList){
                if((vvNodePathBit[pair_pathID]&uncoverKW).count()>0){
                    newPath.clear();
                    newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                    newPath.insert(newPath.end(),vvNodePath[pair_pathID].begin()+1,vvNodePath[pair_pathID].end());
                    newPathBit=(vvNodePathBit[pair_pathID]|temp_pathBit)&QueryBit;
                    newPathLength=vvNodePathLength[pair_pathID]+temp_path_Length;
                    s_newPath_t_1 = newPathLength + nodeST[newPath[0]].first + nodeST[newPath.back()].second;
                    s_newPath_t_2 = newPathLength + nodeST[newPath.back()].first + nodeST[newPath[0]].second;
                    if(maxPathLegnth!=0&&min(s_newPath_t_1,s_newPath_t_2)>maxPathLegnth){
                        continue;
                    }
                    vvNodePath.push_back(newPath);
                    vvNodePathLength.push_back(newPathLength);
                    vvNodePathBit.push_back(newPathBit);
                    if(s_newPath_t_1<=s_newPath_t_2){
                        vvPathLocation.push_back(1);// path is s->poi_node->poi->t;
                    } else {
                        vvPathLocation.push_back(0);//path is s->poi->poi_node->t;
                    }
                    nodeStartMap[newPath[0]].push_back(vvNodePath.size()-1);
                    nodeEndMap[newPath.back()].push_back(vvNodePath.size()-1);
                    double ita=0;
                    ita=vvNodePathLength.back()/(newPathBit.count()*1.0);
                    //int maxLB=getMaxLB(newPairPath,pairPathPosition.back(),CPOI_list,newPairPathBit,s,t);
                    //qPath.update(vvNodePath.size() - 1,ita*1000);
                    double disCount=(min(s_newPath_t_1, s_newPath_t_2)-sp)/(sp*1.0);
                    double kwCount=(QueryBit.count()-newPathBit.count())/(QueryBit.count()*1.0);
                    qPath.update(vvNodePath.size() - 1, min(s_newPath_t_1, s_newPath_t_2));
                    //qPath.update(vvNodePath.size() - 1,(disCount*K1+kwCount*K2)*10000);
                }
            }
        }
        else{
            bool front=false;
            bool end=false;
            int front_dis;
            int end_dis;
            for(auto kw:QueryWord){
                if(uncoverKW.test(kw)){
                    for(auto poi:CPOI_RSP[kw]){
                        front_dis= distanceQuery(poi+1,temp_path[0]+1);
                        end_dis= distanceQuery(temp_path.back()+1,poi+1);
                        if(front_dis+temp_path_Length<=(temp_path.size()*aveDis*para)&&front_dis>aveDis*para){
                            front= true;
                        }
                        if(end_dis+temp_path_Length<=(temp_path.size()*aveDis*para)&&front_dis>aveDis*para){
                            end= true;
                        }
                        if(!front&&!end){
                            continue;
                        }
                        if(front&&end){
                            if(front_dis<end_dis){
                                newPath.clear();
                                newPath.push_back(poi);
                                newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                                newPathLength=temp_path_Length+front_dis;
                            } else{
                                newPath.clear();
                                newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                                newPath.push_back(poi);
                                newPathLength=temp_path_Length+end_dis;
                            }
                        } else if(front && !end){
                            newPath.clear();
                            newPath.push_back(poi);
                            newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                            newPathLength=temp_path_Length+front_dis;
                        }else if(!front&& end){
                            newPath.clear();
                            newPath.insert(newPath.end(),temp_path.begin(),temp_path.end());
                            newPath.push_back(poi);
                            newPathLength=temp_path_Length+end_dis;
                        }
                        newPathBit=(NodesBit[poi]|temp_pathBit)&QueryBit;
                        vvNodePath.push_back(newPath);
                        vvNodePathLength.push_back(newPathLength);
                        vvNodePathBit.push_back(newPathBit);
                        s_newPath_t_1 = newPathLength + nodeST[newPath[0]].first + nodeST[newPath.back()].second;
                        s_newPath_t_2 = newPathLength + nodeST[newPath.back()].first + nodeST[newPath[0]].second;
                        if(maxPathLegnth!=0&&min(s_newPath_t_1,s_newPath_t_2)>maxPathLegnth&&newPathBit.count()<QueryBit.count()){
                            continue;
                        }
                        if(s_newPath_t_1<=s_newPath_t_2){
                            vvPathLocation.push_back(1);// path is s->poi_node->poi->t;
                        } else {
                            vvPathLocation.push_back(0);//path is s->poi->poi_node->t;
                        }
                        nodeStartMap[newPath[0]].push_back(vvNodePath.size()-1);
                        nodeEndMap[newPath.back()].push_back(vvNodePath.size()-1);
                        double ita=0;
                        ita=vvNodePathLength.back()/(newPathBit.count()*1.0);
                        //int maxLB=getMaxLB(newPairPath,pairPathPosition.back(),CPOI_list,newPairPathBit,s,t);
                        //qPath.update(vvNodePath.size() - 1,ita*1000);
                        double disCount=(min(s_newPath_t_1, s_newPath_t_2)-sp)/(sp*1.0);
                        double kwCount=(QueryBit.count()-newPathBit.count())/(QueryBit.count()*1.0);
                        qPath.update(vvNodePath.size() - 1, min(s_newPath_t_1, s_newPath_t_2));
                        //qPath.update(vvNodePath.size() - 1,(disCount*K1+kwCount*K2)*10000);
                    }
                }
            }
        }
    }
    if(kPathLength.empty()){
        cout<<"error"<<endl;
        kPathLength.push_back(GreedyPath);
    }


}
int Graph::getMaxLB(vector<int> currentPath,int location,vector<int> IG_poiList,bitset<KIND> current_pathBit,int ID1,int ID2){
    bitset<KIND> uncoverKW;
    int dis1,dis2;
    int maxLB=0;
        //ID1 currentPath[0]
        uncoverKW = (QueryBit ^ current_pathBit) & QueryBit;
        for(auto kw:QueryWord) {
            int min_dis1 = INT_MAX;
            int min_dis2 = INT_MAX;
            int min_dis = INT_MAX;
            if (uncoverKW.test(kw)) {
                for (auto poi: IG_poiList) {
                    if (NodesBit[poi].test(kw)) {
                        if(location==1){
                            dis1 = distanceQuery(ID1 + 1, poi + 1) + distanceQuery(poi + 1, currentPath[0] + 1);
                            dis2 = distanceQuery(ID2 + 1, poi + 1) + distanceQuery(poi + 1, currentPath.back() + 1);
                        }else {
                            dis1 = distanceQuery(ID2 + 1, poi + 1) + distanceQuery(poi + 1, currentPath[0] + 1);
                            dis2 = distanceQuery(ID1 + 1, poi + 1) + distanceQuery(poi + 1, currentPath.back() + 1);
                        }
                    }
                    min_dis1 = min(min_dis1, dis1);
                    min_dis2 = min(min_dis2, dis2);
                    min_dis = min(min_dis,min_dis1+min_dis2);
                }
                maxLB = max(maxLB, min_dis);
            }
        }
    return maxLB;
}
void Graph::TWE_Base(int ID1, int ID2, vector<int> QueryList, vector<vector<int>> &kPaths, vector<int> &kPathLength,ofstream& ofst){
    //sp path has cover all kw
    vector<int> H2Hpath;
    vector<bitset<KIND>> H2HPathBit;
    int ODDis=H2HPath(ID1, ID2, H2Hpath, H2HPathBit);
    if ((H2HPathBit.back() & QueryBit).count() == QueryBit.count()) {
        kPaths.push_back(H2Hpath);
        kPathLength.push_back(ODDis);
        return;
    }
    vector<pair<int,int>> nodeST(nodeNum);
    vector<int> CPOIByIGTree;
    vector<int> CPOIByTriangle;
    vector<int> CPOIByGreedy;
    map<int,vector<int>> nodeIGNodeMap;
    bitset<SIZE1> subGraphBit;
    vector<vector<int>> vvNodePath;
    vector<int> vvNodePathLength;
    vector<int> vvPathLCAMap;
    vector<int> vvPathLocation;// 0 is left ; 1 is right
    vector<bitset<KIND>> vvNodePathBit;
    int globalUB;
    int aveDis=0;
    getPOICoordinate(ID1,ID2,QueryWord,CPOIByIGTree);
    int LCANode= getLowestCommonAncestor(ID1,ID2);
    TriangleFiltering(ID1,ID2,CPOIByIGTree,CPOIByTriangle);
    vector<int> GPath1,GPath2;
    int maxDisA=0;
    int maxDisB=0;
    int GreedyPath1= buildGreedyPath(ID1,ID2,CPOIByTriangle,GPath1);
    int GreedyPath2= buildGreedyPath(ID2,ID1,CPOIByTriangle,GPath2);
    int GreedyPath=min(GreedyPath1,GreedyPath2);
    //set IGpath as UB
    if(GreedyPath1<GreedyPath2){
        globalUB=GreedyPath1- distanceQuery(ID1+1,GPath1[1]+1)- distanceQuery(GPath1[GPath1.size()-2]+1,ID2+1);
        aveDis=globalUB/(GPath1.size()-3);
    } else{
        globalUB=GreedyPath2- distanceQuery(ID2+1,GPath2[1]+1)- distanceQuery(GPath2[GPath2.size()-2]+1,ID1+1);
        aveDis=globalUB/(GPath2.size()-3);
    }
    bitset<KIND> greedyBit;
    int s_node,node_t,lbDis;
    for(auto node:CPOIByTriangle){
        s_node=distanceQuery(ID1+1,node+1);
        node_t=distanceQuery(node+1,ID2+1);
        lbDis=s_node+node_t;
        if(lbDis<GreedyPath){
            nodeST[node]= make_pair(s_node,node_t);
            CPOIByGreedy.push_back(node);
            greedyBit|=NodesBit[node]&QueryBit;
            nodeIGNodeMap[gtreepath[node].back()].push_back(node);
            subGraphBit.set(gtreepath[node].back());
        }
    }
    if(greedyBit.count()!=QueryBit.count()){
        CPOIByGreedy.clear();
        CPOIByGreedy=CPOIByTriangle;
    }
    vector<int> newPath;
    bitset<KIND> newPathBit;
    benchmark::heap<2, int, int> qPath(nodeNum * 10);
    vector<vector<int>> returnPath;

    vector<int> newpath;
    for(int i=0;i<CPOIByGreedy.size();i++) {
        bitset<KIND> uncoverKW = (QueryBit ^ NodesBit[CPOIByGreedy[i]]) & QueryBit;
        for (int j = i + 1; j < CPOIByGreedy.size(); j++) {
//            if(CPOIByGreedy[i]==145576&&CPOIByGreedy[j]==145344){
//                cout<<"11111"<<endl;
//            }
            if ((NodesBit[CPOIByGreedy[j]] & uncoverKW).count() > 0) {
                int pair_dis = distanceQuery(CPOIByGreedy[i] + 1, CPOIByGreedy[j] + 1);
                if (pair_dis < aveDis) {
                    newPath.clear();
                    newPath.push_back(CPOIByGreedy[j]);
                    newPath.push_back(CPOIByGreedy[i]);
                    //std::sort(newPath.begin(), newPath.end());
                    newPathBit.reset();
                    newPathBit |= (NodesBit[CPOIByGreedy[i]] | NodesBit[CPOIByGreedy[j]]) & QueryBit;
                    if((newPathBit&NodesBit[CPOIByGreedy[i]])==newPathBit||(newPathBit&NodesBit[CPOIByGreedy[j]])==newPathBit){
                        continue;
                    }
                    vvNodePath.push_back(newPath);
                    vvNodePathLength.push_back(pair_dis);
                    vvNodePathBit.push_back(newPathBit);
                    int s_newPath_t_1 =
                            pair_dis + nodeST[CPOIByGreedy[i]].first + nodeST[CPOIByGreedy[j]].second;
                    int s_newPath_t_2 =
                            pair_dis + nodeST[CPOIByGreedy[i]].second + nodeST[CPOIByGreedy[j]].first;
                    qPath.update(vvNodePathLength.size() - 1, min(s_newPath_t_1, s_newPath_t_2));
                }
            }
        }
    }
    cout<<"Init path size is "<<vvNodePath.size()<<endl;
    initPairSum+=vvNodePath.size();
    int topPathID;
    int topPathAdis;
    qPath.extract_min(topPathID,topPathAdis);
    ofst<<topPathID<<" "<<topPathAdis<<endl;
    cout<<topPathID<<" "<<topPathAdis<<endl;
}







