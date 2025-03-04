//
// Created by xiongxing on 8/3/24.
//
#include "graph.h"
void Graph::Top_k_SKORP_V1(int s,int t,int K,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                           vector<int> &resultPathLength,vector<vector<int>> &resultPath) {
    map<int, vector<int>> nodePathTable;
    //map<int,vector<int>> pathTable;//save path that end is node
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::high_resolution_clock::time_point tt1;
    std::chrono::high_resolution_clock::time_point tt2;
    std::chrono::duration<double> time_span1(0.0);
    std::chrono::duration<double> time_span2(0.0);
    t1 = std::chrono::high_resolution_clock::now();
    tt1 = std::chrono::high_resolution_clock::now();
    int OD = distanceQuery(s + 1, t + 1);
    vector<int> H2Hpath;
    vector<int> H2HEdge;
    vector<bitset<KIND>> H2HPathBit;
    H2HPath(s, t, H2Hpath, H2HPathBit);
    resultPathLength.clear();
    resultPath.clear();
    int feasiblePathLength = GLength;
    int UB = GLength;
    if ((H2HPathBit.back() & QueryBit).count() == QueryBit.count()) {
        resultPath.push_back(H2Hpath);
        resultPathLength.push_back(distanceQuery(s + 1, t + 1));
        return;
    }
    vector<vector<int>> vvNodePath;
    vector<bool> nodeLocation(nodeNum);
    vector<int> pathParent(nodeNum);
    vector<int> vvPPathLength;
    vector<bitset<KIND>> vvPPathBit;
    unordered_map<int, vector<int>> keyNodeMap;
    keyNodeMap.clear();
    for (auto node: POICandidate) {
        for (auto kw: Nodes[node].kwList) {
            if (QueryBit.test(kw)) {
                keyNodeMap[kw].push_back(node);
            }
        }
    }
    vector<bool> isCP(nodeNum, false);
    benchmark::heap<2, int, int> qPath(nodeNum * 10);
    vector<int> vsmallDis(nodeNum);
    int count = 0;
    int topPathID = 0;
    int topPathDis = 0;
    //qPathGrid.extract_min(topPathID ,topPathDis);
    int pNum = 0;
    //H2HPath(s,t,H2Hpath,H2HEdge,H2HPathBit);
    for (int i = 0; i < QueryWord.size(); i++) {
        for (int j = i + 1; j < QueryWord.size(); j++) {
            if (uncover.test(QueryWord[i]) && uncover.test(QueryWord[j])) {
                for (auto node1: keyNodeMap[QueryWord[i]]) {
                    int c = 1;
                    int POIDis = 0;
                    int lastDis = 0;
                    for (auto node2: keyNodeMap[QueryWord[j]]) {
                        if (c == 1) {
                            POIDis = distanceQuery(node1 + 1, node2 + 1);
                            lastDis = POIDis;
                            c++;
                        } else {
                            int Edis = (int) EuclideanDistance(node1, node2);
                            Edis = Edis / 10;
                            if (Edis > lengthThreshold) {
                                continue;
                            } else {
                                if (Edis > lastDis) {
                                    continue;
                                } else {
                                    POIDis = distanceQuery(node1 + 1, node2 + 1);
                                    lastDis = MAX(lastDis, POIDis);
                                    //lastDis=POIDis;
                                }
                            }
                        }
                            if(std::find(vvPPathLength.begin(), vvPPathLength.end(),POIDis)==vvPPathLength.end()){
                        if (POIDis < lengthThreshold && POIDis != 0) {
                            vector<int> PPath;
                            if (distanceQuery(s + 1, node1 + 1) + POIDis + distanceQuery(node2 + 1, t + 1) <
                                distanceQuery(s + 1, node2 + 1) + POIDis + distanceQuery(node1 + 1, t + 1)) {
                                PPath.push_back(node1);
                                PPath.push_back(node2);
                            } else {
                                PPath.push_back(node2);
                                PPath.push_back(node1);
                            }
                            bitset<KIND> pathBit = (NodesBit[node1] & uncover | NodesBit[node2] & uncover);
                            if (std::find(vvPPathLength.begin(), vvPPathLength.end(), POIDis) != vvPPathLength.end()) {
                                continue;
                            }
                            vvNodePath.push_back(PPath);
                            vvPPathBit.push_back(pathBit);
                            vvPPathLength.push_back(POIDis);
                            int dis1 = POIDis + distanceQuery(PPath[0] + 1, s + 1) +
                                       distanceQuery(PPath[PPath.size() - 1] + 1, t + 1);
                            int dis2 = POIDis + distanceQuery(PPath[0] + 1, t + 1) +
                                       distanceQuery(PPath[PPath.size() - 1] + 1, s + 1);
                            int smallDis = (dis1 < dis2) ? dis1 : dis2;
                            double ita = smallDis /pathBit.count()*1.0;
                            qPath.update(vvNodePath.size() - 1, smallDis);
                            //qPath.update(vvNodePath.size() - 1, ita);
                            vsmallDis[vvNodePath.size() - 1] = smallDis;
                            pathParent[vvNodePath.size() - 1] = -1;
                            count++;
                            pNum++;
                        }
                        }
                    }
                }
            }
        }
    }
    int NUM = 0;
    int num = 0;
    //qPath.extract_min(topPathID, topPathDis);
    //qPath.extract_min(topPathID, topPathDis);
    t2 = std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    cout << "Init time is " << time_span1.count() << endl;
    vector<int> disVector;
    int sum=0;
//    vector<vector<bool>> nodeEnd(POICandidate.size(),std::vector<bool>(POICandidate.size(), false));
//    vector<vector<bool>> nodeBegin(POICandidate.size(),std::vector<bool>(POICandidate.size(), false));
    while (!qPath.empty() && resultPath.size() < K) {
        vector<int> H2HPA;
        vector<int> H2HPB;
        vector<bitset<KIND>> H2HPABit;
        vector<bitset<KIND>> H2HPBBit;
        NUM++;
        qPath.extract_min(topPathID, topPathDis);
        vector<int> currentPath(vvNodePath[topPathID]);
        bitset<KIND> currentPathBit(vvPPathBit[topPathID]);
        int pathBegin = vvNodePath[topPathID][0];
        int pathEnd = vvNodePath[topPathID].back();
        int H2HD1;
        int H2HD2;
        int H2HD3;
        int H2HD4;
        int SBBeginIndex, SBEndIndex;
        int SEBeginIndex, SEEndIndex;
        int TBBeginIndex, TBEndIndex;
        int TEBeginIndex, TEEndIndex;
        t1 = std::chrono::high_resolution_clock::now();
        vector<int> H2HSBegin;
        vector<bitset<KIND>> H2HSBeginBit;
        vector<int> H2HSEnd;
        vector<bitset<KIND>> H2HSEndBit;
        vector<int> H2HTBegin;
        vector<bitset<KIND>> H2HTBeginBit;
        vector<int> H2HTEnd;
        vector<bitset<KIND>> H2HTEndBit;
        H2HD1 = H2HPathPlusBit(s, pathBegin, H2HSBegin, H2HSBeginBit, SBBeginIndex, SBEndIndex, currentPathBit.count());
        H2HD2 = H2HPathPlusBit(s, pathEnd, H2HSEnd, H2HSEndBit, SEBeginIndex, SEEndIndex, currentPathBit.count());
        H2HD3 = H2HPathPlusBit(pathBegin, t, H2HTBegin, H2HTBeginBit, TBBeginIndex, TBEndIndex, currentPathBit.count());
        H2HD4 = H2HPathPlusBit(pathEnd, t, H2HTEnd, H2HTEndBit, TEBeginIndex, TEEndIndex, currentPathBit.count());
        bitset<KIND> spBit1 = QueryBit & (H2HSBeginBit.back() | H2HTEndBit.back());
        bitset<KIND> spBit2 = QueryBit & (H2HSEndBit.back() | H2HTBeginBit.back());
        int PPathLength1 = H2HD1 + vvPPathLength[topPathID] +
                           H2HD4;
        int PPathLength2 = H2HD2 + vvPPathLength[topPathID] +
                           H2HD3;
        vector<int> feasiblePath;
        feasiblePath.clear();
        if (PPathLength1 < PPathLength2) {
            //s--->a,b,c--->t is shorter
            if (((spBit1 | currentPathBit).count() == QueryBit.count())) {
                //cout<<"s--->a,b,c--->t is shorter"<<endl;
                int startNode = vvNodePath[topPathID][0];
                int endNode = vvNodePath[topPathID][vvNodePath[topPathID].size() - 1];
                if (spBit1.count() == currentPathBit.count()) {
                    //cout<<"sp feasible"<<endl;
                    feasiblePath.push_back(s);
                    feasiblePath.insert(feasiblePath.end(), H2HSBegin.begin() + 1, H2HSBegin.end() - 1);
                    feasiblePath.push_back(currentPath[0]);
                    feasiblePath.push_back(currentPath.back());
                    feasiblePath.insert(feasiblePath.end(), H2HTEnd.begin() + 1, H2HTEnd.end() - 1);
                    feasiblePath.push_back(t);
                    feasiblePathLength = PruneRepeatedPoiPath(feasiblePath);
                } else {
                    //cout<<"current path feasible"<<endl;
                    feasiblePath.push_back(s);
                    feasiblePath.insert(feasiblePath.end(), H2HSBegin.begin() + 1, H2HSBegin.end() - 1);
                    feasiblePath.insert(feasiblePath.end(), currentPath.begin(), currentPath.end());
                    feasiblePath.insert(feasiblePath.end(), H2HTEnd.begin() + 1, H2HTEnd.end() - 1);
                    feasiblePath.push_back(t);
                    feasiblePathLength = PruneRepeatedPoiPath(feasiblePath);
                }
//                if (std::find(resultPathLength.begin(), resultPathLength.end(), feasiblePathLength) !=
//                    resultPathLength.end()) {
//                    continue;
//                }
                resultPath.push_back(feasiblePath);
                resultPathLength.push_back(feasiblePathLength);
                continue;
            }
        }
        if (PPathLength2 < PPathLength1) {
            // s---(c,b,a) --->t is shorter
            if (((spBit2 | currentPathBit).count() == QueryBit.count())) {
                //cout<<"s--->c,b,a--->t is shorter"<<endl;
                int startNode = vvNodePath[topPathID][0];
                int endNode = vvNodePath[topPathID][vvNodePath[topPathID].size() - 1];
                std::reverse(currentPath.begin(), currentPath.end());
                if (spBit2.count() == currentPathBit.count()) {
                    //cout<<"sp feasible"<<endl;
                    feasiblePath.push_back(s);
                    feasiblePath.insert(feasiblePath.end(), H2HSEnd.begin() + 1, H2HSEnd.end() - 1);
                    feasiblePath.push_back(currentPath[0]);
                    feasiblePath.push_back(currentPath.back());
                    feasiblePath.insert(feasiblePath.end(), H2HTBegin.begin() + 1, H2HTBegin.end() - 1);
                    feasiblePath.push_back(t);
                    feasiblePathLength = PruneRepeatedPoiPath(feasiblePath);
                } else {
                    //cout<<"current path feasible"<<endl;
                    feasiblePath.clear();
                    feasiblePath.push_back(s);
                    feasiblePath.insert(feasiblePath.end(), H2HSEnd.begin() + 1, H2HSEnd.end() - 1);
                    feasiblePath.insert(feasiblePath.end(), currentPath.begin(), currentPath.end());
                    feasiblePath.insert(feasiblePath.end(), H2HTBegin.begin() + 1, H2HTBegin.end() - 1);
                    feasiblePath.push_back(t);
                    feasiblePathLength = PruneRepeatedPoiPath(feasiblePath);
                }
//                if (std::find(resultPathLength.begin(), resultPathLength.end(), feasiblePathLength) !=
//                    resultPathLength.end()) {
//                    continue;
//                }
                resultPath.push_back(feasiblePath);
                resultPathLength.push_back(feasiblePathLength);
                continue;
            }
        }
        bitset<KIND> uncoverKW = currentPathBit ^ QueryBit;
        int kMD = 5;
        vector<int> h2hPath;
        vector<bitset<KIND>> H2HBit;
        for (auto kw: QueryWord) {
            if (uncoverKW.test(kw)) {
                for (auto node: keyNodeMap[kw]) {
                    num = vvNodePath[topPathID].size();
                    vvPPathLength[topPathID] = PruneRepeatedPoiPath(vvNodePath[topPathID]);
                    int dis1 = distanceQuery(node + 1, vvNodePath[topPathID][0] + 1) + vvPPathLength[topPathID];
                    int dis2 = distanceQuery(vvNodePath[topPathID].back() + 1, node + 1) + vvPPathLength[topPathID];
                    if (dis1 < (lengthThreshold * num) || dis2 < (lengthThreshold * num)) {
                        int insertStart = 0;
                        if (dis1 < dis2) {
                            insertStart = 1;
                        }
                        vector<int> newPPath;
                        bitset<KIND> newPPahBit(vvPPathBit[topPathID]);
                        if (insertStart == 1) {
                            newPPath.push_back(node);
                            newPPath.insert(newPPath.end(), currentPath.begin(), currentPath.end());
                        } else {
                            newPPath.insert(newPPath.end(), currentPath.begin(), currentPath.end());
                            newPPath.push_back(node);
                        }
                        newPPahBit |= NodesBit[node] & uncoverKW;
                        if (std::find(vvNodePath.begin(), vvNodePath.end(), newPPath) == vvNodePath.end()) {
                            vvPPathLength.push_back(PruneRepeatedPoiPath(newPPath));
                            vvNodePath.push_back(newPPath);
                            vvPPathBit.push_back(newPPahBit);
                            int startNode = vvNodePath[vvNodePath.size() - 1][0];
                            int endNode = vvNodePath[vvNodePath.size() - 1].back();
                            int newSPDis1 = distanceQuery(s + 1, startNode + 1) + vvPPathLength.back() +
                                            distanceQuery(endNode + 1, t + 1);
                            int newSPDis2 = distanceQuery(startNode + 1, t + 1) + vvPPathLength.back() +
                                            distanceQuery(endNode + 1, s + 1);
                            double ita = 0;
                            if (newSPDis1 < newSPDis2) {
                                ita = vvPPathLength.back() / ((newPPahBit.count()) * 1.0);
                            } else {
                                ita = vvPPathLength.back() / ((newPPahBit.count()) * 1.0);
                            }
                            int newSPDis = (newSPDis1 < newSPDis2) ? newSPDis1 : newSPDis2;
                            //cout<<newSPDis<<endl;
                            int LB = 0;
                            qPath.update(vvNodePath.size() - 1, newSPDis);
                            //qPath.update(vvNodePath.size() - 1, ita);
                        }
                    }
                    //}
                }
            }
        }
    }
    int Sum = 0;
    if (resultPath.empty()) {
        resultPathLength.push_back(GLength);
    }
    int pathSum = 0;
    tt2 = std::chrono::high_resolution_clock::now();
    time_span2 = std::chrono::duration_cast<std::chrono::duration<double> >(tt2 - tt1);
    //cout << "Method time is " << time_span2.count() << endl;
}


void Graph::Top_k_SKORP_Plus(int s,int t,int K,int lengthThreshold,bitset<KIND> uncover,vector<int> POICandidate,int GLength,
                           vector<int> &resultPathLength,vector<vector<int>> &resultPath){
    map<int,vector<int>> pathTable;
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::high_resolution_clock::time_point tt1;
    std::chrono::high_resolution_clock::time_point tt2;
    std::chrono::duration<double> time_span1(0.0);
    std::chrono::duration<double> time_span2(0.0);
    t1 = std::chrono::high_resolution_clock::now();
    tt1 = std::chrono::high_resolution_clock::now();
    int OD= distanceQuery(s+1,t+1);

    //int OD= (int)EuclideanDistance(s,t);
    //int OD=GLength;
    vector<int> H2Hpath;
    vector<int> H2HEdge;
    vector<bitset<KIND>> H2HPathBit;
    H2HPath(s,t,H2Hpath,H2HPathBit);
    resultPathLength.clear();
    resultPath.clear();
    int UB=GLength;
    if((H2HPathBit.back()&QueryBit).count()==QueryBit.count()){
        //cout<<"shortest path is feasible path"<<endl;
        resultPath.push_back(H2Hpath);
        resultPathLength.push_back(distanceQuery(s+1,t+1));
        return;
    }
    vector<vector<int>> vvNodePath;
    vector<int> vvPPathLength;
    vector<int> pathParent;
    vector<bitset<KIND>> vvPPathBit;
    map<int,vector<int>> keyNodeMap;

    keyNodeMap.clear();
    unordered_map<int,vector<int>> PartialPath;//value save path id in vvNodePath
//    for(auto node:POICandidate){
//        PartialPath[(NodesBit[node]&QueryBit).count()].push_back({NodesBit[node]&QueryBit,{node}});
//    }
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

    //benchmark::heap<2,int,int> qPathGrid(nodeNum);
    int count=0;

    int topPathID=0;
    int topPathDis=0;
    //qPathGrid.extract_min(topPathID ,topPathDis);
    int pNum=0;
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
                            if(std::find(vvNodePath.begin(), vvNodePath.end(),PPath)==vvNodePath.end()){
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
                                PartialPath[pathBit.count()].push_back(vvNodePath.size()-1);
                                pathParent.push_back(-1);
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
    }
    int NUM=0;
    int num=0;
    t2 = std::chrono::high_resolution_clock::now();
    time_span1 = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    cout<<"Init time is "<<time_span1.count()<<endl;
    vector<int> disVector;
    while(!qPath.empty()&&resultPath.size()<K){
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
        bitset<KIND> spBit1=QueryBit&(H2HSBeginBit.back()|H2HTEndBit.back());
        bitset<KIND> spBit2=QueryBit&(H2HSEndBit.back()|H2HTBeginBit.back());

        int PPathLength1=H2HD1+vvPPathLength[topPathID]+
                         H2HD4;
        int PPathLength2=H2HD2+vvPPathLength[topPathID]+
                         H2HD3;

        if(PPathLength1<PPathLength2&&PPathLength1){
            //cout<<"Case1"<<endl;
            //current path is a feasible path
            if(((spBit1|currentPathBit).count()==QueryBit.count())){
                int  startNode=vvNodePath[topPathID][0];
                int endNode=vvNodePath[topPathID][vvNodePath[topPathID].size()-1];
                if(spBit1.count()==currentPathBit.count()){
                    //cout<<"sp feasible"<<endl;
                    vector<int> feasiblePath;
                    feasiblePath.push_back(s);
                    feasiblePath.insert(feasiblePath.end(),H2HSBegin.begin()+1,H2HSBegin.end()-1);
                    feasiblePath.push_back(currentPath[0]);
                    feasiblePath.push_back(currentPath.back());
                    feasiblePath.insert(feasiblePath.end(),H2HTEnd.begin()+1,H2HTEnd.end()-1);
                    feasiblePath.push_back(t);
                    int feasiblePathLength= PruneRepeatedPoiPath(feasiblePath);
                    //if(std::find(resultPathLength.begin(), resultPathLength.end(),feasiblePathLength)==resultPathLength.end()){
                    //PruneRepeatedPoiPath(feasiblePath);
                    resultPath.push_back(feasiblePath);
                    resultPathLength.push_back(feasiblePathLength);
                    if(feasiblePathLength<UB){
                        UB=feasiblePathLength;
                    }
                    //}
                } else{
                    //cout<<"current path feasible"<<endl;
                    vector<int> feasiblePath;
                    feasiblePath.push_back(s);
                    feasiblePath.insert(feasiblePath.end(),H2HSBegin.begin()+1,H2HSBegin.end()-1);
                    feasiblePath.insert(feasiblePath.end(),currentPath.begin(),currentPath.end());
                    feasiblePath.insert(feasiblePath.end(),H2HTEnd.begin()+1,H2HTEnd.end()-1);
                    feasiblePath.push_back(t);
                    int feasiblePathLength=getTempPathDistance(feasiblePath);
                    //if(std::find(resultPathLength.begin(), resultPathLength.end(),feasiblePathLength)==resultPathLength.end()){
                    //PruneRepeatedPoiPath(feasiblePath);
                    resultPath.push_back(feasiblePath);
                    resultPathLength.push_back(feasiblePathLength);
                    if(feasiblePathLength<UB){
                        UB=feasiblePathLength;
                    }
                    //}
                    continue;
                }

            }
        }
        if(PPathLength2<PPathLength1){
            //cout<<"Case 2"<<endl;
            if(((spBit2|currentPathBit).count()==QueryBit.count())){
                int  startNode=vvNodePath[topPathID][0];
                int endNode=vvNodePath[topPathID][vvNodePath[topPathID].size()-1];
                std::reverse(currentPath.begin(), currentPath.end());
                if(spBit2.count()==currentPathBit.count()){
                    //cout<<"sp feasible"<<endl;
                    vector<int> feasiblePath;
                    feasiblePath.push_back(s);
                    feasiblePath.insert(feasiblePath.end(),H2HSEnd.begin()+1,H2HSEnd.end()-1);
                    feasiblePath.push_back(currentPath[0]);
                    feasiblePath.push_back(currentPath.back());
                    feasiblePath.insert(feasiblePath.end(),H2HTBegin.begin()+1,H2HTBegin.end()-1);
                    feasiblePath.push_back(t);
                    int feasiblePathLength= PruneRepeatedPoiPath(feasiblePath);
                    if(std::find(resultPathLength.begin(), resultPathLength.end(),feasiblePathLength)==resultPathLength.end()){
                        //PruneRepeatedPoiPath(feasiblePath);
                        resultPath.push_back(feasiblePath);
                        resultPathLength.push_back(feasiblePathLength);
                        if(feasiblePathLength<UB){
                            UB=feasiblePathLength;
                        }
                    }
                } else{
                    //cout<<"current path feasible"<<endl;
                    vector<int> feasiblePath;
                    feasiblePath.clear();
                    feasiblePath.push_back(s);
                    feasiblePath.insert(feasiblePath.end(),H2HSEnd.begin()+1,H2HSEnd.end()-1);
                    feasiblePath.insert(feasiblePath.end(),currentPath.begin(),currentPath.end());
                    feasiblePath.insert(feasiblePath.end(),H2HTBegin.begin()+1,H2HTBegin.end()-1);
                    feasiblePath.push_back(t);
                    int feasiblePathLength= PruneRepeatedPoiPath(feasiblePath);
                    if(std::find(resultPathLength.begin(), resultPathLength.end(),feasiblePathLength)==resultPathLength.end()){
                        //PruneRepeatedPoiPath(feasiblePath);
                        resultPath.push_back(feasiblePath);
                        resultPathLength.push_back(feasiblePathLength);
                        if(feasiblePathLength<UB){
                            UB=feasiblePathLength;
                        }
                    }
                }
                continue;
            }
        }
        if(min(PPathLength1,PPathLength2)>GLength){
            //cout<<"GLength Pruned"<<endl;
            continue;
        }
        if(min(PPathLength1,PPathLength2)>UB){
            //cout<<"UB Pruned"<<endl;
            continue;
        }
        bitset<KIND> uncoverKW=currentPathBit^QueryBit;
        //cout<<"uncover count is "<<uncoverKW.count()<<endl;

        int uncoverNum=uncoverKW.count();
        //first ti find in map
        int firstCount=0;
        if(!PartialPath[uncoverNum].empty()){
            for(auto data:PartialPath[uncoverKW.count()]){
                if(((vvPPathBit[data]|currentPathBit)&QueryBit).count()==QueryBit.count()){
                    vector<int> pathA;
                    pathA.insert(pathA.end(),currentPath.begin(),currentPath.end());
                    pathA.insert(pathA.end(),vvNodePath[data].begin(),vvNodePath[data].end());
                    vector<int> pathB;
                    pathB.insert(pathB.end(),vvNodePath[data].begin(),vvNodePath[data].end());
                    pathB.insert(pathB.end(),currentPath.begin(),currentPath.end());
                   int ADis= getTempPathDistance(pathA)+ distanceQuery(s+1,pathA[0]+1)+distanceQuery(pathA.back()+1,t+1);
                   int BDis= getTempPathDistance(pathB)+distanceQuery(s+1,pathB[0]+1)+distanceQuery(pathB.back()+1,t+1);
                   if(min(ADis,BDis)>GLength){
                       continue;
                   }
                   if(ADis<BDis){
                       resultPath.push_back(pathA);
                       resultPathLength.push_back(ADis);
                       if(ADis<UB){
                           UB=ADis;
                       }
                   } else{
                       resultPath.push_back(pathB);
                       resultPathLength.push_back(BDis);
                       if(BDis<UB){
                           UB=BDis;
                       }
                   }
                }
            }
        }
        if(resultPath.size()==K){
            break;
        }
            for (auto kw: QueryWord) {
                if (uncoverKW.test(kw)) {
                    for (auto node: keyNodeMap[kw]) {
                        num = vvNodePath[topPathID].size();
                        vvPPathLength[topPathID] = PruneRepeatedPoiPath(vvNodePath[topPathID]);
                        int dis1 = distanceQuery(node + 1, vvNodePath[topPathID][0] + 1) + vvPPathLength[topPathID];
                        int dis2 = distanceQuery(vvNodePath[topPathID].back() + 1, node + 1) + vvPPathLength[topPathID];
                        if (dis1 < (lengthThreshold * num * 0.3) || dis2 < (lengthThreshold * num * 0.3)) {
                            int insertStart = 0;
                            if (dis1 < dis2) {
                                insertStart = 1;
                            }
                            vector<int> newPPath;
                            bitset<KIND> newPPahBit(vvPPathBit[topPathID]);
                            if (insertStart == 1) {
                                newPPath.push_back(node);
                                newPPath.insert(newPPath.end(), currentPath.begin(), currentPath.end());
                                vvPPathLength.push_back(PruneRepeatedPoiPath(newPPath));
                            } else {
                                newPPath.insert(newPPath.end(), currentPath.begin(), currentPath.end());
                                newPPath.push_back(node);
                                vvPPathLength.push_back(PruneRepeatedPoiPath(newPPath));
                            }
                            newPPahBit |= NodesBit[node] & uncoverKW;
                            if (std::find(vvNodePath.begin(), vvNodePath.end(), newPPath) == vvNodePath.end()) {
                                vvNodePath.push_back(newPPath);
                                vvPPathBit.push_back(newPPahBit);
                                int startNode = vvNodePath[vvNodePath.size() - 1][0];
                                int endNode = vvNodePath[vvNodePath.size() - 1].back();
//                        int newSPDis1=getMaxLB(newPPath,newPPahBit,startNode,s)+vvPPathLength.back()+
//                                getMaxLB(newPPath,newPPahBit,endNode,t);
//                        int newSPDis2= getMaxLB(newPPath,newPPahBit,endNode,s)+vvPPathLength.back()+
//                                getMaxLB(newPPath,newPPahBit,startNode,t);
                                int newSPDis1 = distanceQuery(s + 1, startNode + 1) + vvPPathLength.back() +
                                                distanceQuery(endNode + 1, t + 1);
                                int newSPDis2 = distanceQuery(startNode + 1, t + 1) + vvPPathLength.back() +
                                                distanceQuery(endNode + 1, s + 1);
                                double ita = 0;
                                if (newSPDis1 < newSPDis2) {
                                    ita = vvPPathLength.back() / ((newPPahBit.count()) * 1.0);
                                } else {
                                    ita = vvPPathLength.back() / ((newPPahBit.count()) * 1.0);
                                }
                                int newSPDis = (newSPDis1 < newSPDis2) ? newSPDis1 : newSPDis2;
                                if(newSPDis>UB&&newPPahBit.count()!=QueryBit.count()){
                                    continue;
                                }
                                //cout<<newSPDis<<endl;

                                //method1
                                //qPath.update(vvNodePath.size()-1,ita*1000);
                                //method2
                                //disVector.push_back(newSPDis);
                                int LB=0;
                                if(newSPDis1==newSPDis){
                                     LB=getMaxLB(newPPath,newPPahBit,keyNodeMap,newPPath[0],s);
                                    qPath.update(vvNodePath.size()-1,LB+vvPPathLength.back()+distanceQuery(endNode + 1, t + 1));
                                } else{
                                    LB=getMaxLB(newPPath,newPPahBit,keyNodeMap,newPPath.back(),t);
                                    qPath.update(vvNodePath.size()-1,LB+vvPPathLength.back()+distanceQuery(startNode + 1, s + 1));
                                }

                                //qPath.update(vvNodePath.size() - 1, newSPDis);
                                PartialPath[newPPahBit.count()].push_back(vvNodePath.size()-1);
                            }
                        }
                    }
                }
        }
    }
    //cout<<topPathID<<" "<<topPathDis<<endl;

    int Sum=0;
//    if(!resultPath.empty()){
//        //resultPathLength=getTempPathDistance(resultPath);
//        resultPathLength=PruneRepeatedPoiPath(resultPath);
//    } else{
//        resultPathLength=GLength;
//        cout<<"error"<<endl;
//    }
    if(resultPath.empty()){
        resultPathLength.push_back(GLength);
    }
    int pathSum=0;
    for(int i=0;i<resultPathLength.size();i++){
        pathSum+=resultPathLength[i];
    }
    //resultPathLength.push_back(pathSum/resultPath.size());
    tt2 = std::chrono::high_resolution_clock::now();
    time_span2 = std::chrono::duration_cast<std::chrono::duration<double> >(tt2 - tt1);
    cout<<"Method time is "<<time_span2.count()<<endl;
}