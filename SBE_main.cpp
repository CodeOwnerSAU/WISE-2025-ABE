#include <iostream>
#include "graph.h"
string keysFile,graphFile;
const char *H2HIndexFile;
const char *edgeFile;
const char *nodeFile;
const char *IGTreeIndex;
const char *IGTreeMIND;
const char * IGTreePath;
string testFile,resultFile;
string QueryTest;
vector<int> test={12,56,135,198,234,289,341,389,434,489};
vector<int> Qu={1,2,3,4,5,6};
#define K 10
#define QueryNum 6
#define TestNum 200
string DisFile;
void readRoadNetwork(){
    edgeFile="../USA-road-d.NY.gr";
    nodeFile="../NY_keys_20_uniform_500.node";
    H2HIndexFile ="../NY.index";
    IGTreeIndex = "../IGTree/NY.IGtree";
    IGTreeMIND = "../IGTree/NY.minds";
    IGTreePath = "../IGTree/NY.paths";
    testFile="../query.node";
    DisFile="/home/xiongxing/reviewProject/SY_SKORP/SKORP_Dis_OD.result";
}
void TWE(){
    readRoadNetwork();
    Graph G=Graph(edgeFile,nodeFile);
    G.Init_H2H_Index(H2HIndexFile, edgeFile);
    G.Init_IGTree_Index(IGTreeIndex,IGTreeMIND,IGTreePath);
    //G.createGridIndex();
    ifstream testOD(testFile);
    ofstream result(resultFile);
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::high_resolution_clock::time_point t3;
    std::chrono::high_resolution_clock::time_point t4;
    std::chrono::duration<double> time_span(0.0);
    std::chrono::duration<double> time_totol(0.0);
    double aveTotalLength;
    vector<int> Query(test.begin(),test.begin()+QueryNum);
    int s,t,dis;
    int count=0;
    for(auto q:Query){
        cout<<q<<" ";
    }
    cout<<endl;
    vector<vector<int>> kPaths;
    ofstream ofs("../result.txt");
    vector<int> kPathLength;
    long long pairSum=0;
    int total_dis=0;
   int temp;
    while(testOD>>s>>t){
        total_dis=0;
        kPaths.clear();
        kPathLength.clear();
        G.initPairSum=0;
        if(count==TestNum){
            break;
        }
        G.InitQueryInfo(s,t,Query);
        G.sp=G.distanceQuery(s+1,t+1);
        cout<<count<<": "<<s<<" "<<t<<" "<<G.distanceQuery(s+1,t+1)<<endl;
        ofs<<count<<": "<<s<<" "<<t<<" "<<G.distanceQuery(s+1,t+1)<<endl;
        t1 = std::chrono::high_resolution_clock::now();
        G.KTWE(s,t,K,Query,kPaths,kPathLength,ofs);
        t2 = std::chrono::high_resolution_clock::now();
        for(auto path:kPathLength){
            total_dis+=path;
            cout<<path<<endl;
        }
        cout<<"ave dis is "<<total_dis/kPathLength.size()<<endl;
        ofs<<"ave dis is "<<total_dis/kPathLength.size()<<endl;
        time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
        time_totol+=time_span;
        cout<<time_span.count()<<endl;
        pairSum+=G.initPairSum;
        count++;
        aveTotalLength+=(total_dis*1.0)/kPathLength.size();
    }
    cout<<"Query path num is "<<K<<endl;
    cout<<"total ave length is "<<aveTotalLength/(count*10000)<<endl;
    cout<<"ave time is "<<(time_totol.count()/count)*1000<<" ms"<<endl;
}
int main() {
    TWE();
    return 0;
}
