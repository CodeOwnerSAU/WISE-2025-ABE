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
#define K 1
#define QueryNum 6
#define TestNum 100
string DisFile;
void DJMT(){
    graphFile = "/home/xiongxing/reviewProject/Test_MAP/NY/USA-road-d.NY.gr";
    H2HIndexFile ="/home/xiongxing/reviewProject/Test_MAP/NY/H2H/NY.index";
    edgeFile="/home/xiongxing/reviewProject/Test_MAP/NY/USA-road-d.NY.gr";
    //nodeFile="/home/xiongxing/reviewProject/Test_MAP/CAL/POI/D8/CAL_keys_20_uniform_500.node";
    nodeFile="/home/xiongxing/reviewProject/buildSubNet/DJM/djm_keywords_4_8/NY_keys_20_uniform_500.node";
    QueryTest="/home/xiongxing/reviewProject/buildSubNet/DJM/djm_keywords_4_8/keywords_4/query.txt";
}
void readRoadNetwork(){
    //data file
    edgeFile="/home/xiongxing/reviewProject/Test_MAP/NY/USA-road-d.NY.gr";
    nodeFile="/home/xiongxing/reviewProject/Test_MAP/NY/POI/D8/NY_keys_20_uniform_500.node";
    //nodeFile="/home/xiongxing/reviewProject/buildSubNYt/DJM/POI/NY_keys_20_uniform_100.node";
    H2HIndexFile ="/home/xiongxing/reviewProject/Test_MAP/NY/H2H/NY.index";
    IGTreeIndex = "/home/xiongxing/reviewProject/Test_MAP/NY/IGTree/NY.IGtree";
    IGTreeMIND = "/home/xiongxing/reviewProject/Test_MAP/NY/IGTree/NY.minds";
    IGTreePath = "/home/xiongxing/reviewProject/Test_MAP/NY/IGTree/NY.paths";
    //testFile
    testFile="/home/xiongxing/reviewProject/Test_MAP/NY/query.node";
    resultFile="/home/xiongxing/reviewProject/SY_SKORP/SKORP.result";
    DisFile="/home/xiongxing/reviewProject/SY_SKORP/SKORP_Dis_OD.result";
    //QueryTest="/home/xiongxing/reviewProject/buildSubNet/DJM/IG-Tree/IG-Tree_5_H2H_3.node";
}
void bug(){
    ifstream TWE_Base("/home/xiongxing/reviewProject/SY_SKORP/TWE_Base.init");
    ifstream TWE("/home/xiongxing/reviewProject/SY_SKORP/TWE.init");
    vector<pair<int,int>> TWEBase;
    vector<pair<int,int>> TWEIG;
    int a,b;
    while (TWE_Base>>a>>b){
        pair<int,int> A={a,b};
        if(std::find(TWEBase.begin(), TWEBase.end(),A)!=TWEBase.end()){
            cout<<A.first<<" "<<A.second<<endl;
        }
        TWEBase.emplace_back(a,b);
    }
    while (TWE>>a>>b){
        pair<int,int> A={a,b};
        if(std::find(TWEIG.begin(), TWEIG.end(),A)!=TWEIG.end()){
            cout<<A.first<<" "<<A.second<<endl;
        }
        TWEIG.emplace_back(a,b);
    }
    cout<<TWEBase.size()<<endl;
    cout<<TWEIG.size()<<endl;
    int count=0;
    for(auto path:TWEBase){
        if(std::find(TWEIG.begin(), TWEIG.end(),path)==TWEIG.end()){
            cout<<path.first<<" "<<path.second<<endl;
        }
    }
    cout<<count<<endl;
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
    int aveTotalLength;
    vector<int> Query(test.begin(),test.begin()+QueryNum);
    int s,t,dis;
    int count=0;
    for(auto q:Query){
        cout<<q<<" ";
    }
    cout<<endl;
    vector<vector<int>> kPaths;
    ofstream ofs("/home/xiongxing/reviewProject/SY_SKORP/TWE+maxLB.init");
    vector<int> kPathLength;
    long long pairSum=0;
    int total_dis=0;
    while(testOD>>s>>t){
//        s=35173;
//        t=176941;
        total_dis=0;
        kPaths.clear();
        kPathLength.clear();
        G.initPairSum=0;
        if(count==100){
            break;
        }
//        s=116841;
//        t=167192;
        G.InitQueryInfo(s,t,Query);
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
        aveTotalLength+=total_dis/kPathLength.size();

        //cout<<"current ave length is "<<aveTotalLength/count<<endl;
    }
    cout<<"Query path num is "<<K<<endl;
    cout<<"total ave length is "<<aveTotalLength/count<<endl;
    cout<<"ave time is "<<(time_totol.count()/count)*1000<<" ms"<<endl;
    //cout<<"ave pair sum "<<pairSum/count<<endl;
}
void SKORP(){
    readRoadNetwork();
    Graph G=Graph(edgeFile,nodeFile);
    G.Init_IGTree_Index(IGTreeIndex,IGTreeMIND,IGTreePath);
    G.Init_H2H_Index(H2HIndexFile, edgeFile);
    ifstream testOD(testFile);
    ofstream result(resultFile);
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::high_resolution_clock::time_point t3;
    std::chrono::high_resolution_clock::time_point t4;
    std::chrono::duration<double> time_span(0.0);
    std::chrono::duration<double> IG_base_time_totol(0.0);
    vector<int> Query(test.begin(),test.begin()+QueryNum);
    int GPathSum=0;
    int SYPathSum=0;
    int s,t,dis;
    int count=0;
    for(auto q:Query){
        cout<<q<<" ";
    }
    cout<<endl;
    G.test=0;
    //fstream testFileDJM(QueryTest);
    ofstream ofs(resultFile);
    ofstream CPOI("/home/xiongxing/reviewProject/SY_SKORP/CPOI.result");
    int C=0;
    while(testOD>>s>>t){
//        s=10641;
//        t=4395;
        if(count==TestNum){
            break;
        }
        G.InitQueryInfo(s,t,Query);
        cout<<count<<": "<<s<<" "<<t<<" "<<G.distanceQuery(s+1,t+1)<<endl;
        vector<int> POICandidate;
        t1 = std::chrono::high_resolution_clock::now();
        vector<vector<int>> Kpath;
        vector<int> KpathDis;
        //G.TWE(s,t,Query,Kpath,KpathDis);
        G.getPOICoordinate(s,t,Query,POICandidate);
        vector<int> Gpath1;
        vector<int> Gpath2;
        vector<int> Gpath;
        vector<int> POICandidate_Greedy_Filter;
        //G.GreedyFiler(s,t,POICandidate,POICandidate_Greedy_Filter);
        vector<int> POICandidate_Triangle_Filter;
        //POICandidate_Triangle_Filter=G.RPOI;
        G.TriangleFiltering(s,t,POICandidate,POICandidate_Triangle_Filter);
        Gpath1.clear();
        Gpath2.clear();
        int A1;
        int B2;
        int Length1=G.buildGreedyPath(s,t,POICandidate_Triangle_Filter,Gpath1,A1);
        int Length2=G.buildGreedyPath(t,s,POICandidate_Triangle_Filter,Gpath2,B2);
        int pathEdge;int Length;
        if(Length1<Length2){
            pathEdge=Gpath1.size()-1;
            Gpath=Gpath1;
            Length=Length1;
        }else{
            pathEdge=Gpath2.size()-1;
            Gpath=Gpath2;
            Length=Length2;
        }
        POICandidate.clear();
        for(auto node:POICandidate_Triangle_Filter){
            if(G.distanceQuery(s+1,node+1)+G.distanceQuery(node+1,t+1)<Length){
                POICandidate.push_back(node);
            }
        }
        double aveEdgeDis=0;
        if(Gpath.size()>3){
            //aveEdgeDis=Length/pathEdge;
            aveEdgeDis=(Length
                        -G.distanceQuery(Gpath[0]+1,Gpath[1]+1)
                        -G.distanceQuery(Gpath[Gpath.size()-2]+1,Gpath[Gpath.size()-1]+1))/(pathEdge-2);
        }else{
            aveEdgeDis=Length/2;
        }
        vector<vector<int>> PPathSet;
        vector<bitset<KIND>> PPathBit;
        vector<int> resultPath;
        int resultPathLength;
        bitset<KIND> resultPathBit;
        //t1 = std::chrono::high_resolution_clock::now();
        vector<vector<int>> top_K_result;
        vector<int> top_K_result_Length;
        int KPathSum=0;
        //G.GreedyFiler(s,t,POICandidate_Triangle_Filter,POICandidate_Greedy_Filter);
        G.Top_k_SKORP_V1(s,t,K,aveEdgeDis,G.QueryBit,POICandidate,Length,top_K_result_Length,top_K_result);
        //cout<<top_K_result_Length[0]<<endl;
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
        IG_base_time_totol+=time_span;
        bitset<KIND> resultBit;
        ofs<<s<<" "<<t<<" "<<G.distanceQuery(s+1,t+1)<<" "<<"time "<<time_span.count()<<endl;
        for(int i=0;i<top_K_result.size();i++){
            ofs<<"Dis : "<<top_K_result_Length[i]<<endl;
            for(auto node:top_K_result[i]){
                if((G.NodesBit[node]&G.QueryBit).count()>0||node==s||node==t)
                    ofs<<node<<" ";
            }
            ofs<<endl;
            KPathSum+=top_K_result_Length[i];
            cout<<top_K_result_Length[i]<<endl;
        }
        ofs<<endl;
        int topKSum=0;
        count++;
        SYPathSum+=KPathSum/top_K_result_Length.size();
        cout<<"Top-k ave dis:"<<KPathSum/top_K_result_Length.size()<<endl;
        cout<<"time is "<<time_span.count()<<endl;
        cout<<"total time is "<<IG_base_time_totol.count()<<endl;
        cout<<"--------------"<<endl;
        cout<<endl;
        //cout<<"SKORP: "<<SYPathSum<<endl;
        cout<<"current ave is "<<SYPathSum/count<<endl;
    }
    //cout<<CELLSIZE<<endl;
    cout<<"Total Dis "<<SYPathSum<<endl;
    cout<<"ave Dis is "<<SYPathSum/TestNum<<endl;
    cout<<"ave Time is "<<IG_base_time_totol.count()/TestNum<<endl;
    cout<<"Query Num "<<Query.size()<<endl;
    cout<<nodeFile<<endl;
}
void DJM(){
    DJMT();
    std::cout << "Hello, World!" << std::endl;
    //readRoadNetwork();
    Graph G=Graph(edgeFile,nodeFile);
    G.createGridIndex();
    G.Init_H2H_Index(H2HIndexFile, edgeFile);
    //G.Init_IGTree_Index(IGTreeIndex, IGTreeMIND, IGTreePath);
    ifstream testOD(testFile);
    ofstream result(resultFile);
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::high_resolution_clock::time_point t3;
    std::chrono::high_resolution_clock::time_point t4;
    std::chrono::duration<double> time_span(0.0);
    std::chrono::duration<double> IG_base_time_totol(0.0);
    //vector<int> Query(test.begin(),test.begin()+6);
    int GPathSum=0;
    int SYPathSum=0;
    int s,t;
    int count=0;
    cout<<endl;
    G.test=0;
    fstream testFileDJM(QueryTest);
    string line;
    int ID1, ID2,H2HNode;
    while(getline(testFileDJM,line)){
        count++;
        std::istringstream iss(line);
        iss>>s>>t>>H2HNode;
        std::string temp;
        iss >> temp;
        int pos=temp.find(":");
        int Q =stoi(temp.substr(pos+1,temp.length()));
        vector<int> Query;
        Query.push_back(Q);
        int number;
        for (int i = 0; i < QueryNum-1; ++i) {
            iss >> number;
            Query.push_back(number);
        }
        for(auto kw:Query){
            cout<<kw<<",";
        }
        cout<<endl;
//        s=261110;
//        t=1168;
//        vector<int> Temp={0,189,194,308,411,422,430,455};
//        Query=Temp;
        G.InitQueryInfo(s,t,Query);
        cout<<count<<" "<<s<<" "<<t<<" "<<G.distanceQuery(s+1,t+1)<<endl;
        vector<int> POICandidate;
        t1 = std::chrono::high_resolution_clock::now();
        G.getPOICoordinate(s,t,Query,POICandidate);
        G.getPOICoordinateByGrid(s,t,Query,POICandidate);
        //G.getPOICoordinate(s,t,Query,POICandidate);
        vector<int> Gpath1;
        vector<int> Gpath2;
        vector<int> Gpath;
        vector<int> POICandidate_Greedy_Filter;
        //G.GreedyFiler(s,t,POICandidate,POICandidate_Greedy_Filter);
        vector<int> POICandidate_Triangle_Filter;
        G.TriangleFiltering(s,t,POICandidate,POICandidate_Triangle_Filter);
        Gpath1.clear();
        Gpath2.clear();
        vector<vector<int>> Kpath;
        vector<int> KpathDis;
        int D1;
        int D2;
        int Length1=G.buildGreedyPath(s,t,POICandidate_Triangle_Filter,Gpath1,D1);
        int Length2=G.buildGreedyPath(t,s,POICandidate_Triangle_Filter,Gpath2,D2);
        int pathEdge;int Length;
        if(Length1<Length2){
            pathEdge=Gpath1.size()-1;
            Gpath=Gpath1;
            Length=Length1;
        } else{
            pathEdge=Gpath2.size()-1;
            Gpath=Gpath2;
            Length=Length2;
        }
        double aveEdgeDis=0;
        if(Gpath.size()>3){
            aveEdgeDis=(Length
                        -G.distanceQuery(Gpath[0]+1,Gpath[1]+1)
                        -G.distanceQuery(Gpath[Gpath.size()-2]+1,Gpath[Gpath.size()-1]+1))/(pathEdge-2);
        }else{
            aveEdgeDis=Length/2;
        }
        vector<vector<int>> PPathSet;
        vector<bitset<KIND>> PPathBit;
        vector<vector<int>> resultPath;
        vector<int> resultPathLength;
        bitset<KIND> resultPathBit;
        t1 = std::chrono::high_resolution_clock::now();
        //G.Top_k_SKORP_V1(s,t,K,aveEdgeDis,G.QueryBit,POICandidate_Triangle_Filter,Length,top_K_result_Length,top_K_result);
        G.Top_k_SKORP_V1(s,t,K,aveEdgeDis,G.QueryBit,POICandidate_Triangle_Filter,Length,resultPathLength,resultPath);
//        for(auto node:resultPath){
//            resultPathBit|=G.NodesBit[node]&G.QueryBit;
//        }
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
        IG_base_time_totol+=time_span;
        //cout<<"Grid time is "<<time_span.count()<<endl;
        cout<<resultPathLength[0]<<endl;
        //count++;
        SYPathSum+=resultPathLength[0];
        cout<<"current ave dis is "<<SYPathSum/count<<endl;
        result<<s<<" "<<t<<" "<<G.distanceQuery(s+1,t+1)<<endl;
        result<<"result path length is "<<resultPathLength[0]<<endl;
//        for(auto node:resultPath[0]){
//            result<<node<<" ";
//        }
        result<<endl;
        result<<"--------------------------------------------------"<<endl;
        result<<endl;
    }
    cout<<CELLSIZE<<endl;
    cout<<"ave Dis is "<<SYPathSum/1000<<endl;
    cout<<"ave Time is "<<IG_base_time_totol.count()/1000<<endl;
    cout<<nodeFile<<endl;
}
int main() {
    TWE();
    return 0;
}
