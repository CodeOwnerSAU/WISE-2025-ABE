//
// Created by xiongxing on 7/23/24.
//
#include "graph.h"
void Graph::createGridIndex() {
    cout<<"Start creat Grid Index"<<endl;
    //biggest is a MBR that can include all node
    MBR* biggest=new MBR();
    for(int i=0;i<nodeNum;i++){
        biggest->unionWith(PointList[i]);
    }
    this->extend=biggest;
    //build Index by graph
    this->buildIndex();
    printf("Finish the building of grid! \n");
}

void Graph::buildIndex() {
    POIGridMap.resize(nodeNum);
    if(PointList.size()<MINPOINT){
        return;
    }
    calculateGridSize(CELLSIZE);
    colSize = ceil((double)extend->width() / colLength);
    cout<<"colSize = " << colSize <<endl;
    rowSize = ceil((double)extend->height() / rowLength);
    cout<< "rowSize = " << rowSize << endl;
    //all samll grid num is numNode
    int numNode = colSize * rowSize;//获取当前最大矩形可以划分的网格总数
    cout<< "numNode = " << numNode << endl;
    cellContain.resize(numNode);//cellContain保存了每一个网格中包含POI的数量信息
    //初始化所有网格中的关键字信息和POI信息
    for(int i = 0; i < numNode; i++)
    {
        gridNode item;
        gridNodeMap.push_back(item);
        cellContain[i] = 0;
    }
    //get each POI point 遍历每一个POI
    for (int pid = 0; pid < POIObjects.size(); pid++)
    {
        //从点集合中获取当前遍历到的POI,通过getCell函数获取他所处在的网格
        int gridid = getCell(PointList[POIObjects[pid]]);
        POIGridMap[POIObjects[pid]]=gridid;
//        cout<<"gridId is"<<gridid<<endl;
        //获取当前POI的id
        int node =PointList[POIObjects[pid]].id;
        //更新当前POI所处的网格中所包含的关键字信息
        gridNodeMap[gridid].gridNodeBits |= NodesBit[node];
        //将当前POI点放入当前网格
        gridNodeMap[gridid].POI.push_back(node);
        //网格存储的POI点数量+1
        cellContain[gridid]++;
        if(gridid >= numNode)
            printf("Error: Grid::build() occurs an error!\n");
    }
}

void Graph::calculateGridSize(const int decimeter) {
    //1000 decimeter is 100m
    int minx, miny, maxx, maxy;
    minx = extend->minx;
    miny = extend->miny;
    maxx = extend->maxx;
    maxy = extend->maxy;
    Point leftDown;
    leftDown.coord[0] = minx;
    leftDown.coord[1] = miny;
    Point leftUp;
    leftUp.coord[0] = minx;
    leftUp.coord[1] = maxy;
    Point rightDown;
    rightDown.coord[0] = maxx;
    rightDown.coord[1] = miny;
    //get biggest MBR xLength and yLength
    //单位是分米
    double xLength = leftDown.Dist(rightDown);
    double yLength = leftDown.Dist(leftUp);

    //列
    colLength = (maxx - minx) * ((double)decimeter/xLength);
    //行 length
    rowLength = (maxy - miny) * ((double)decimeter/yLength);
    cout << "minx = " << minx << "\n";
    cout << "maxx = " << maxx << "\n";
    cout << "miny = " << miny << "\n";
    cout << "maxy = " << maxy << "\n";
    cout << "xLength = " << xLength << "\n";
    cout << "yLength = " << yLength << "\n";
    cout << "colLength = " << colLength << "\n";
    cout << "rowLength = " << rowLength << "\n";
    cout << "decimeter = " << decimeter << "\n";
}

int Graph::getCell(Point &point) {
    int col = abs((point.coord[0] - extend->minx)) / colLength;
    int row = abs((point.coord[1] - extend->miny)) / rowLength;
    return  row*colSize + col;
}
int Graph::rightOf(int cid)
{
    if(cid  == -1)
        return -1;
    else
    {
        int result = cid + 1;
        return (result%colSize == 0) ? -1 : result;
    }
}

int Graph::leftOf(int cid)
{
    if(cid  == -1)
        return -1;
    else
    {
        int result = cid - 1;
        return (cid%colSize == 0) ? -1 : result;
    }
}

int Graph::upperOf(int cid)
{
    if(cid  == -1)
        return -1;
    else
    {
        int result = cid + colSize;
        return (result >= rowSize*colSize) ? -1 : result;
    }
}

int Graph::belowOf(int cid)
{
    if(cid  == -1)
        return -1;
    else
    {
        int result = cid - colSize;
        return (result<0) ? -1 : result;
    }
}
void Graph::get_POI_Straight_From_StoT(int s, int t, vector<int> &cellList) {
    // straight from s ti yo is a function as  f(x)=ax+b
    //y=kx+b;
    int dx=PointList[t].coord[0]-PointList[s].coord[0];
    int dy=PointList[t].coord[1]-PointList[s].coord[1];

    if(dx==0||dy==0){
        if(dx==0){
            dx=1;
        }else if(dy==0){
            dy=1;
        }
    }
    int increx=dx/abs(dx);
    int increy=dy/abs(dy);
    vector<int> tmp;

    double k=(dy*1.0)/(dx*1.0);
    double b=-k*(PointList[t].coord[0]*1.0)+PointList[t].coord[1]*1.0;
    bitset<KIND> cellbit;
    //bitset<KIND> uncoverbit(graph.QueryBit);
    //1 get cell straight from s to t

    for(int i=0;i<abs(dx);i+=colLength/3){
        //graph.GridextendNum++;
        int x=PointList[s].coord[0]+(increx*i);
        int y=(int)(k*x+b);
        int col=((x-extend->minx))/colLength;
        int row=((y-extend->miny))/rowLength;
        int cell=row*colSize+col;
        if(cell<0){
            cout<<cell<<endl;
        }
        tmp.push_back(cell);
        if(cellContain[cell]>0){
            //this means that cell contain POI that have keywords from queryKey
            if((gridNodeMap[cell].gridNodeBits&QueryBit).count()>0){
                cellbit|=gridNodeMap[cell].gridNodeBits&QueryBit;
                cellList.push_back(cell);
                //STCell.push_back(cell);
            }
        }
    }
    //uncoverbit^=cellbit;
    for(int i=0;i<abs(dy);i+=rowLength/3){
        //.GridextendNum++;
        int y=PointList[s].coord[1]+(increy*i);
        int x=(int)((y-PointList[s].coord[1])/k)+PointList[s].coord[0];
        int col=((x-extend->minx))/colLength;
        int row=((y-extend->miny))/rowLength;
        int cell=row*colSize+col;
        if(cell<0){
            cout<<cell<<endl;
        }
        tmp.push_back(cell);
        if(cellContain[cell]>0){
            //this means that cell contain POI that have keywords from queryKey
            if((gridNodeMap[cell].gridNodeBits&QueryBit).count()>0){
                cellbit|=gridNodeMap[cell].gridNodeBits&QueryBit;
                cellList.push_back(cell);
                //STCell.push_back(cell);
            }
        }
    }
    //uncoverbit^=cellbit;
    sort(tmp.begin(),tmp.end());
    tmp.erase(unique(tmp.begin(),tmp.end()),tmp.end());

    //2 expend from straight from s to T until all keywords has covered
    vector<int> temp;//temp save all test cell
    int flag=0;
    bitset<KIND> uncover=((cellbit^QueryBit)&QueryBit);
    while(uncover.count()!=0&&(cellbit^QueryBit).count()!=0){
        temp.clear();
        for(int i=0;i<tmp.size();i++){
            //below
            int cell= belowOf(tmp[i]);
            //.GridextendNum++;
            if(cell!=-1){
                temp.push_back(cell);
                if(cellContain[cell]>0){
                    if((gridNodeMap[cell].gridNodeBits&uncover).count()>0){
                        cellbit|=gridNodeMap[cell].gridNodeBits&uncover;
                        cellList.push_back(cell);
                    }
                }
            }
            //upper
            cell= upperOf(tmp[i]);

            //.GridextendNum++;
            if(cell!=-1){
                temp.push_back(cell);
                if(cellContain[cell]>0){
                    if((gridNodeMap[cell].gridNodeBits&uncover).count()>0){
                        cellbit|=gridNodeMap[cell].gridNodeBits&uncover;
                        cellList.push_back(cell);
                    }
                }
            }

            //left
            cell= leftOf(tmp[i]);
            //.GridextendNum++;
            if(cell!=-1){
                temp.push_back(cell);
                if(cellContain[cell]>0){
                    if((gridNodeMap[cell].gridNodeBits&uncover).count()>0){
                        cellbit|=gridNodeMap[cell].gridNodeBits&uncover;
                        cellList.push_back(cell);
                    }
                }
            }

            //right
            cell= rightOf(tmp[i]);
            //.GridextendNum++;
            if(cell!=-1){
                temp.push_back(cell);
                if(cellContain[cell]>0){
                    if((gridNodeMap[cell].gridNodeBits&uncover).count()>0){
                        cellbit|=gridNodeMap[cell].gridNodeBits&uncover;
                        cellList.push_back(cell);
                    }
                }
            }
        }
        sort(temp.begin(),temp.end());
        //delete same cell
        temp.erase(unique(temp.begin(),temp.end()),temp.end());
        tmp.clear();
        //change tmp ,tmp=temp in order to next iterator
        tmp.insert(tmp.begin(),temp.begin(),temp.end());
        sort(cellList.begin(),cellList.end());
        cellList.erase(unique(cellList.begin(),cellList.end()),cellList.end());
    }
    sort(cellList.begin(),cellList.end());
    cellList.erase(unique(cellList.begin(),cellList.end()),cellList.end());
    //cout<<"get cell candidate finished"<<endl;
}
void Graph::getPOICoordinateByGrid(int s,int t,vector<int> Query,vector<int> &POICandidate){
   POICandidate.clear();
    //delete Keywords that in ID1 and ID2
    //save cell that include keywords in queryKey
    vector<int> cellList;
    cellList.clear();
    get_POI_Straight_From_StoT(s,t,cellList);
    //save all poi node in CellList
    for(int i=0;i<cellList.size();i++){
        for(auto node:gridNodeMap[cellList[i]].POI){
            if((NodesBit[node]&QueryBit).count()>0){
                POICandidate.push_back(node);
            }
        }
        //POICandidate.insert(POICandidate.begin(),gridNodeMap[cellList[i]].POI.begin(),gridNodeMap[cellList[i]].POI.end());
    }
}

