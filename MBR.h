//
// Created by xiongxing on 3/25/24.
//

#ifndef GRIDINDEX_MBR_H
#define GRIDINDEX_MBR_H


#include "Point.h"
#include <cmath>
#include <sys/param.h>
class MBR {
 public:
    int minx = 999999999;
    int miny = 999999999;
    int maxx = -999999999;
    int maxy = -999999999;
    MBR(){};
    MBR(int _minx, int _miny, int _maxx, int _maxy);
    MBR(Point & minPoint, Point & maxPoint);
    void unionWith(Point _point);
    int area();
    int width();
    int height();
    double diagonal();
    ~MBR(){};
};


#endif //GRIDINDEX_MBR_H
