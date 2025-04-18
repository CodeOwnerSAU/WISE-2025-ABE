//
// Created by xiongxing on 3/25/24.
//

#ifndef GRIDINDEX_POINT_H
#define GRIDINDEX_POINT_H


class Point {
public:
    int id;
    int coord[2];
    int nodeId = -1;
    Point();
    Point(int _x, int _y);
    Point(int _id, int _x, int _y);
    ~Point();
    double Dist(Point _point);
    double deg2rad(double deg);
    double rad2deg(double rad);
    double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d);

};


#endif //GRIDINDEX_POINT_H
