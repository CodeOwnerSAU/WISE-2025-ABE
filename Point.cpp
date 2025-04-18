//
// Created by xiongxing on 3/25/24.
//

#include <cmath>
#include "Point.h"
Point::Point() {
    id = -1;
}

Point::Point(int _x, int _y) {
    coord[0] = _x; coord[1] = _y;
}

Point::Point(int _id, int _x, int _y) { _id = _id; coord[0] = _x; coord[1] = _y; }

double Point::Dist(Point point)
{
    double div = 1000000.0;
    double x1lon = ((double) coord[0]) / div;
    double y1lat = ((double) coord[1]) / div;
    double x2lon = ((double) point.coord[0]) / div;
    double y2lat = ((double) point.coord[1]) / div;
    return distanceEarth(y1lat, x1lon, y2lat, x2lon) * 10000;
}

double Point::deg2rad(double deg) {
    return (deg * M_PI / 180);
}

double Point::rad2deg(double rad) {
    return (rad * 180 / M_PI);
}

double Point::distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(lat1d);
    lon1r = deg2rad(lon1d);
    lat2r = deg2rad(lat2d);
    lon2r = deg2rad(lon2d);
    u = sin((lat2r - lat1r)/2);
    v = sin((lon2r - lon1r)/2);
    return 2.0 * 6371.0 * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}

Point::~Point()
{
}