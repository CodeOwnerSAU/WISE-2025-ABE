//
// Created by xiongxing on 3/25/24.
//

#include "MBR.h"

void MBR::unionWith(Point p) {
    minx = MIN(minx, p.coord[0]);
    maxx = MAX(maxx, p.coord[0]);
    miny = MIN(miny, p.coord[1]);
    maxy = MAX(maxy, p.coord[1]);
}
int MBR::area()
{
    return (maxx - minx) * (maxy - miny);
}
int MBR::width()
{
    return maxx - minx;
}
int MBR::height()
{
    return maxy - miny;
}
double MBR::diagonal()
{
    return sqrt((maxy - miny)*(maxy - miny) + (maxx - minx)*(maxx - minx));
}
