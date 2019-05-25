#ifndef MAIN_H
#define MAIN_H

#include <ctime>
#include <cmath>

using real = double;

using integer = unsigned int;

struct Coor {
    real x;
    real y;
    real z;
    Coor() {}
    Coor(real x, real y, real z): x(x), y(y), z(z) {}
    real norm() { return sqrt(x*x + y*y + z*z); }
    real norm2() { return x*x + y*y + z*z; }
    friend Coor operator- (const Coor& coor1, const Coor& coor2) {
        return Coor(coor1.x-coor2.x, coor1.y-coor2.y, coor1.z-coor2.z);
    }
};

#endif