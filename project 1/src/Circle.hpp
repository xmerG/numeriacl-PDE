#ifndef _CIRCLE_
#define _CIRCLE_
#include<iostream>
#include<cmath>
using namespace std;
class Circle{
private:
    double x0=0.0;
    double y0=0.0;
    double radius-0.0;
public:
    Circle();
    Circle(double x, double y);
    bool inCircle(double x, double y) const;

    double get_radius() const;

    double getX() const;

    double getY()const;

    double distance(double x, double y) const;

    double x_distance_to_circle(double x, double y) const;

    double y_distance_to_circle(double x, double y) const;
};


#endif