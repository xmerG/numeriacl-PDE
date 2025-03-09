#ifndef _CIRCLE_
#define _CIRCLE_
#include<iostream>
#include<cmath>
using namespace std;
class Circle{
private:
    double x0=0.0;
    double y0=0.0;
    double radius=0.0;
public:
    Circle();
    Circle(double x, double y, double r);
    bool inCircle(double x, double y) const;

    double get_radius() const;

    double getX() const;

    double getY()const;

    double distance(double x, double y) const;

    double x_distance_to_circle(double x, double y) const;  //有向距离

    double y_distance_to_circle(double x, double y) const; //有向距离

    double angle_x_direction(double x) const; 

    double angle_y_direction(double y) const;
};


#endif