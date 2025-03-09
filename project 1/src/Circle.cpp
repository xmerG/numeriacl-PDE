#include"Circle.hpp"
using namespace std;

Circle::Circle(){}

Circle::Circle(double x, double y, double r):x0(x), y0(y), radius(r){}

double Circle::get_radius() const{
    return radius;
}

double Circle::getX() const{
    return x0;
}

double Circle::getY() const{
    return y0;
}

bool Circle::inCircle(double x, double y) const{
    if(!(this->distance(x, y) <radius)){
        return false;
    }
    else{
        return true;
    }
}

double Circle::distance(double x, double y) const{
    return sqrt(pow(x0-x, 2) + pow(y0-y, 2));
}

double Circle::x_distance_to_circle(double x, double y) const{
    if(this->inCircle(x,y)){
        cerr<<"can't calculate cause the point is in circle"<<endl;
        return -1;
    }
    else if(y>=y0-radius && y<=y0+radius){
        if(x<x0){
            return x0-x-sqrt(pow(radius,2)-pow(y-y0, 2));
        }
        else{
            return x0-x+sqrt(pow(radius, 2)-pow(y-y0, 2));
        }
    }
    else{
        return 2.0;
    }
}

double Circle::y_distance_to_circle(double x, double y) const{
    if(this->inCircle(x,y)){
        cerr<<"can't calculate cause the point is in circle"<<endl;
        return -1;
    }
    else if(x>=x0-radius && x<=x0+radius){
        if(y<y0){
            return y0-y-sqrt(pow(radius, 2)-pow(x-x0, 2));
        }
        else{
            return y0-y+sqrt(pow(radius, 2)-pow(x-x0, 2));
        }
    }
    else{
        return 2.0;
    }
}

double Circle::angle_x_direction(double x) const{
    return abs(x-x0)/radius;
}

double Circle::angle_y_direction(double y) const{
    return abs(y-y0)/radius;
}