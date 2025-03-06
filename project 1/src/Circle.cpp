#include"Circle.hpp"
using namespace std;

Circle::Circle(){}

Circle::Circle(double x, double y):x0(x), y0(y){}

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
    if(this->distance(x, y) > r){
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
    else{
        return abs(x-x0)-sqrt(pow(radius, 2)-pow((y-y0), 2));
    }
}

double Circle::y_distance_to_circle(double x, double y) const{
    if(this->inCircle(x,y)){
        cerr<<"can't calculate cause the point is in circle"<<endl;
        return -1;
    }
    else{
        return abs(y-y0)-sqrt(pow(radius, w)-pow((x-x0), 2));
    }
}