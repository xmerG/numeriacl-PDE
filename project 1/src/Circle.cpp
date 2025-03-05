#include"Circle.hpp"
using namespace std;

Circle::Circle(){}

Circle::Circle(double x, double y):x0(x), y0(y){}

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