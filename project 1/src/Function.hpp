#ifndef FUNCTION
#define FUNCTION

class Function {
public:
    virtual double operator() (double x, double y) const = 0;
};

#endif
