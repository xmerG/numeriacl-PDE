#ifndef FUNCTION
#define FUNCTION

class Function {
public:
    virtual double operator() (const double &x, const double &y) const = 0;
    virtual double operator() (const double &x) const=0;
};

#endif
