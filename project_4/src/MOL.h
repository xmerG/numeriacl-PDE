#ifndef _MOL_
#define _MOL_
#include"Triangle_Matrix.h"
#include"Function.h"
#include"Vector.h"
#include <memory>
#include<functional>
using namespace std;

class MOL{
protected:
    double r, t_begin, t_end, T, x_begin, x_end, k, h;
    void set_matrix();
    int N, M; //M 表示空间离散的点数目， N表示时间的
    Triangle_Matrix A=Triangle_Matrix(M);
    vector<Vector> U;
    vector<vector<double>> bc_value;
public:
    MOL(const double &v, const double &h, const double &k, const double &t1, const double &t2, const double &x1, const double &x2);
    void get_bc(const bc_Function &);
    Vector get_initial(const Cauchy_FUnction &);
    virtual void solve(const bc_Function &f, const Cauchy_FUnction &g)=0;
};

class FTCS:public MOL{
public:
    FTCS(const double &v, const double &h, const double &k, const double &t1, 
         const double &t2, const double &x1, const double &x2)
        : MOL(v, h, k, t1, t2, x1, x2) {}
    void solve(const bc_Function &f, const Cauchy_FUnction &g);
};

class Crank_Nicolson:public MOL{
public:
    Crank_Nicolson(const double &v, const double &h, const double &k, const double &t1, 
         const double &t2, const double &x1, const double &x2)
        : MOL(v, h, k, t1, t2, x1, x2) {}
    void solve(const bc_Function &f, const Cauchy_FUnction &g);
};

class BTCS:public MOL{
public:
    BTCS(const double &v, const double &h, const double &k, const double &t1, 
         const double &t2, const double &x1, const double &x2)
        : MOL(v, h, k, t1, t2, x1, x2) {}
    void solve(const bc_Function &f, const Cauchy_FUnction &g);
};


class MOLFactory {
public:
    using CreateMOLCallback = std::function<std::unique_ptr<MOL>(const double&, const double&, 
                                                     const double&, const double&, 
                                                     const double&, const double&, 
                                                     const double&)>;
private:
    using CallbackMap = map<string, CreateMOLCallback>;
    
    CallbackMap callbacks_;
    
    MOLFactory() = default;
    MOLFactory(const MOLFactory&) = delete;
    MOLFactory& operator=(const MOLFactory&) =delete;
    ~MOLFactory()=default;
    
public:
    static MOLFactory& getInstance() {
        static MOLFactory instance;
        return instance;
    }
    
    void registerMOL(const std::string& methodName, CreateMOLCallback createFn) {
        callbacks_[methodName] = createFn;
    }
    
    std::unique_ptr<MOL> createMOL(const std::string& methodName, 
                                  const double& v, const double& h, 
                                  const double& k, const double& t1, 
                                  const double& t2, const double& x1, 
                                  const double& x2) {
        if (!callbacks_.count(methodName)) {
            std::cerr << "MOLFactory: No such method called '" << methodName << "'." << std::endl;
            return nullptr;
        }
        return callbacks_[methodName](v, h, k, t1, t2, x1, x2);
    }
};

void RegisterMOLMethods();

#endif