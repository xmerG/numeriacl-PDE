#include<iostream>
#include"EquationSolver.hpp"
#include<cmath>
#include"TestFunctionA.hpp"
#include <nlohmann/json.hpp>
#include"Circle.hpp"
#include<fstream>
#include<string>
#include<memory>
using namespace std;


BoundaryCondition str2BC(const string  &str){
    if(str=="Dirichlet"){return BoundaryCondition::Dirichlet;}
    else if(str=="Neumann"){return BoundaryCondition::Neumann;}
    else{return BoundaryCondition::Mixed;}
}

Domain str2D(const string &str){
    if(str=="regular"){return Domain::regular;}
    else {return Domain::irregular;}
}

unique_ptr<Function> createFunction(const string& bc, const string &d,Circle* c = nullptr, const vector<int> &v=vector<int>{}) {
    if (bc == "Dirichlet" && d=="regular") {
        return make_unique<DirichletF>();
    } else if (bc == "Neumann" && d=="regular") {
        return make_unique<NeumannF>();
    } else if (bc == "Mixed" && d=="regular") {
        return make_unique<Mixed>(v);
    } else if (bc == "Dirichlet" && d=="irregular") {
        return make_unique<primitive>();
    } else if (bc == "Neumann" && d=="irregular") {
        return make_unique<irreNeumann>(c);
    } else {
        return make_unique<irreMixed>(c,v);
    }
}


int main(){
    Laplacian f1;
    primitive f0;
   /* fstream ifs("../input/input.json");
    nlohmann::json json;
    ifs>>json;
    for(const auto &item : json){
        //string functionName=json["function"];
        Circle *c=nullptr;
        vector<int> v;
        string BC=item["boundary_condition"];
        BoundaryCondition bc=str2BC(BC);
        string D=item["domain"];
        Domain d=str2D(D);
        int grid_number=item["grid_number"];
        if(d==Domain::irregular){
            if(item.contains("circle") && item["circle"].is_array()){
                auto data=item["circle"];
                double x=data[0];
                double y=data[1];
                double r=data[2];
                c=new Circle(x, y,r);
            }
        }
        if(bc==BoundaryCondition::Mixed){
            if(item.contains("mixed") && item["mixed"].is_array()){
                vector<int> newv=item["mixed"];
                v=newv;
            }
        }
        cout<<"-------------running test "<<BC<<" "<<D<<" "<<"for grid number "<<grid_number<<"--------------------"<<endl;
        unique_ptr<Function> func = createFunction(BC, D, c, v);
        if (d == Domain::regular && bc == BoundaryCondition::Dirichlet) {
            EquationSolver<Domain::regular, BoundaryCondition::Dirichlet> solver(grid_number, f1);
            solver.solveEquation(*func);
            solver.norm_error(f0,"../output/error1.json");
            solver.print("../output/output1.json", f0);
        } else if (d == Domain::regular && bc == BoundaryCondition::Neumann) {
            EquationSolver<Domain::regular, BoundaryCondition::Neumann> solver(grid_number, f1);
            solver.solveEquation(*func,f0(1.0/grid_number, 1.0/grid_number));
            solver.norm_error(f0,"../output/error1.json");
            solver.print("../output/output1.json", f0);
        } else if (d == Domain::regular && bc == BoundaryCondition::Mixed) {
            EquationSolver<Domain::regular, BoundaryCondition::Mixed> solver(grid_number, f1);
            solver.solveEquation(*func,0.0,v);
            solver.norm_error(f0,"../output/error1.json");
            solver.print("../output/output1.json", f0);
        } else if (d == Domain::irregular && bc == BoundaryCondition::Dirichlet) {
            EquationSolver<Domain::irregular, BoundaryCondition::Dirichlet> solver(grid_number, f1, c);
            solver.solveEquation(*func);
            solver.norm_error(f0,"../output/error1.json");
            solver.print("..output/output1.json", f0);
        } else if (d == Domain::irregular && bc == BoundaryCondition::Neumann) {
            EquationSolver<Domain::irregular, BoundaryCondition::Neumann> solver(grid_number, f1, c);
            solver.solveEquation(*func, f0(1.0/grid_number, 1.0/grid_number));
            solver.norm_error(f0,"../output/error1.json");
            solver.print("../output/output1.json", f0);
        } else{
            EquationSolver<Domain::irregular, BoundaryCondition::Mixed> solver(grid_number, f1, c);
            solver.solveEquation(*func,0.0,v);
            solver.norm_error(f0,"../output/error1.json");
            solver.print("../output/output1.json", f0);
        }
        cout<<"------------------------------test ends----------------------------------"<<endl;
    }

    return 0;
}*/
    //Laplacian f1;
    //DirichletF g;
    //primitive f0;

    /*EquationSolver<Domain::regular, BoundaryCondition::Dirichlet> solver1(32, f1);
    solver1.solveEquation(g);
    solver1.norm_error(f0);
    solver1.print("test.json", f0);*/

    /*vector<int> neumann=vector<int>{1,1,1,1};
    NeumannF g1;
    EquationSolver<Domain::regular, BoundaryCondition::Neumann> solver2(4,f1);
    solver2.solveEquation(g1,f0(1.0/4.0, 1.0/4.0),neumann);
    solver2.norm_error(f0);
    solver2.print("test.json", f0);*/

    /*vector<int> mixed=vector<int>{1,0,0,1};
    Mixed h1;
    EquationSolver<Domain::regular, BoundaryCondition::Mixed> solver5(32, f1);
    solver5.solveEquation(h1, 0.0,mixed);
    solver5.norm_error(f0);
    solver5.print("test.json", f0);*/



    vector<double> D{f0(0.0,0.0), f0(1.0,0.0), f0(0.0,1.0), f0(1.0,1.0)};
    Circle *c = new Circle(0.5, 0.5, 0.2);
    EquationSolver<Domain::irregular,BoundaryCondition::Dirichlet> solver3(8,f1,c);
    solver3.solveEquation(f0);
    solver3.norm_error(f0,"1.json");
    solver3.print("test.json", f0);
    delete c;

    /*irreMixed m1(c,vector<int>{0,0,1,1,1});
    EquationSolver<Domain::irregular, BoundaryCondition::Mixed> solver6(8, f1, c);
    solver6.solveEquation(m1,0.0, vector<int>{0,0,1,1,1});
    solver6.norm_error(f0);
    solver6.print("../output/test1.json", f0);*/

    /*irreNeumann g1(c);
    EquationSolver<Domain::irregular, BoundaryCondition::Neumann> solver4(64,f1,c);
    solver4.solveEquation(g1, f0(1.0/64.0,1.0/64.0));
    solver4.print("test.json", f0);
    solver4.norm_error(f0);*/
    return 0;
}