#pragma once
#include <iostream>
#include <string>
#include <map>

class EquationSolverFactory {
public:
    using CreateEquationSolverCallback = EquationSolver* (*)(int);
private:
    using CallbackMap = std::map<std::string, CreateEquationSolverCallback>;
public:
    void registerEquationSolver(const std::string& ID, CreateEquationSolverCallback createFn) {
        callbacks_[ID] = createFn;
    }
    EquationSolver* createEquationSolver(const std::string &ID) {
        if(!callbacks_.count(ID)){
            std::cerr << "EquationSolver:: No such Equation Solver called '" << ID << "'." << std::endl;
            return nullptr;
        }
        return callbacks_[ID];
    }
private:
    CallbackMap callbacks_;
private:
    EquationSolverFactory() = default;
    EquationSolverFactory(const EquationSolverFactory&) = default;
    EquationSolverFactory& operator = (const EquationSolverFactory&) = default;
    ~EquationSolverFactory() = default;
public:
    static EquationSolverFactory& getInstance() {
        static EquationSolverFactory instance;
        return instance;
    }
};

void RegisterAllEquationSolvers() {
    auto& fac = EquationSolverFactory::getInstance();
    fac.createEquationSolver("Bisection", BisectionSolver());
    fac.createEquationSolver("Newton", NewtonSolver());
    fac.createEquationSolver("Secant", SecantSolver());
}