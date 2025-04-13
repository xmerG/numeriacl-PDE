#ifndef _TEST_
#define _TEST_
#include"Multigrid.h"
using namespace std;
void test(const string &filename, const string &output);

unique_ptr<Function> Laplacian(int dim);

unique_ptr<Function> Primitive(int dim);

unique_ptr<Function> BC_Function(int dim, string bc, vector<int> mix=vector<int>{0,0,0,0});

void write_test_results_to_csv(
    int n,
    const std::string& BC,
    const std::string& restriction_opr,
    const std::string& prolongation_opr,
    const std::string& Cycle,
    int dimension,
    const std::chrono::milliseconds& duration
);
#endif