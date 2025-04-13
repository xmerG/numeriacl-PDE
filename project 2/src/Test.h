#ifndef _TEST_
#define _TEST_
#include"Multigrid.h"

void test(const string &filename);

unique_ptr<Function> Laplacian(int dim);

unique_ptr<Function> Primitive(int dim);

unique_ptr<Function> BC_Function(int dim, string bc, vector<int> mix=vector<int>{0,0,0,0});
#endif