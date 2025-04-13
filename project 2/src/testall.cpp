#include"Test.h"

int main(){
    test("../input/linear+full_weighting.json", "../output/linear_fullweighting.json");
    test("../input/linear+injection.json", "../output/linear_injection.json");
    test("../input/quadratic+full_weighting.json", "../output/quadratic_full_weighting.json");
    test("../input/quadratic+injection.json", "../output/quadratic_injection.json");
    return 0;
}