#include <cstdio>
#include <iostream>

#include "matrix.hpp"
#include "parsing.hpp"
#include "simplex.hpp"

int main() {
    std::ifstream input("input.txt");
    prs::parser parser(std::move(input));

    auto task = parser.parse();

    auto status = smp::simplex_main(task);

    if(status.code == smp::FINISH) {
        //mat::print_matrix(status.x);
        mat::vector<double> res(5);
        res(0) = status.x(0);
        res(1) = status.x(1);
        res(2) = status.x(2);
        res(3) = status.x(3) - status.x(4);
        res(4) = status.x(5) - status.x(6);
        mat::print_matrix(res);
        std::cout << '\n' << -(task.c * status.x) << std::endl;
    }
    if(status.code == smp::UNBOUNDED) {
        std::cout << "SYSTEM IS UNBOUNDED" << std::endl;
    }
    if(status.code == smp::EMPTY) {
        std::cout << "TARGET IS EMPTY" << std::endl;
    }

    return 0;
}
