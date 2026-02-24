#ifndef PARSING_HPP
#define PARSING_HPP
#include <fstream>

#include "simplex.hpp"

namespace prs {
class parser {
private:
    std::ifstream file_;
    smp::task task_;

public:
    explicit parser(std::ifstream file) : file_(std::move(file)) {}

    std::string parse_string() {
        std::string name;
        while(std::isalpha(file_.peek()) || std::isdigit(file_.peek())) {
            name.push_back(file_.get());
        }

        return name;
    }

    void skip_spaces() {
        while(std::isspace(file_.peek()) && file_.peek() != '\n') {
            file_.get();
        }
    }

    int parse_integer() {
        std::string data;

        while(std::isdigit(file_.peek())) {
            data += file_.get();
        }

        return std::stoi(data);
    }

    void to_next_row() {
        while(std::isspace(file_.peek()) || file_.peek() == EOF) {
            file_.get();
        }
    }

    double parse_double() {
        bool negative = false;
        std::string data;

        // Handle sign
        if (file_.peek() == '-') {
            negative = true;
            data += file_.get();
        } else if (file_.peek() == '+') {
            file_.get();
        }

        // Parse integer part
        while (isdigit(file_.peek())) {
            data += file_.get();
        }

        // Parse decimal part
        if (file_.peek() == '.') {
            data += file_.get();
            while (isdigit(file_.peek())) {
                data += file_.get();
            }
        }

        return std::stod(data);
    }

    mat::matrix<double> parse_matrix(int rows, int cols) {
        mat::matrix<double> result(rows, cols);

        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j) {
                skip_spaces();
                result(i, j) = parse_double();
            }
            to_next_row();
        }

        return result;
    }

    mat::vector<double> parse_vector(int rows) {
        mat::vector<double> result(rows);

        for(int i = 0; i < rows; ++i) {
            skip_spaces();
            result(i) = parse_double();
        }

        return result;
    }

    void parse_section() {
        //section starts w/ section name, then params until \n or EOF
        std::string name = parse_string();

        if(name == "A") {
            skip_spaces();
            auto rows = parse_integer();
            skip_spaces();
            auto cols = parse_integer();
            to_next_row();
            auto matrix = parse_matrix(rows, cols);

            task_.A = matrix;
        } else if(name == "b") {
            skip_spaces();
            auto rows = parse_integer();
            to_next_row();
            auto b = parse_vector(rows);
            task_.b = b;
        } else if(name == "c") {
            skip_spaces();
            auto rows = parse_integer();
            to_next_row();
            auto c = parse_vector(rows);
            task_.c = c;
        }
    }

    smp::task parse() {
        while(file_.peek() != EOF) {
            if(std::isspace(file_.peek())) {
                file_.get();
                continue;
            }
            if(std::isalpha(file_.peek())) {
                parse_section();
                continue;
            }
            throw std::runtime_error("parsing error");
        }

        return task_;
    }
};
}

#endif //PARSING_HPP
