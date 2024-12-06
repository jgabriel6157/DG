#ifndef PARSER_H
#define PARSER_H

#include "muParser.h"
#include <functional>
#include <memory>
#include <string>

class Parser {
private:
    mu::Parser parser;                // muParser instance
    std::shared_ptr<double> x;        // Pointer to the x variable
    std::shared_ptr<double> v;        // Pointer to the v variable

public:
    // Constructor
    Parser();

    // Sets the expression to be parsed
    void setExpression(const std::string& expression);

    // Returns a callable function to evaluate the parsed expression
    std::function<double(double, double)> getFunction();
};

#endif // PARSER_H
