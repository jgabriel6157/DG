#ifndef PARSER_H
#define PARSER_H

#include "muParser.h"
#include <functional>
#include <memory>
#include <string>

class Parser 
{
private:
    mu::Parser parser;                // muParser instance
    std::shared_ptr<double> x;        // Pointer to the x variable
    std::shared_ptr<double> vx;        // Pointer to the vx variable
    std::shared_ptr<double> vy;        // Pointer to the vy variable
    std::shared_ptr<double> vz;        // Pointer to the vz variable

public:
    // Constructor
    Parser();

    // Sets the expression to be parsed
    void setExpression(const std::string& expression);

    // Returns a callable function to evaluate the parsed expression
    std::function<double(double, double, double, double)> getFunction();
};

#endif // PARSER_H
