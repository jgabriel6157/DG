#include "Parser.hxx"

// Constructor
Parser::Parser() : x(std::make_shared<double>(0.0)), v(std::make_shared<double>(0.0)) {
    parser.DefineVar("x", x.get());
    parser.DefineVar("v", v.get());
}

// Sets the expression
void Parser::setExpression(const std::string& expression) {
    parser.SetExpr(expression);
}

// Returns a callable function to evaluate the parsed expression
std::function<double(double, double)> Parser::getFunction() {
    return [this](double x_val, double v_val) {
        *x = x_val;
        *v = v_val;
        return parser.Eval();
    };
}
