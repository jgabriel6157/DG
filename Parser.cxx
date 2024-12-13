#include "Parser.hxx"

// Constructor
Parser::Parser() : x(std::make_shared<double>(0.0)), vx(std::make_shared<double>(0.0)), vy(std::make_shared<double>(0.0)), vz(std::make_shared<double>(0.0)) 
{
    parser.DefineVar("x", x.get());
    parser.DefineVar("vx", vx.get());
    parser.DefineVar("vy", vy.get());
    parser.DefineVar("vz", vz.get());
}

// Sets the expression
void Parser::setExpression(const std::string& expression) 
{
    parser.SetExpr(expression);
}

// Returns a callable function to evaluate the parsed expression
std::function<double(double, double, double, double)> Parser::getFunction() 
{
    return [this](double x_val, double vx_val, double vy_val, double vz_val) 
    {
        *x = x_val;
        *vx = vx_val;
        *vy = vy_val;
        *vz = vz_val;
        return parser.Eval();
    };
}
