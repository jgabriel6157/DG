#include "Parser.hxx"

// Constructor
Parser::Parser() : x(0.0), vx(0.0), vy(0.0), vz(0.0)
{
    parser.DefineVar("x", &x);
    parser.DefineVar("vx", &vx);
    parser.DefineVar("vy", &vy);
    parser.DefineVar("vz", &vz);
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
        x = x_val;
        vx = vx_val;
        vy = vy_val;
        vz = vz_val;
        return parser.Eval();
    };
}
