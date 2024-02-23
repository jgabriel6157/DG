#include "FunctionMapper.hxx"

// Map initialization
std::map<std::string, FunctionMapper::FunctionType1> FunctionMapper::functionMap1;
std::map<std::string, FunctionMapper::FunctionType2> FunctionMapper::functionMap2;

void FunctionMapper::initializeMap() 
{
    functionMap1["legendre"] = SpecialFunctions::legendre;
    functionMap1["legendreDerivative"] = SpecialFunctions::legendreDerivative;
    functionMap1["legendreOrthonormal"] = SpecialFunctions::legendreOrthonormal;
    functionMap1["legendreOrthonormalDerivative"] = SpecialFunctions::legendreOrthonormalDerivative;
    functionMap1["quadratic"] = SpecialFunctions::quadratic;
    functionMap1["quadraticDerivative"] = SpecialFunctions::quadraticDerivative;
    functionMap1["linear"] = SpecialFunctions::linear;
    functionMap1["linearDerivative"] = SpecialFunctions::linearDerivative;

    functionMap2["topHat"] = SpecialFunctions::topHat;
    functionMap2["pulse"] = SpecialFunctions::gaussianPulse;
    functionMap2["sin"] = [](double x) { return sin(x); };
}