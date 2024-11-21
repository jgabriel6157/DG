#include "FunctionMapper.hxx"

// Map initialization
std::map<std::string, FunctionMapper::FunctionType1> FunctionMapper::functionMap1;
std::map<std::string, FunctionMapper::FunctionType2> FunctionMapper::functionMap2;
std::map<std::string, FunctionMapper::FunctionType1> FunctionMapper::derivativeMap;

void FunctionMapper::initializeMap() 
{
    functionMap1["legendre"] = SpecialFunctions::legendre;
    functionMap1["legendreOrthonormal"] = SpecialFunctions::legendreOrthonormal;
    functionMap1["quadratic"] = SpecialFunctions::quadratic;
    functionMap1["linear"] = SpecialFunctions::linear;

    derivativeMap["legendre"] = SpecialFunctions::legendreDerivative;
    derivativeMap["legendreOrthonormal"] = SpecialFunctions::legendreOrthonormalDerivative;
    derivativeMap["quadratic"] = SpecialFunctions::quadraticDerivative;
    derivativeMap["linear"] = SpecialFunctions::linearDerivative;

    functionMap2["topHat"] = SpecialFunctions::topHat;
    functionMap2["pulse"] = SpecialFunctions::gaussianPulse;
    functionMap2["twinPulses"] = SpecialFunctions::twinGaussianPulse;
    functionMap2["sin"] = [](double x) { return sin(x); };
}

// Retrieve the derivative of a given function
FunctionMapper::FunctionType1 FunctionMapper::getDerivative(const std::string& functionName) 
{
    auto iterator = derivativeMap.find(functionName);
    if (iterator != derivativeMap.end()) 
    {
        return iterator->second;
    } 
    else 
    {
        std::cerr << "Derivative not found for the provided function.\n";
        return [](int, double) { return 0.0; }; // Default zero derivative
    }
}

std::string FunctionMapper::getFunctionName(const FunctionType1& function) 
{
    for (const auto& entry : functionMap1) 
    {
        if (entry.second.target<FunctionType1::result_type>() == function.target<FunctionType1::result_type>()) 
        {
            return entry.first;
        }
    }
    return "unknown";  // Default if no match found
}