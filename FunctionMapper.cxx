// #include "FunctionMapper.hxx"
// #include <cmath>
// #include <iostream>

// // Map initialization
// std::map<std::string, FunctionMapper::FunctionType1> FunctionMapper::functionMap1;
// std::map<std::string, FunctionMapper::FunctionType2> FunctionMapper::functionMap2;

// void FunctionMapper::initializeMap() {
//     functionMap1["legendre"] = SpecialFunctions::legendre;
//     functionMap1["legendreDerivative"] = SpecialFunctions::legendreDerivative;
//     functionMap1["legendreOrthonormal"] = SpecialFunctions::legendreOrthonormal;
//     functionMap1["legendreOrthonormalDerivative"] = SpecialFunctions::legendreOrthonormalDerivative;
//     functionMap1["quadratic"] = SpecialFunctions::quadratic;
//     functionMap1["quadraticDerivative"] = SpecialFunctions::quadraticDerivative;
//     functionMap1["linear"] = SpecialFunctions::linear;
//     functionMap1["linearDerivative"] = SpecialFunctions::linearDerivative;

//     functionMap2["topHat"] = SpecialFunctions::topHat;
//     functionMap2["pulse"] = SpecialFunctions::gaussianPulse;
//     functionMap2["sin"] = sin;
// }

// template<typename FunctionType>
// auto& FunctionMapper::getFunctionMap() {
//     if constexpr (std::is_same_v<FunctionType, FunctionType1>) 
//     {
//         return functionMap1;
//     } 
//     else if constexpr (std::is_same_v<FunctionType, FunctionType2>) 
//     {
//         return functionMap2;
//     }
// }

// template<typename FunctionType, typename... Args>
// FunctionType FunctionMapper::getFunction(const std::string& functionName) {
//     auto iterator = getFunctionMap<FunctionType>().find(functionName);
//     if (iterator != getFunctionMap<FunctionType>().end()) 
//     {
//         return iterator->second;
//     } 
//     else 
//     {
//         std::cerr << "Unknown function: " << iterator->First << "\n";
//         return 1;
//     }
// }

#include "FunctionMapper.hxx"

// Map initialization
std::map<std::string, FunctionFactory::FunctionType1> FunctionFactory::functionMap1;
std::map<std::string, FunctionFactory::FunctionType2> FunctionFactory::functionMap2;

void FunctionFactory::initializeMaps() {
    functionMap1["function1"] = [](int x, double y) { return x * y; };
    functionMap1["default1"] = [](int x, double y) { return x + y; };
    functionMap2["function2"] = [](double x) { return x * x; };
    functionMap2["default2"] = [](double x) { return x + 1.0; };
}

template<typename FunctionType>
auto& FunctionFactory::getFunctionMap() {
    if constexpr (std::is_same_v<FunctionType, FunctionType1>) {
        return functionMap1;
    } else if constexpr (std::is_same_v<FunctionType, FunctionType2>) {
        return functionMap2;
    }
}

template<typename FunctionType, typename... Args>
FunctionType FunctionFactory::createFunction(const std::string& functionName) {
    auto it = getFunctionMap<FunctionType>().find(functionName);
    if (it != getFunctionMap<FunctionType>().end()) {
        return it->second;
    } else {
        return getFunctionMap<FunctionType>().at("default");
    }
}

