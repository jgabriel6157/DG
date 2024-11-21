#ifndef FUNCTIONMAPPERHEADERDEF
#define FUNCTIONMAPPERHEADERDEF

#include <functional>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include "SpecialFunctions.hxx"

class FunctionMapper 
{
public:
    // Define the type of function
    using FunctionType1 = std::function<double(int, double)>;
    using FunctionType2 = std::function<double(double)>;

    // Initialize the map with function names and corresponding functions
    static void initializeMap();

    // Template functions must be defined in header file
    // Function to create function based on function name
    template<typename FunctionType, typename... Args>
    static FunctionType getFunction(const std::string& functionName) 
    {
        auto& functionMap = getFunctionMap<FunctionType>();
        auto iterator = functionMap.find(functionName);
        if (iterator != functionMap.end()) 
        {
            return iterator->second;
        } 
        else 
        {
            std::cerr << "Unknown function: " << functionName << "\n";
            // Return a default value based on the function type
            if constexpr (std::is_same_v<FunctionType, FunctionType1>) 
            {
                return SpecialFunctions::linear; // Return a default function for FunctionType1
            } 
            else 
            {
                return [](double x) { return 1.0; }; // Return a default function for FunctionType2
            }
        }
    }

    // Retrieve the derivative of a given function
    static FunctionType1 getDerivative(const std::string& functionName);

    // Function to map from function to its string representation
    static std::string getFunctionName(const FunctionType1& function);

private:
    // Map to store function names and corresponding functions
    static std::map<std::string, FunctionType1> functionMap1;
    static std::map<std::string, FunctionType1> derivativeMap;
    static std::map<std::string, FunctionType2> functionMap2;

    // Template functions must be defined in header file
    template<typename FunctionType>
    static auto& getFunctionMap() 
    {
        if constexpr (std::is_same_v<FunctionType, FunctionType1>) 
        {
            return functionMap1;
        } 
        else if constexpr (std::is_same_v<FunctionType, FunctionType2>) 
        {
            return functionMap2;
        }
    }
};

#endif
