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
    using FunctionType3 = std::function<double(double, double)>;

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
            else if constexpr (std::is_same_v<FunctionType, FunctionType2>)
            {
                return [](double x) { return 1.0; }; // Return a default function for FunctionType2
            }
            else
            {
                return [](double x, double vx) { return 1.0; }; // Return a default function for FunctionType3
            }
        }
    }

private:
    // Map to store function names and corresponding functions
    static std::map<std::string, FunctionType1> functionMap1;
    static std::map<std::string, FunctionType2> functionMap2;
    static std::map<std::string, FunctionType3> functionMap3;

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
        else if constexpr (std::is_same_v<FunctionType, FunctionType3>) 
        {
            return functionMap3;
        }
    }
};

#endif
