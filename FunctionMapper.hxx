// #ifndef FUNCTION_MAPPER_HPP
// #define FUNCTION_MAPPER_HPP

// #include <functional>
// #include <string>
// #include <map>
// #include "SpecialFunctions.hxx"

// class FunctionMapper {
// public:
//     // Define the type of function
//     using FunctionType1 = std::function<double(int, double)>;
//     using FunctionType2 = std::function<double(double)>;

//     // Initialize the map with function names and corresponding functions
//     static void initializeMap();

//     // Function to create function based on function name
//     template<typename FunctionType, typename... Args>
//     static FunctionType getFunction(const std::string& functionName);

// private:
//     // Map to store function names and corresponding functions
//     static std::map<std::string, FunctionType1> functionMap1;
//     static std::map<std::string, FunctionType2> functionMap2;

//     // Function to get the appropriate map based on function type
//     template<typename FunctionType>
//     static auto& getFunctionMap();
// };

// #endif

#ifndef FUNCTION_FACTORY_HPP
#define FUNCTION_FACTORY_HPP

#include <functional>
#include <string>
#include <map>

class FunctionFactory {
public:
    // Define the type of function
    using FunctionType1 = std::function<double(int, double)>;
    using FunctionType2 = std::function<double(double)>;

    // Initialize the map with function names and corresponding functions
    static void initializeMaps();

    // Function to create function based on function name
    template<typename FunctionType, typename... Args>
    static FunctionType createFunction(const std::string& functionName);

private:
    // Map to store function names and corresponding functions
    static std::map<std::string, FunctionType1> functionMap1;
    static std::map<std::string, FunctionType2> functionMap2;

    // Function to get the appropriate map based on function type
    template<typename FunctionType>
    static auto& getFunctionMap();
};

#endif

