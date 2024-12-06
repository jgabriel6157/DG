CC = g++
CFLAGS = -std=c++17 -Wall -fopenmp -I$(MUPARSER_DIR)

# muParser directory and source files
MUPARSER_DIR = muparser-2.3.4/src
MUPARSER_SRCS = $(MUPARSER_DIR)/muParserBase.cpp $(MUPARSER_DIR)/muParserBytecode.cpp \
                $(MUPARSER_DIR)/muParserCallback.cpp $(MUPARSER_DIR)/muParser.cpp \
                $(MUPARSER_DIR)/muParserDLL.cpp $(MUPARSER_DIR)/muParserError.cpp \
                $(MUPARSER_DIR)/muParserInt.cpp $(MUPARSER_DIR)/muParserTokenReader.cpp

MUPARSER_OBJS = $(MUPARSER_SRCS:.cpp=.o)

# Your source files
SRCS = Vector.cxx Matrix.cxx NewtonCotes.cxx NewtonSolver.cxx SpecialFunctions.cxx \
       GaussianQuadrature.cxx FunctionMapper.cxx Parser.cxx Mesh.cxx Solver.cxx \
       dg1dAdvection.cxx

OBJS = $(SRCS:.cxx=.o)

EXEC = dg1d

# Default target
all: $(EXEC)

# Compile source files
%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files into the executable
$(EXEC): $(OBJS) $(MUPARSER_OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(MUPARSER_OBJS) -o $@

# Clean up object files and the executable
clean:
	rm -f $(OBJS) $(MUPARSER_OBJS) $(EXEC)
