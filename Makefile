CC = g++
CFLAGS = -std=c++17 -Wall

# Source files
SRCS = Vector.cxx Matrix.cxx SpecialFunctions.cxx GaussianQuadrature.cxx FunctionMapper.cxx Mesh.cxx Solver.cxx dg1dAdvection.cxx

# Object files
OBJS = $(SRCS:.cxx=.o)

# Executable name
EXEC = dg1d

# Default target
all: $(EXEC)

# Compile source files
%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files to create the executable
$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC)

# Clean target to remove object files and executable
clean:
	rm -f $(OBJS) $(EXEC)
