# Compiler and Flags
CC = gcc            # Define the compiler
CFLAGS = -Wall -g -O2  # Enable warnings, debugging info, and optimize for speed

# Libraries
LIBS = -lm          # Link math library

# Executable and Object Files
TARGET = dynamics   # Name of the final executable
OBJS = src/dynamics.o src/utils.o src/error.o  # List of object files

# Default Target: Build the executable
all: $(TARGET)

# Rule to link object files into the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Compile dynamics.c into dynamics.o
dynamics.o: dynamics.c utils.h error.h
	$(CC) $(CFLAGS) -c src/dynamics.c $(LIBS)

# Compile utils.c into utils.o
utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c src/utils.c $(LIBS)

# Compile error.c into error.o
error.o: error.c error.h
	$(CC) $(CFLAGS) -c src/error.c $(LIBS)

# Clean target: Remove build artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Run target: Build and execute the program
run: all
	./$(TARGET)

