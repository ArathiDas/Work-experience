# Compiler and Flags
CC = gcc
CFLAGS = -Wall -g -O2

LIBS = -lm
TREX = -ltrexio

# Executable and Object Files
TARGET = mp2_energy
OBJS   = src/mp2_energy.o src/utils.o

# Default Target
all: $(TARGET)

# Rule to Build the Executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(TREX) $(LIBS)

# Compile Source Files into Object Files
mp2_energy.o: mp2_energy.c
	$(CC) $(CFLAGS) -c mp2_energy.c

# Clean Target to Remove Build Artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Run the Executable
run: all
	./$(TARGET)

