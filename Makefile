# Makefile for project in ./src directory
# Compiles everything into ./build
# Includes depend on make options
# --------------------------------------------------------------------------------------------------
# make all: Compiles everything
# make clean: Cleans up build files
# make run: Compiles everything and runs executable
# make test: Compiles everything into test environment and runs test executable
# --------------------------------------------------------------------------------------------------

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -pedantic -g -fopenmp -Iinclude -MMD -MP
LDFLAGS = -fopenmp -lm
SRC_DIR = src
OBJ_DIR = build
BIN = moca

SRCS = $(wildcard $(SRC_DIR)/*.c)
SRC_WITHOUT_MAIN := $(filter-out $(SRC_DIR)/main.c, $(wildcard $(SRC_DIR)/*.c))
OBJS = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRCS))
DEPS = $(OBJS:.o=.d)

# Make all compiles everything in ./src and dependencies from ./include
all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Monte Carlo application compiled OK"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

-include $(DEPS)

# Make clean is make clean
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d $(BIN)
	@echo "Cleaned up"

# Make run compiles everything and runs
run: all
	./moca -omp 8 input

# Testing environment. Compiles everything in ./test into one binary and runs
# This replaces the main() with the one in tests.c

TEST_SRC := $(wildcard test/*.c) $(SRC_WITHOUT_MAIN)
TEST_OBJ := $(patsubst test/%.c, build/%.o, $(TEST_SRC))
TEST_BIN := build/tests

test: all $(TEST_BIN)
	@./$(TEST_BIN)

$(TEST_BIN): $(TEST_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

build/%.o: test/%.c
	$(CC) $(CFLAGS) -Iinclude -c $< -o $@
