# Makefile for no_node_theory
# Generates no-node LST classification

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# Eigen path - adjust as needed
# macOS (Homebrew): /opt/homebrew/include/eigen3 or /usr/local/include/eigen3
# Linux: /usr/include/eigen3
EIGEN_PATH ?= /usr/include/eigen3

# If using Homebrew on Apple Silicon
ifeq ($(shell uname -s),Darwin)
    ifeq ($(shell uname -m),arm64)
        EIGEN_PATH = /opt/homebrew/include/eigen3
    else
        EIGEN_PATH = /usr/local/include/eigen3
    endif
endif

INCLUDES = -I$(EIGEN_PATH) -I.

# Source files
SRCS = no_node_theory_v2.cpp Tensor.C
OBJS = $(SRCS:.cpp=.o)
OBJS := $(OBJS:.C=.o)

# Target
TARGET = no_node_theory_v2

# Output files
TEX_OUTPUT = no_node_LSTs.tex
TXT_OUTPUT = no_node_LSTs.txt

.PHONY: all clean run pdf

all: $(TARGET)

$(TARGET): no_node_theory_v2.o Tensor.o
	$(CXX) $(CXXFLAGS) -o $@ $^

no_node_theory_v2.o: no_node_theory_v2.cpp Theory_enhanced.h Tensor.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

Tensor.o: Tensor.C Tensor.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

pdf: run
	pdflatex $(TEX_OUTPUT)
	@echo "Generated: no_node_LSTs.pdf"

clean:
	rm -f $(TARGET) *.o $(TEX_OUTPUT) $(TXT_OUTPUT) *.aux *.log *.pdf

# Debug build
debug: CXXFLAGS += -g -DDEBUG
debug: clean all

# Help
help:
	@echo "Usage:"
	@echo "  make          - Build the program"
	@echo "  make run      - Build and run"
	@echo "  make pdf      - Build, run, and generate PDF"
	@echo "  make clean    - Remove all generated files"
	@echo "  make debug    - Build with debug symbols"
	@echo ""
	@echo "Configuration:"
	@echo "  EIGEN_PATH    - Path to Eigen headers (default: auto-detected)"
	@echo ""
	@echo "Example:"
	@echo "  make EIGEN_PATH=/path/to/eigen3 run"
