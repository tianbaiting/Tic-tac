ROOT := ../..
BUILD_DIR := $(ROOT)/build
TARGET := check_3nf_normalization

CXX ?= g++
CPPFLAGS := -Wall -Wno-sign-compare -std=c++17 -O3 -ggdb -fopenmp
CPPFLAGS += -I$(ROOT)/include
CPPFLAGS += -I$(ROOT)/src -I$(ROOT)/src/interactions -I$(ROOT)/src/utils

LDLIBS := -Wl,--no-as-needed -lgomp -lpthread -lgsl -lgslcblas -lm

TOOL_SRC := $(wildcard *.cpp)
TOOL_OBJ := $(TOOL_SRC:.cpp=.o)

# Object files we need from src/: CMake puts them in CMakeFiles/Tic-tac.dir/
CMAKE_OBJ_DIR := $(BUILD_DIR)/CMakeFiles/Tic-tac.dir/src
SRC_OBJS := $(CMAKE_OBJ_DIR)/utils/spin_isospin_algebra.cpp.o \
            $(CMAKE_OBJ_DIR)/utils/error_management.cpp.o \
            $(CMAKE_OBJ_DIR)/utils/coupling_coefficients.cpp.o \
            $(CMAKE_OBJ_DIR)/interactions/three_nucleon_force_model.cpp.o

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(TOOL_OBJ) $(SRC_OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

$(SRC_OBJS):
	$(MAKE) -C $(ROOT) tic-tac

clean:
	rm -f *.o $(TARGET)
