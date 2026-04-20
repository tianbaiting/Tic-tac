ROOT := ../..
BUILD_DIR := $(ROOT)/build
TARGET := check_3nf_normalization

CXX ?= g++
CPPFLAGS := -Wall -Wno-sign-compare -std=c++17 -O3 -ggdb -fopenmp
CPPFLAGS += -I$(ROOT)/include
CPPFLAGS += -I$(ROOT)/src -I$(ROOT)/src/interactions -I$(ROOT)/src/utils

LDLIBS := -Wl,--no-as-needed -lgomp -lpthread -lm

TOOL_SRC := $(wildcard *.cpp)
TOOL_OBJ := $(TOOL_SRC:.cpp=.o)

# Object files we need from src/: populated in Task 4+ as deps emerge.
SRC_OBJS :=

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(TOOL_OBJ) $(SRC_OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET)
