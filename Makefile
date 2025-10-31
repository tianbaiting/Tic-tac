# Tic-tac 三体核物理计算程序 Makefile
# 重构版本 - 2024

# 编译器配置
CXX := g++
FORTRAN := gfortran

# 目标文件名
TARGET := tic-tac

# 目录结构
SRC_DIR := src
INCLUDE_DIR := include
BUILD_DIR := build
DATA_DIR := data

# 源文件目录
CORE_DIR := $(SRC_DIR)/core
UTILS_DIR := $(SRC_DIR)/utils
IO_DIR := $(SRC_DIR)/io
CONFIG_DIR := $(SRC_DIR)/config
INTERACTIONS_DIR := $(SRC_DIR)/interactions

# 查找所有C++源文件
CPP_SOURCES := $(shell find $(SRC_DIR) -name "*.cpp")
FORTRAN_SOURCES := $(shell find $(SRC_DIR) -name "*.f90")
FORTRAN77_SOURCES := $(shell find $(SRC_DIR) -name "*.f")

# 目标文件
CPP_OBJECTS := $(CPP_SOURCES:%.cpp=$(BUILD_DIR)/%.o)
FORTRAN_OBJECTS := $(FORTRAN_SOURCES:%.f90=$(BUILD_DIR)/%.o)
FORTRAN77_OBJECTS := $(FORTRAN77_SOURCES:%.f=$(BUILD_DIR)/%.o)

ALL_OBJECTS := $(CPP_OBJECTS) $(FORTRAN_OBJECTS) $(FORTRAN77_OBJECTS)

# 依赖文件
DEP_FILES := $(CPP_SOURCES:%.cpp=$(BUILD_DIR)/%.d)

# 编译选项
CPPFLAGS := -Wall -Wno-sign-compare -Wno-unused-variable -std=c++17 -O3 -ggdb -fopenmp
CPPFLAGS += -I$(INCLUDE_DIR)
CPPFLAGS += -I$(SRC_DIR) -I$(CORE_DIR) -I$(UTILS_DIR) -I$(IO_DIR) -I$(CONFIG_DIR) -I$(INTERACTIONS_DIR)
CPPFLAGS += -I$(CORE_DIR)/faddeev_solver -I$(CORE_DIR)/state_space -I$(CORE_DIR)/potential -I$(CORE_DIR)/resolvent

# Fortran编译选项
FORTFLAGS_90 := -O3 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC
FORTFLAGS_77 := -O3

# 链接库
LDLIBS := -Wl,--no-as-needed -lgomp -lgsl -lpthread -lm -ldl -lgfortran 
LDLIBS += -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_serial_hl -lhdf5_serial 
LDLIBS += -lstdc++fs -llapacke -llapack -lblas

# 包含依赖文件
-include $(DEP_FILES)

# 默认目标
all: $(TARGET)

# 主要目标：链接所有目标文件
$(TARGET): $(ALL_OBJECTS)
	@echo "链接目标文件生成可执行文件..."
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDLIBS)
	@echo "构建完成: $(TARGET)"

# C++源文件编译规则
$(BUILD_DIR)/%.o: %.cpp
	@echo "编译 C++ 文件: $<"
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) -MMD -MP -c $< -o $@

# Fortran 90源文件编译规则
$(BUILD_DIR)/%.o: %.f90
	@echo "编译 Fortran 90 文件: $<"
	@mkdir -p $(dir $@)
	$(FORTRAN) $(FORTFLAGS_90) -c $< -o $@

# Fortran 77源文件编译规则
$(BUILD_DIR)/%.o: %.f
	@echo "编译 Fortran 77 文件: $<"
	@mkdir -p $(dir $@)
	$(FORTRAN) $(FORTFLAGS_77) -c $< -o $@

# 清理目标文件
clean:
	@echo "清理目标文件..."
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET)

# 清理所有生成文件
cleanall: clean
	@echo "清理所有生成文件..."
	rm -f *.mod *.exe

# 运行程序
run: $(TARGET)
	@echo "运行程序..."
	./$(TARGET) --input data/input.txt

# 调试运行
debug: $(TARGET)
	@echo "调试模式运行..."
	gdb ./$(TARGET)

# 安装目标
install: $(TARGET)
	@echo "安装程序到系统..."
	cp $(TARGET) /usr/local/bin/

# 显示帮助信息
help:
	@echo "Tic-tac 构建系统"
	@echo "可用目标:"
	@echo "  all      - 构建程序 (默认)"
	@echo "  clean    - 清理目标文件"
	@echo "  cleanall - 清理所有生成文件"
	@echo "  run      - 构建并运行程序"
	@echo "  debug    - 调试模式运行"
	@echo "  install  - 安装到系统"
	@echo "  help     - 显示此帮助信息"

# 声明伪目标
.PHONY: all clean cleanall run debug install help

# 显示配置信息
info:
	@echo "=== Tic-tac 构建配置 ==="
	@echo "编译器: $(CXX)"
	@echo "Fortran编译器: $(FORTRAN)"
	@echo "目标文件: $(TARGET)"
	@echo "源文件目录: $(SRC_DIR)"
	@echo "包含目录: $(INCLUDE_DIR)"
	@echo "构建目录: $(BUILD_DIR)"
	@echo "C++编译选项: $(CPPFLAGS)"
	@echo "链接库: $(LDLIBS)"
	@echo "========================"