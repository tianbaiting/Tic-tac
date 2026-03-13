# Tic-tac 三体核物理计算程序 Makefile
# 单一源码入口：默认使用 src/include 构建，可通过变量覆盖输出路径。

# 编译器配置
CXX ?= g++
FORTRAN ?= gfortran

# 目标与目录
TARGET ?= tic-tac
SRC_DIR ?= src
INCLUDE_DIR ?= include
BUILD_DIR ?= build
DATA_DIR ?= data
MODULE_DIR := $(BUILD_DIR)/modules

# 源文件目录
CORE_DIR := $(SRC_DIR)/core
UTILS_DIR := $(SRC_DIR)/utils
IO_DIR := $(SRC_DIR)/io
CONFIG_DIR := $(SRC_DIR)/config
INTERACTIONS_DIR := $(SRC_DIR)/interactions

# 查找所有C++/Fortran源文件
CPP_SOURCES := $(shell find $(SRC_DIR) -name "*.cpp" | sort)
FORTRAN_SOURCES := $(shell find $(SRC_DIR) -name "*.f90" | sort)
FORTRAN77_SOURCES := $(shell find $(SRC_DIR) -name "*.f" | sort)

# 目标文件
CPP_OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(CPP_SOURCES))
FORTRAN_OBJECTS := $(patsubst %.f90,$(BUILD_DIR)/%.o,$(FORTRAN_SOURCES))
FORTRAN77_OBJECTS := $(patsubst %.f,$(BUILD_DIR)/%.o,$(FORTRAN77_SOURCES))
ALL_OBJECTS := $(CPP_OBJECTS) $(FORTRAN_OBJECTS) $(FORTRAN77_OBJECTS)

# 依赖文件
DEP_FILES := $(CPP_OBJECTS:.o=.d)

# 显式Fortran模块依赖，避免并行构建时 `use idaho_chiral_potential` 随机失败
CHP_MODULE_OBJECT := $(BUILD_DIR)/$(INTERACTIONS_DIR)/chp/chiral-twobody-potentials.o
CHP_PRESET_OBJECT := $(BUILD_DIR)/$(INTERACTIONS_DIR)/chp/chp-set.o

# 编译选项
CPPFLAGS := -Wall -Wno-sign-compare -Wno-unused-variable -std=c++17 -O3 -ggdb -fopenmp
CPPFLAGS += -I$(INCLUDE_DIR)
CPPFLAGS += -I$(SRC_DIR) -I$(CORE_DIR) -I$(UTILS_DIR) -I$(IO_DIR) -I$(CONFIG_DIR) -I$(INTERACTIONS_DIR)
CPPFLAGS += -I$(CORE_DIR)/faddeev_solver -I$(CORE_DIR)/state_space -I$(CORE_DIR)/potential -I$(CORE_DIR)/resolvent

# Fortran编译选项
FORTFLAGS_90 := -O3 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC
FORTFLAGS_90 += -J$(MODULE_DIR) -I$(MODULE_DIR)
FORTFLAGS_77 := -O3 -I$(MODULE_DIR)

# Optional conda toolchain integration:
# If CONDA_PREFIX is set, prefer headers/libs from that environment.
PKG_CONFIG ?= pkg-config
HDF5_HAVE_SERIAL_PC := $(shell $(PKG_CONFIG) --exists hdf5-serial && echo yes 2>/dev/null)
CONDA_INCLUDE :=
CONDA_LDFLAGS :=
HDF5_LIBS :=
ifneq ($(strip $(CONDA_PREFIX)),)
CONDA_INCLUDE := -I$(CONDA_PREFIX)/include
CONDA_LDFLAGS := -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
HDF5_LIBS := -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
else ifneq ($(strip $(HDF5_HAVE_SERIAL_PC)),)
HDF5_LIBS := -lhdf5_serial_hl_cpp -lhdf5_serial_cpp -lhdf5_serial_hl -lhdf5_serial
else
HDF5_LIBS := -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
endif

CPPFLAGS += $(CONDA_INCLUDE)
LDFLAGS += $(CONDA_LDFLAGS)

# 链接库
LDLIBS := -Wl,--no-as-needed -lgomp -lgsl -lpthread -lm -ldl -lgfortran
LDLIBS += $(HDF5_LIBS)
LDLIBS += -lstdc++fs -llapacke -llapack -lblas -lcurl

# 包含依赖文件
-include $(DEP_FILES)

# 默认目标
all: $(TARGET)

# 链接所有目标文件
$(TARGET): $(ALL_OBJECTS)
	@echo "链接目标文件生成可执行文件..."
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)
	@echo "构建完成: $(TARGET)"

# 显式Fortran模块编译顺序
$(CHP_PRESET_OBJECT): $(CHP_MODULE_OBJECT)

# C++源文件编译规则
$(BUILD_DIR)/%.o: %.cpp
	@echo "编译 C++ 文件: $<"
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) -MMD -MP -c $< -o $@

# Fortran 90源文件编译规则
$(BUILD_DIR)/%.o: %.f90
	@echo "编译 Fortran 90 文件: $<"
	@mkdir -p $(dir $@) $(MODULE_DIR)
	$(FORTRAN) $(FORTFLAGS_90) -c $< -o $@

# Fortran 77源文件编译规则
$(BUILD_DIR)/%.o: %.f
	@echo "编译 Fortran 77 文件: $<"
	@mkdir -p $(dir $@) $(MODULE_DIR)
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
	./$(TARGET) $(DATA_DIR)/input.txt

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
.PHONY: all clean cleanall run debug install help info

# 显示配置信息
info:
	@echo "=== Tic-tac 构建配置 ==="
	@echo "编译器: $(CXX)"
	@echo "Fortran编译器: $(FORTRAN)"
	@echo "目标文件: $(TARGET)"
	@echo "源文件目录: $(SRC_DIR)"
	@echo "包含目录: $(INCLUDE_DIR)"
	@echo "构建目录: $(BUILD_DIR)"
	@echo "模块目录: $(MODULE_DIR)"
	@echo "C++编译选项: $(CPPFLAGS)"
	@echo "Fortran 90 编译选项: $(FORTFLAGS_90)"
	@echo "链接库: $(LDLIBS)"
	@echo "========================"
