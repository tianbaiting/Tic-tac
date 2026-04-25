// ===============================================================
// 抽取自仓库 [current]: Makefile
// 行号区段：52..75
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
