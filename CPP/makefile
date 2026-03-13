ROOT_DIR := ..
ROOT_MAKEFILE := $(ROOT_DIR)/Makefile
TARGET := $(CURDIR)/run
BUILD_DIR := $(CURDIR)/build

.PHONY: all clean cleanall info help

all:
	$(MAKE) -C $(ROOT_DIR) -f $(notdir $(ROOT_MAKEFILE)) all TARGET=$(TARGET) BUILD_DIR=$(BUILD_DIR)

clean:
	$(MAKE) -C $(ROOT_DIR) -f $(notdir $(ROOT_MAKEFILE)) clean TARGET=$(TARGET) BUILD_DIR=$(BUILD_DIR)

cleanall:
	$(MAKE) -C $(ROOT_DIR) -f $(notdir $(ROOT_MAKEFILE)) cleanall TARGET=$(TARGET) BUILD_DIR=$(BUILD_DIR)

info:
	$(MAKE) -C $(ROOT_DIR) -f $(notdir $(ROOT_MAKEFILE)) info TARGET=$(TARGET) BUILD_DIR=$(BUILD_DIR)

help:
	@echo "Compatibility wrapper for CPP/run"
	@echo "This makefile delegates to ../Makefile so CPP/run is compiled from src/include."
