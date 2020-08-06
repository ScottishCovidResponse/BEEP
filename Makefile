CXX := mpicxx
CXXFLAGS := -g -O3 -W -Wall -std=c++11
# -B flag forces compilation of all files
BUILD_DIR := ./build
MKDIR_P ?= mkdir -p
SHELL = bash

ifeq (${DEBUG},1)
CXXFLAGS += -D_GLIBCXX_DEBUG # -D_LIBCPP_DEBUG (this is for libc++ / macOS, but fails to link)
endif

ifeq (${COVERAGE},1)
LDFLAGS += --coverage
CXXFLAGS += --coverage
endif

DATA_PIPELINE := 0

ifeq (${DATA_PIPELINE},1)
CXXFLAGS += -DUSE_DATA_PIPELINE

PYTHON           := python3
PYTHON_CONFIG    := $(PYTHON)-config
PYBIND11_INCLUDE := -I$(shell $(PYTHON) -c 'import os; import pybind11; print(os.path.dirname(pybind11.__file__)+"/include")')
PYTHON_CONFIG_FLAGS := $(shell $(PYTHON) -c 'import sys; print("" if (sys.version_info.major < 3 or sys.version_info.minor < 8) else "--embed")')
PYTHON_CFLAGS := $(PYTHON_EXTRA_CFLAGS) $(PYBIND11_INCLUDE) $(shell $(PYTHON_CONFIG) $(PYTHON_CONFIG_FLAGS) --cflags) -UNDEBUG -fPIC
PYTHON_LDFLAGS := $(PYTHON_EXTRA_LDFLAGS) $(shell $(PYTHON_CONFIG) $(PYTHON_CONFIG_FLAGS) --ldflags)

DATA_PIPELINE_DIR := ../data_pipeline_api/bindings/cpp
DATA_PIPELINE_LDFLAGS := -L$(DATA_PIPELINE_DIR)/build -ldatapipeline
DATA_PIPELINE_CPPFLAGS := -I$(DATA_PIPELINE_DIR)

exe_deps := $(DATA_PIPELINE_DIR)/build/libdatapipeline.a

CXXFLAGS := $(CXXFLAGS) $(PYTHON_CFLAGS) $(DATA_PIPELINE_CPPFLAGS)
LDFLAGS := $(LDFLAGS) $(PYTHON_LDFLAGS) $(DATA_PIPELINE_LDFLAGS)

endif

# Please keep these in lexicographic order to aid merging
srcs := \
 chain.cc \
 combinetrace.cc \
 data.cc \
 details.cc \
 generateQ.cc \
 gitversion.cc \
 inputs.cc \
 main.cc \
 mcmc.cc \
 model.cc \
 obsmodel.cc \
 output.cc \
 pack.cc \
 poptree.cc \
 simulate.cc \
 timers.cc \
 utils.cc

objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps := $(objs:.o=.d)

CPPFLAGS := $(CPPFLAGS_EXTRA) -MMD -MP -I$(BUILD_DIR)

# Disable warning for old versions of GCC in files that include toml11/toml.hpp
CXXFLAGS_inputs.cc := -Wno-unused-parameter -O

# Main executable
exe := beepmbp

SRC_DIR := .

# Test executable
TEST_EXEC_NAME := runtests
TEST_NAMES := test_data.cc test_utils.cc
TEST_EXEC := $(BUILD_DIR)/$(TEST_EXEC_NAME)
TEST_EXEC_SRCS := $(SRC_DIR)/$(TEST_EXEC_NAME).cc $(filter-out main.cc,$(srcs)) $(TEST_NAMES:%=$(SRC_DIR)/codetests/%)
TEST_EXEC_OBJS := $(TEST_EXEC_SRCS:%=$(BUILD_DIR)/%.o)

$(exe): $(objs) $(exe_deps)
	$(CXX)  $(objs) $(LDFLAGS) -o $@

$(BUILD_DIR)/%.cc.o: %.cc | gitversion
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_$<) $(CPPFLAGS) -c $< -o $@

# Link the test executable from its corresponding object files
$(TEST_EXEC): $(BUILD_DIR)/% : $(TEST_EXEC_OBJS)
	$(MKDIR_P) $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS)

# $(TARGET_ARCH)

.PHONY : gitversion
gitversion:
	@$(MKDIR_P) $(BUILD_DIR)
	bash ./gitversion.sh $(BUILD_DIR)/gitversion.hh

.PHONY : all
all: $(exe)

.PHONY : clean
clean:
	rm -rf $(exe) $(BUILD_DIR)

-include $(deps)

.PHONY : test
test: $(exe)
	external/regtests/bin/run-all-regression-tests

.PHONY : codetest
codetest: $(TEST_EXEC)
	$(TEST_EXEC)

.PHONY : test-update
test-update:
	set -e; for test in tests/*; do external/regtests/bin/store-regression-test-results regression_test_results/$${test##*/} tests/$${test##*/}; done

.PHONY : test-commit
test-commit:
	set -e; git add tests/*/refdata && git commit -m "Update test reference data"

.PHONY : coverage
coverage :
	gcovr --object-directory $(BUILD_DIR) --exclude toml11 --exclude codetests --exclude catch.hpp --print-summary $(GCOVRFLAGS)
