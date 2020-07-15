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

srcs := generateQ.cc MBP.cc MBPCHAIN.cc analysis.cc data.cc model.cc obsmodel.cc output.cc pack.cc poptree.cc simulate.cc timers.cc utils.cc

objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps := $(objs:.o=.d)

CPPFLAGS := $(CPPFLAGS_EXTRA) -MMD -MP -I$(BUILD_DIR)

# Disable warning for old versions of GCC in files that include toml11/toml.hpp
CXXFLAGS_analysis.cc := -Wno-unused-parameter
CXXFLAGS_model.cc    := -Wno-unused-parameter

exe := run

$(exe): $(objs)
	$(CXX)  $(objs) $(LDFLAGS) -o $@

$(BUILD_DIR)/%.cc.o: %.cc | gitversion
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_$<) $(CPPFLAGS) -c $< -o $@

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

.PHONY : test-update
test-update:
	set -e; for test in tests/*; do external/regtests/bin/store-regression-test-results regression_test_results/$${test##*/} tests/$${test##*/}; done

.PHONY : coverage
coverage :
	gcovr --object-directory $(BUILD_DIR) --exclude toml11 --print-summary
