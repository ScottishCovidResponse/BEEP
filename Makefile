CXX := mpicxx
CXXFLAGS := -g -O3 -W -Wall -std=c++11 -Werror=vla
# -B flag forces compilation of all files
BUILD_DIR := ./build
MKDIR_P ?= mkdir -p
SHELL = bash

ifeq (${DEBUG},1)
CXXFLAGS += -D_GLIBCXX_DEBUG # -D_LIBCPP_DEBUG (this is for libc++ / macOS, but fails to link)
endif

srcs := generateQ.cc MBP.cc MBPCHAIN.cc analysis.cc data.cc model.cc obsmodel.cc output.cc pack.cc poptree.cc simulate.cc timers.cc utils.cc

objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps := $(objs:.o=.d)

CPPFLAGS := $(CPPFLAGS_EXTRA) -MMD -MP -I$(BUILD_DIR)

# The TOML parser causes very long compile times, so we compile
# analysis.cc without optimisation
CXXFLAGS_analysis.cc := -O0

exe := run

$(exe): $(objs)
	$(CXX) $(CXXFLAGS) $(objs) -o $@

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
