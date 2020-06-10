MPICXX := mpicxx
CXXFLAGS := -g -O3 -std=c++11

srcs := MBP.cc MBPCHAIN.cc PART.cc PMCMC.cc analysis.cc data.cc model.cc obsmodel.cc output.cc pack.cc poptree.cc simulate.cc timers.cc utils.cc
hdrs := MBP.hh MBPCHAIN.hh PART.hh PMCMC.hh consts.hh data.hh model.hh obsmodel.hh output.hh pack.hh poptree.hh simulate.hh timers.hh utils.hh
ifeq (${DEBUG},1)
CXXFLAGS += -D_GLIBCXX_DEBUG # -D_LIBCPP_DEBUG (this is for libc++ / macOS, but fails to link)
endif

hdrs += gitversion.hh

exe := run

$(exe): $(srcs) $(hdrs) | gitversion
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(srcs) -o $(exe)

.PHONY : gitversion
gitversion:
	bash ./gitversion.sh gitversion.hh

gitversion.hh : gitversion

.PHONY : all
all: $(exe)

.PHONY : clean
clean:
	rm -f $(exe) gitversion.hh
