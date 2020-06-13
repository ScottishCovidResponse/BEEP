MPICXX := mpicxx
CXXFLAGS := -g -O3

srcs :=    output.cc model.cc  MBP.cc obsmodel.cc MBPCHAIN.cc PMCMC.cc analysis.cc data.cc  pack.cc PART.cc  timers.cc utils.cc simulate.cc poptree.cc 
#srcs :=  pack.cc PART.cc PMCMC.cc data.cc output.cc analysis.cc timers.cc utils.cc model.cc simulate.cc poptree.cc 

hdrs := $(wildcard *.h) $(wildcard *.hh)

exe := run

$(exe): $(srcs) $(hdrs)
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(srcs) -o $(exe)

.PHONY : all
all: $(exe)

.PHONY : clean
clean:
	rm -f $(exe)
