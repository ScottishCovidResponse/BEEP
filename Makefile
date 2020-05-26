MPICXX := mpicxx
CXXFLAGS := -g -O3

srcs := PMCMC.cc data.cc PART.cc output.cc analysis.cc timers.cc utils.cc model.cc simulate.cc poptree.cc 
hdrs := $(wildcard *.h) $(wildcard *.hh)

exe := analysis

$(exe): $(srcs) $(hdrs)
	$(MPICXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(srcs) -o $(exe)

.PHONY : all
all: $(exe)

.PHONY : clean
clean:
	rm -f $(exe)
