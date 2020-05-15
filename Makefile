
CXXFLAGS := -g -O3

srcs := analysis.cc var.cc functions.cc model.cc simulate.cc PMCMC.cc PART.cc init.cc
hdrs := $(wildcard *.h) $(wildcard *.hh)

exe := analysis

$(exe): $(srcs) $(hdrs)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(srcs) -o $(exe)

.PHONY : all
all: $(exe)

.PHONY : clean
clean:
	rm -f $(exe)
