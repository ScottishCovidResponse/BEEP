
CXXFLAGS := -g -O3

srcs := analysis.cc timers.cc utils.cc model.cc simulate.cc PMCMC.cc PART.cc poptree.cc
hdrs := PART.hh PMCMC.hh consts.hh model.hh poptree.hh simulate.hh timers.hh utils.hh
hdrs += gitversion.hh

exe := analysis

$(exe): $(srcs) $(hdrs) | gitversion
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(srcs) -o $(exe)

.PHONY : gitversion
gitversion:
	./gitversion.sh gitversion.hh

gitversion.hh : gitversion

.PHONY : all
all: $(exe)

.PHONY : clean
clean:
	rm -f $(exe) gitversion.hh
