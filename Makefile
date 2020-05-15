
CXXFLAGS := -g -Wall

srcs := analysis.cc
hdrs := $(wildcard *.h)

exe := analysis

$(exe): $(srcs) $(hdrs)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(srcs) -o $(exe)

.PHONY : all
all: $(exe)

.PHONY : clean
clean:
	rm -f $(exe)
