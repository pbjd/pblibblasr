CXXOPTS := -I . -pedantic -c -fmessage-length=0 -MMD -MP

sources := $(wildcard algorithms/alignment/*.cpp) \
		   $(wildcard algorithms/alignment/sdp/*.cpp) \
		   $(wildcard datastructures/alignment/*.cpp) \
		   $(wildcard datastructures/matrix/*.cpp) \
		   $(wildcard datastructures/reads/*.cpp) \
		   $(wildcard datastructures/anchoring/*.cpp) \
		   $(wildcard qvs/*.cpp) \
		   $(wildcard tuples/*.cpp) \
		   $(wildcard utils/*.cpp) \
		   $(wildcard *.cpp) \

objects := $(sources:.cpp=.o)

all : GCCOPTS = -O3

debug : GCCOPTS = -g

profile : GCCOPTS = -Os -pg

all debug profile: libblasr.a

libblasr.a: libblasr.a($(objects))
	$(AR) $(ARFLAGS)s $@ $%

%.o: %.cpp
	$(CXX) $(GCCOPTS) $(CXXOPTS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -o $@ $<

.INTERMEDIATE: $(objects)

clean: 
	rm -f libblasr.a

-include $(sources:.cpp=.d)
