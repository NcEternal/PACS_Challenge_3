CXX	 = mpic++
LDFLAGS	 = -L../lib -Wl,-R../lib
LINK.o 	 = $(CXX) $(LDFLAGS) -fopenmp
CXXFLAGS = -std=c++20 -O3 -fopenmp -fPIC
CPPFLAGS = -Wall -I../include -I../include/muParser
LDLIBS   = -lstdc++ -lm

ifeq ($(wildcard ../test/cache.dat),)
	REPS = 1
else
	REPS = $(shell cat "../test/cache.dat")
endif

DEPEND = make.dep

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXEC = main performance_test 

main_INCLUDE = main.o Mesh.o Parsers.o Solver.o Input_Output.o Parameter_Communication.o
main: LIBS = -lmuparser
performance_test_INCLUDE = performance_test.o Mesh.o Solver.o  
performance_test.o: CPPFLAGS += -DREPS=$(REPS)


.PHONY = all $(EXEC) $(OBJS) clean distclean $(DEPEND)

all: $(DEPEND) $(EXEC)

.SECONDEXPANSION:
$(EXEC): %: $$(filter $$($$*_INCLUDE),$(OBJS))
	$(LINK.o) $(filter $($@_INCLUDE),$(OBJS)) -o $@ $(LIBS)

$(OBJS): %.o: %.cpp

clean:
	$(RM) $(DEPEND)
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC)
	$(RM) *.csv *.out *.bak *~

$(DEPEND): $(SRCS)
	@$(RM) $(DEPEND)
	@for file in $(SRCS); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

-include $(DEPEND)
