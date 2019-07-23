CC = mpic++
exe = run
sdr = src
bdr = build
LDFLAGS = -Llib/hypre-2.16.0/src/hypre/lib -lHYPRE
CXXFLAGS = -Ilib -Ilib/hypre-2.16.0/src/hypre/include -std=c++11 -O3 -march=native -Wall -D VERBOSE

src = $(wildcard $(sdr)/*.cpp)
srn = $(src:src/%=%)
obj = $(srn:%.cpp=$(bdr)/%.o)
dep = $(obj:.o=.d)

$(exe): $(obj)
	$(CC) $^ $(LDFLAGS) $(CXXFLAGS) -o $@

$(bdr)/%.o: $(sdr)/%.cpp | $(bdr)
	$(CC) $(CXXFLAGS) -c $< -o $@ 
	
-include $(dep)   #include all dep files in the makefile
$(bdr)/%.d: $(sdr)/%.cpp | $(bdr)
	@$(CPP) $(CXXFLAGS) $< -MM -MT $(@:.d=.o) >$@

$(bdr):
	mkdir -p $@

.PHONY: clean
clean:
	rm -f $(obj) $(exe)
	rm -f $(dep)
