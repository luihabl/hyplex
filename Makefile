exe = run
sdr = src
bdr = build
hypre_dir=hypre-2.16.0

ifeq ($(machine),zoidberg)
  CC = g++
  LDFLAGS = -lHYPRE -lmpi -lhdf5
  CXXFLAGS = -Ilib -std=c++11 -O3 -march=native -Wall -D VERBOSE
else
  # CC = mpic++
  # LDFLAGS = -Llib/$(hypre_dir)/src/hypre/lib -lHYPRE -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a
  # CXXFLAGS = -Ilib -Ilib/$(hypre_dir)/src/hypre/include -I/usr/include/hdf5/serial -lpthread -lsz -lz -ldl -std=c++11 -O3 -march=native -Wall -D VERBOSE
  CC = h5c++
  LDFLAGS = -L/usr/local/lib -Llib/$(hypre_dir)/src/hypre/lib -lmpi -lHYPRE
  CXXFLAGS = -I/usr/local/include -Ilib -Ilib/$(hypre_dir)/src/hypre/include -std=c++11 -O3 -march=native -Wall -D VERBOSE
endif

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
