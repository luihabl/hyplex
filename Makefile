exe = run
sdr = src
bdr = build
hypre_dir=hypre-2.16.0
cpp_std=c++17

ifeq ($(machine),zoidberg)
  CC = h5c++
  LDFLAGS = -lHYPRE -lm
  CXXFLAGS = -Ilib -I$(HYPREINCLUDE) -std=$(cpp_std) -O3 -Wall
else ifeq ($(machine),hopper)
  CC = mpic++
  LDFLAGS = -lm -L/usr/local/lib -Llib/yaml -lyaml-cpp -L$(HYPREDIR) -lHYPRE -lstdc++fs
  CXXFLAGS = -I/usr/local/include -Ilib -Ilib/yaml -Ilibm -I$(HYPREINCLUDE) -std=$(cpp_std) -O3 -Wall -D FS_EXPERIMENTAL
else
  CC = g++
  LDFLAGS = -L/usr/local/lib -Llib/$(hypre_dir)/src/hypre/lib -Llib/yaml -lmpi -lHYPRE -lyaml-cpp
  CXXFLAGS = -I/usr/local/include -Ilib -Ilib/$(hypre_dir)/src/hypre/include -Ilib/yaml -std=$(cpp_std) -O3 -march=native -Wall
endif

# GIT_HASH=`git rev-parse HEAD`
# COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_VERSION=`git describe --tags --always`
# GIT_BRANCH=`git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_VERSION="\"$(GIT_VERSION)\""

src = $(wildcard $(sdr)/*.cpp)
srn = $(src:src/%=%)
obj = $(srn:%.cpp=$(bdr)/%.o)
dep = $(obj:.o=.d)

$(exe): $(obj)
	$(CC) $^ $(LDFLAGS) $(CXXFLAGS) -o $@

$(bdr)/%.o: $(sdr)/%.cpp | $(bdr)
	$(CC) $(VERSION_FLAGS) $(CXXFLAGS) -c $< -o $@ 

gitversion.c: .git/HEAD .git/index
	echo "const char *gitversion = \"$(shell git rev-parse HEAD)\";" > $@
	
-include $(dep)   #include all dep files in the makefile
$(bdr)/%.d: $(sdr)/%.cpp | $(bdr)
	@$(CPP) $(CXXFLAGS) $< -MM -MT $(@:.d=.o) >$@

$(bdr):
	mkdir -p $@

.PHONY: clean
clean:
	rm -f $(obj) $(exe)
	rm -f $(dep)
