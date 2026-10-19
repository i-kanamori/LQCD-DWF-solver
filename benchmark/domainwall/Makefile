.SUFFIXES: .cpp .h .o

BRIDGE_INSTALL_PATH = ../../build

include $(BRIDGE_INSTALL_PATH)/make.inc

CXXFLAGS += -I./

sources = main.cpp \
 job_Utils.cpp \
 test.cpp \
 spectrum_alt.cpp \
 test_alt_QXS.cpp \
 test_alt_Corelib.cpp \
 test_alt_Accel.cpp \
 test_alt_SIMD.cpp \
 spectrum_Domainwall_alt_Corelib.cpp \

objects = $(sources:%.cpp=%.o)

program = bridge.elf

# relink when the libraries are rebuilt, otherwise the executable
# silently keeps the previously linked library code.
libraries = $(BRIDGE_INSTALL_PATH)/libbridge.a \
 $(BRIDGE_INSTALL_PATH)/libbridgecuda.a

# object file to show version (always compiled)
src_show_version =show_version.cpp
obj_show_version =show_version.o

## version info
GIT_COMMIT_AUTHOR := $(shell git log -1 --format='Author: %an <%ae>')
GIT_COMMIT_ID := $(shell git log -1 --format='commit %H')
GIT_COMMIT_DATE := $(shell git log -1 --format='%ad')
#BRIDGE_VERSION := $(shell cat ../../VERSION)

VERSION_INFO= -D_COMMIT_AUTHOR="$(GIT_COMMIT_AUTHOR)"  -D_COMMIT_ID="$(GIT_COMMIT_ID)"  -D_COMMIT_DATE="$(GIT_COMMIT_DATE)"

all: msg $(program)

$(program): $(objects) $(obj_show_version) $(wildcard $(libraries))
	$(LD) $(LDFLAGS) -o $@ $(objects) $(obj_show_version) $(LDLIBS)

.PHONY: $(src_show_version)
$(obj_show_version): $(src_show_version)
	$(CXX) $(CXXFLAGS) $(VERSION_INFO) -c -o $@ $<

#	$(CXX) $(CXXFLAGS) $(VERSION_INFO) -c -o $@ $<

%.o: %.cpp $(headers)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

msg:
	@cat $(BRIDGE_INSTALL_PATH)/message.txt

clean:
	rm -f $(program) $(objects)

