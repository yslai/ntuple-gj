CXX =		g++
CC =		$(CXX)
#CXXFLAGS =      -g -O -s -frounding-math
CXXFLAGS =      -O
CXXFLAGS +=     -Wextra -Wall -Werror -fPIC
CXXFLAGS +=	-Wno-unused-variable -Wno-unused-parameter
CXXFLAGS +=	-Wno-unused-function
ifeq ($(findstring g++,$(CXX)),g++)
    CXXFLAGS += -no-canonical-prefixes -Wno-builtin-macro-redefined \
                -D__DATE__="redacted" -D__TIMESTAMP__="redacted" \
                -D__TIME__="redacted" -U_FORTIFY_SOURCE \
                -D_FORTIFY_SOURCE=1 -fstack-protector
endif
CXXFLAGS +=     -std=c++11
CXXFLAGS +=	$(shell root-config --cflags)
LDFLAGS =	$(shell root-config --ldflags)
LDLIBS =	$(shell root-config --libs)

#HDF5_ROOT =	/opt/local/
HDF5_ROOT = 	/usr/common/usg/software/hdf5/1.8.13/
CXXFLAGS +=	-I$(HDF5_ROOT)/include
LDLIBS +=	 -L$(HDF5_ROOT)/lib -lhdf5_cpp -lhdf5
CXXFLAGS +=	-I$(BOOST_ROOT)/include

#LIDLIBS for other systems
#LDLIBS +=	-L$(HDF5_ROOT)/lib64 -L$(HDF5_ROOT)/lib -lhdf5_cpp -lhdf5


TARGET =	

all:		$(TARGET)

clean:
		/usr/bin/rm -f *~ *.o $(TARGET)
