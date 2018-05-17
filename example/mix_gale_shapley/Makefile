BCXX =		g++
CC =		$(CXX)
<<<<<<< HEAD
CXXFLAGS =      -g -O3 -s -frounding-math
=======
#CXXFLAGS =      -g -O3 -s -frounding-math
>>>>>>> eed1d90c226e1a6035acb9dcd3626dbdc670d08d
CXXFLAGS =       -O3
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
#CXXFLAGS +=     -fopenmp
CXXFLAGS +=	$(shell root-config --cflags)
LDFLAGS =	$(shell root-config --ldflags)
LDLIBS =	$(shell root-config --libs)

CXXFLAGS +=     -I..
LDFLAGS += 	-Wl,-rpath,$(patsubst -L%,%,$(filter \
              -L%,$(shell root-config --libs)))

#Hardcoded attempts:

#LDFLAGS += 	-Wl,-rpath,/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/ROOT/v5-34-30-alice10-8/lib -lssl

#CXXFLAGS += 	/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/AliEn-Runtime/v2-19-le-37/bin/openssl

#LDFLAGS += -Wl,-rpath,/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/ROOT/v5-34-30-alice10-10/lib

#LDFLAGS += 	-Wl,-rpath,/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/ROOT/v5-34-30-alice10-8/lib -lssl

#LDFLAGS += /usr/lib64/libssl.so.0.9.8e

TARGET =	mix_gale_shapley

all:		$(TARGET)

clean:
		/usr/bin/rm -f *~ *.o $(TARGET)

%.cc:		../jet.h
