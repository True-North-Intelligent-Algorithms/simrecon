# @(#)makefile

ifeq ($(shell s=`uname`; echo $${s:0:10}), MINGW32_NT)
    IVE_BASE     = $(HOME)/IVE/win64
else
    IVE_BASE     = $(IVE_EXE)/..
endif

IVE_LIB_DIR  = $(IVE_BASE)/LIB
IVE_INC_DIR  = $(IVE_BASE)/INCLUDE

IWLIBS       =  -limlib

#FFTWTHREAD   = -lfftw3f_threads -lpthread
FFTWLIBS     = -lfftw3f

ifneq ($(shell s=`uname`; echo $${s:0:10}), MINGW32_NT)
  FFTWLIBS += -lfftw3f_threads
endif

IVE_LIB_DIR_FLAG = -L$(IVE_LIB_DIR)

ifeq ($(shell uname),Darwin)
  CC = gcc-mp-4.9
  CXX = gcc-mp-4.9
else
  CC = gcc
endif

CFLAGS      = -I$(IVE_INC_DIR) -Wall -fopenmp -fcx-limited-range
CFLAGS      += -O3 -DNDEBUG
#CFLAGS      += -g
LDFLAGS     = -fopenmp $(IVE_LIB_DIR_FLAG)
ifeq ($(shell uname),Darwin)
    LDFLAGS += -L/opt/local/lib
	CXXFLAGS += -I/opt/local/include
endif

LDLIBS      = $(IWLIBS) $(FFTWTHREAD) $(FFTWLIBS)
ifeq ($(shell uname),Darwin)
    LDLIBS += -framework Accelerate #vecLib deprecated 2015-05-18
    LDFLAGS += -Wl,-rpath,$(IVE_LIB_DIR)
else
	LDLIBS += -llapack -lblas
endif

ifeq ($(shell s=`uname`; echo $${s:0:10}), MINGW32_NT)
  CFLAGS  += -I/usr/local/include
  CXXFLAGS  += -I/usr/local/include
  LDFLAGS += -L/usr/local/lib
  IWLIBS += -live
  LDLIBS += -lgfortran
endif

SRC         = sirecon.c drift_n_phasestep.c small_math_funcs.c helpers.c
#how does one do multiple substitutions? OBJ = $(SRC: .c=.o)

all: sirecon radialft otf2d wiener2d

sirecon: sirecon.o drift_n_phasestep.o small_math_funcs.o helpers.o

radialft: mainfile.o
	$(CC) -o radialft $^ $(LDFLAGS) $(IWLIBS) -lsfftw -lsrfftw

otf2d: radialft2d_Hesper.o tiffhandle.o
	$(CXX) -o otf2d $^ $(LDFLAGS) $(IWLIBS) -lsfftw -lsrfftw -ltiff

wiener2d: wiener2d.o
	$(CC) -o $@ $^ $(LDFLAGS) $(IWLIBS) -lsfftw -lsrfftw
clean:
	$(RM) *.o *.exe *~


sirecon.o: sirecon.h
drift_n_phasestep.o: sirecon.h
small_math_funcs.o: sirecon.h
helpers.o: sirecon.h
# depend:
# 	makedepend -- $(CFLAGS) -- $(SRC)

# DO NOT DELETE

