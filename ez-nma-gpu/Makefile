INCDIR := ./include
SRCDIR := ./src
OBJDIR := ./obj
BINDIR := ./bin

CUDA_PATH := /usr/local/cuda
INCLUDES  += -I$(CUDA_PATH)/include
LIBS := -L$(CUDA_PATH)/lib64
LDFLAGS := -lrt -lm -lcudart

CX_SOURCE := $(wildcard ${SRCDIR}/*.cpp)
CU_SOURCE := $(wildcard ${SRCDIR}/*.cu)
CX_OBJECT := $(patsubst %.cpp,${OBJDIR}/%.o,$(notdir ${CX_SOURCE}))
CU_OBJECT := $(patsubst %.cu,${OBJDIR}/%.o,$(notdir ${CU_SOURCE}))

TARGET   := ez-nma
BIN_TARG := ${BINDIR}/${TARGET}

CXX  := g++
NVCC := nvcc
LINKER := nvcc -std=c++11 -arch=compute_61 -code=compute_61 -O2 -use_fast_math
CXXFLAGS   := -std=c++11 -I${INCDIR} -O2
NVCCFLAGS  := -std=c++11 -arch=compute_61 -code=compute_61 -O2 -use_fast_math

all:${BIN_TARG}

${BIN_TARG}:${CX_OBJECT} ${CU_OBJECT}
	$(LINKER) -o $@ $^ $(INCLUDES) $(LIBS) $(LDFLAGS)

${OBJDIR}/%.o:${SRCDIR}/%.cpp
	$(CXX) ${CXXFLAGS} $(INCLUDES) -c $< -o $@

${OBJDIR}/%.o:${SRCDIR}/%.cu
	$(NVCC) ${NVCCFLAGS} $(INCLUDES) -c $< -o $@

.PHONY:clean
clean:
	find ${OBJDIR} -name *.o -exec rm -rf {} \;
	rm -rf $(BIN_TARG)
