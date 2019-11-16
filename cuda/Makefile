
CUDA_TOOLKIT_DIR = /usr/local/cuda
#CUDA_SDK_DIR = /home/$(USER)/NVIDIA_GPU_Computing_SDK

#generic config
CC = g++
CFLAGS = -O2 -g -Wall

#INC = $(OPENSSL_INCDIR)

#CUDA related config

NVCC = nvcc

NVCFLAGS= -O2 --use_fast_math

NVINC = -I/usr/local/cuda/include -I/usr/local/cuda-10.1/samples/common/inc/

CUDA_LIBDIR = -L$(CUDA_TOOLKIT_DIR)/lib64
CUDA_LIBS = -lcudart

NVCFLAGS += -gencode arch=compute_50,code=sm_50

#NVCFLAGS += -rdc=true --compiler-options '-fPIC'
#NVCFLAGS += --compiler-options '-fPIC'


OBJS_DIR = objs

TARGET_DIR = lib
TARGET_FILE = libasxp.a

TARGET = $(addprefix $(TARGET_DIR)/, $(TARGET_FILE))

.SUFFIXES : .cu .cc .o

CU_SRC_FILES = $(wildcard *.cu)
CC_SRC_FILES = $(wildcard *.cc)
HEADER_FILES = $(wildcard *.hh) $(wildcard *.h)

SRC_FILES = $(CU_SRC_FILES) $(CC_SRC_FILES)
#OBJS_FILE = $(CU_SRC_FILES:.cu=.o) $(CC_SRC_FILES:.cc=.o)
OBJS_FILE = pimage.o

OBJS = $(addprefix $(OBJS_DIR)/, $(OBJS_FILE))
DEPS = Makefile.dep

all: $(TARGET)

$(TARGET): $(DEPS) $(OBJS) | $(TARGET_DIR) $(OBJS_DIR)
	ar rcs $@ $(OBJS)

$(TARGET_DIR):
	mkdir $(TARGET_DIR)

$(OBJS_DIR):
	mkdir $(OBJS_DIR)

$(DEPS): $(SRC_FILES) $(HEADER_FILES)
	$(CC) -MM -MP -x c++ $(CU_SRC_FILES) $(CC_SRC_FILES) | sed 's![^:]*.o:!objs/&!g' > Makefile.dep

$(OBJS_DIR)/%.o : %.cc
	$(CC) $(CFLAGS) $(INC) $(NVINC) -c $< -o $@

$(OBJS_DIR)/%.o: pimage.cu
	$(NVCC) $(NVCFLAGS) $(NVINC) -c pimage.cu -o $(OBJS_DIR)/pimage.o
		


.PHONY : clean doc

clean:
	rm -f $(TARGET) $(OBJS) $(DEPS)


#doc: $(SRC_FILES) doxygen.config
#	doxygen doxygen.config

ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif
