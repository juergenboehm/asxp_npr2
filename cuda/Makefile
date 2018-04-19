GPUCRYPTO_DIR = ./

include Makefile.in

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

$(OBJS_DIR)/%.o: ppoly.cu pimage.cu
	$(NVCC) $(NVCFLAGS) $(NVINC) -c pimage.cu -o $(OBJS_DIR)/pimage.o
#	$(NVCC) $(NVCFLAGS) $(NVINC) -dc pimage.cu -o $(OBJS_DIR)/pimage.o
#	$(NVCC) $(NVCFLAGS) -dlink $(OBJS_DIR)/ppoly.o $(OBJS_DIR)/pimage.o -o phelp.o -L/usr/local/cuda/lib \
#		-L/usr/lib -lcuda -lcudart -lcudadevrt
		


.PHONY : clean doc

clean:
	rm -f $(TARGET) $(OBJS) $(DEPS)
	rm -f phelp.o

doc: $(SRC_FILES) doxygen.config
	doxygen doxygen.config

ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif