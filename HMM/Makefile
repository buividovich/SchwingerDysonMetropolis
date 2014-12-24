# First test which platform we are using
# MinGW+MSYS are assumed for Windows

CC		    :=	gcc
CCFLAGS		:=	-D NATIVE_COMPLEX
CCFLAGS	+=	-std=gnu99 -O3 -Wall -pedantic -Wno-unknown-pragmas
#CCFLAGS	    +=	-std=gnu99 -O0 -g -Wall -Wno-unknown-pragmas

############## COMMON DEFINITIONS FOR ALL APPLICATIONS ##############################

CLUE_DIR  = ../../clue

#Include directories
INC_DIRS := ./ $(CLUE_DIR)/include $(CLUE_DIR)/complex $(CLUE_DIR)/lapack

INCLUDE	 := $(addprefix -I, $(INC_DIRS))

#List of include files
INCLUDE_FILES := $(wildcard $(addsuffix /*.h, $(INC_DIRS)))

#Source directories
SRC_DIRS = ./ $(CLUE_DIR)/src

#List of source files
SRC_FILES := $(wildcard $(addsuffix /*.c, $(SRC_DIRS)))

#Object directory
OBJ_DIR		=	./obj

#List of object files
OBJ_FILES := $(patsubst %.c, $(OBJ_DIR)/%.o, $(notdir $(SRC_FILES)))

vpath %.c $(SRC_DIRS)

$(OBJ_DIR)/%.o: 	%.c $(INCLUDE_FILES)
	$(CC) $(CCFLAGS) $(INCLUDE) -c $< -o $@

ifeq ($(OS), Windows_NT)
LIBFILES = $(CLUE_DIR)/lib/liblapacke.lib $(CLUE_DIR)/lib/liblapack.lib $(CLUE_DIR)/lib/libblas.lib
EXEEXT   = .exe
endif

all: hmm

clean:
	rm -f -v $(OBJ_FILES)
	rm -f -v ./*$(EXEEXT)

############# APPLICATION - SPECIFIC DEFINITIONS #################

HMM_SRC = $(wildcard ./*.c)
HMM_INC = $(wildcard ./*.h)

hmm:	$(OBJ_FILES) $(INCLUDE_FILES) $(HMM_SRC) $(HMM_INC)
	$(CC) $(CCFLAGS) $(INCLUDE) -I./ -L$(CLUE_DIR)/lib $(OBJ_FILES) -lm -larpack -llapacke -llapack -lblas -lgfortran -lm $(LIBFILES) -o ./$@$(EXEEXT)
