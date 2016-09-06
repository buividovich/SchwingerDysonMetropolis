# First test which platform we are using
# MinGW+MSYS are assumed for Windows

CC		    :=	gcc
CCFLAGS		:=	-std=gnu99 -Wall -pedantic -D NATIVE_COMPLEX

ifeq ($(OS),Windows_NT)
 CCFLAGS	+=	-fmax-errors=2
endif

ifeq ($(DEBUG), 1)
 CCFLAGS	+=	 -O0 -g3  -Wno-unknown-pragmas -m32
 SUFF       :=   dbg
else
 CCFLAGS	+=	-O2
 SUFF       :=
endif

ifeq ($(PROFILE), 1)
 CCFLAGS	+=	 -pg
 SUFF       :=   prf
endif

#choice of remote hosts - lrz in munich
RHOST_CODE = di72kov@lxlogin5.lrz.de://home/hpc/b3101/di72kov/

############## COMMON DEFINITIONS FOR ALL APPLICATIONS ##############################

CLUE_DIR  = ../clue

#Include and source directories
INC_DIRS := ./include
CLUE_INC_DIRS := $(CLUE_DIR)/include $(CLUE_DIR)/complex
SRC_DIRS = ./src

#Some additional things which should be there for testing
ifeq ($(TESTS), 1)
 INC_DIRS += ./tests
 SRC_DIRS += ./tests
 CCFLAGS     +=  -DTESTS
endif

ifeq ($(OS),Windows_NT)
 CLUE_INC_DIRS	+=	$(CLUE_DIR)/complex/my_complex_sys
endif

INCLUDE	         := $(addprefix -I, $(INC_DIRS))
CLUE_INCLUDE	 := $(addprefix -I, $(CLUE_INC_DIRS))

#List of include files
INCLUDE_FILES := $(wildcard $(addsuffix /*.h, $(INC_DIRS)))
CLUE_INCLUDE_FILES := $(wildcard $(addsuffix /*.h, $(CLUE_INC_DIRS)))

#List of source files
SRC_FILES := $(wildcard $(addsuffix /*.c, $(SRC_DIRS)))

#Manually add some files from CLUE
CLUE_SRC   = clue_logs.c rand_num_generators.c ranlxd.c square_lattice.c clue_io.c
CLUE_SRC_FILES = $(addprefix $(CLUE_DIR)/src/, $(CLUE_SRC))

#Object directory
OBJ_DIR		=	./obj

#List of object files
OBJ_FILES := $(patsubst %.c, $(OBJ_DIR)/$(SUFF)%.o, $(notdir $(SRC_FILES)))
OBJ_FILES += $(patsubst %.c, $(OBJ_DIR)/$(SUFF)%.o, $(notdir $(CLUE_SRC_FILES)))

#List of libraries to be linked with every executable
LIBS        = -lm

#Binary directory
BIN_DIR     =   ./bin

vpath %.c $(SRC_DIRS) $(CLUE_DIR)/src

clean:
	rm -f -v $(OBJ_DIR)/*.o
	rm -f -v $(BIN_DIR)/*$(EXEEXT)

$(OBJ_DIR)/$(SUFF)%.o: 	%.c $(INCLUDE_FILES) $(CLUE_INCLUDE_FILES)
	$(CC) $(CCFLAGS) $(INCLUDE) $(CLUE_INCLUDE) -c $< -o $@

ifeq ($(OS), Windows_NT)
EXEEXT   = .exe
endif

define app_compile_template
 $(1)_DIR  = ./apps/$(1)
 $(1)_SRC  = $$(wildcard $$($(1)_DIR)/*.c)
 $(1)_INC  = $$(wildcard $$($(1)_DIR)/*.h)
$(1): $$(OBJ_FILES) $$(INCLUDE_FILES) $$(CLUE_INCLUDE_FILES) $$($(1)_SRC) $$($(1)_INC)
	$$(CC) $$(CCFLAGS) $$(INCLUDE) $$(CLUE_INCLUDE) -I$$($(1)_DIR) -L$$(CLUE_DIR)/lib $$(OBJ_FILES) $$($(1)_SRC) $(LIBS) $$(LIBFILES) -o $$(BIN_DIR)/$$@$$(SUFF)$$(EXEEXT)
endef

define app_upload_template
 $(1)_DIR  = ./apps/$(1)
 $(1)_SRC  = $$(wildcard $$($(1)_DIR)/*.c)
 $(1)_INC  = $$(wildcard $$($(1)_DIR)/*.h)
$(1)_upload: clue_upload $$(INCLUDE_FILES) $$($(1)_SRC) $$($(1)_INC)
	ssh-pageant rsync -v -R -z ./Makefile $$(SRC_FILES) $$(INCLUDE_FILES) $$($(1)_SRC) $$($(1)_INC) $(RHOST_CODE)/sd_metropolis/
endef

CLUE_INCLUDE_LIST = $(subst ../clue/, ./, $(CLUE_INCLUDE_FILES))

clue_upload: $(CLUE_SRC_FILES) $(CLUE_INCLUDE_FILES)
	cd $(CLUE_DIR); ssh-pageant rsync -v -R -z ./os_profile.pl ./ $(addprefix ./src/, $(CLUE_SRC)) $(CLUE_INCLUDE_LIST) $(RHOST_CODE)/clue/; cd ../sd_metropolis/

APPS          = $(notdir $(shell find ./apps/* -type d))

$(foreach app, $(APPS), $(eval $(call app_compile_template,$(app))))
$(foreach app, $(APPS), $(eval $(call app_upload_template,$(app))))

all: $(APPS)

upload:    $(foreach app, $(APPS), $(app)_upload)
