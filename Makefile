# First test which platform we are using
# MinGW+MSYS are assumed for Windows

CC		    :=	gcc
CCFLAGS		:=	-D NATIVE_COMPLEX

ifeq ($(DEBUG), 1)
 CCFLAGS	+=	-std=gnu99 -O0 -g3 -Wall -pedantic -Wno-unknown-pragmas
else
 CCFLAGS	+=	-std=gnu99 -O3 -Wall -pedantic -Wno-unknown-pragmas
endif

############## COMMON DEFINITIONS FOR ALL APPLICATIONS ##############################

CLUE_DIR  = ../clue

#Include directories
INC_DIRS := ./include $(CLUE_DIR)/include $(CLUE_DIR)/complex $(CLUE_DIR)/lapack

INCLUDE	 := $(addprefix -I, $(INC_DIRS))

#List of include files
INCLUDE_FILES := $(wildcard $(addsuffix /*.h, $(INC_DIRS)))

#Source directories
SRC_DIRS = ./src $(CLUE_DIR)/src

#List of source files
SRC_FILES := $(wildcard $(addsuffix /*.c, $(SRC_DIRS)))

#Object directory
OBJ_DIR		=	./obj

#List of object files
OBJ_FILES := $(patsubst %.c, $(OBJ_DIR)/%.o, $(notdir $(SRC_FILES)))

#List of libraries to be linked with every executable
LIBS        = -lm -larpack -llapacke -llapack -lblas -lgfortran -lm

#Binary directory
BIN_DIR     =   ./bin

vpath %.c $(SRC_DIRS)

clean:
	rm -f -v $(OBJ_FILES)
	rm -f -v $(BIN_DIR)/*$(EXEEXT)

$(OBJ_DIR)/%.o: 	%.c $(INCLUDE_FILES)
	$(CC) $(CCFLAGS) $(INCLUDE) -c $< -o $@

ifeq ($(OS), Windows_NT)
LIBFILES = $(CLUE_DIR)/lib/liblapacke.lib $(CLUE_DIR)/lib/liblapack.lib $(CLUE_DIR)/lib/libblas.lib
EXEEXT   = .exe
endif

define app_compile_template
 $(1)_DIR  = ./apps/$(1)
 $(1)_SRC  = $$(wildcard $$($(1)_DIR)/*.c)
 $(1)_INC  = $$(wildcard $$($(1)_DIR)/*.h)
$(1): $$(OBJ_FILES) $$(INCLUDE_FILES) $$($(1)_SRC) $$($(1)_INC)
	$$(CC) $$(CCFLAGS) $$(INCLUDE) -I$$($(1)_DIR) -L$$(CLUE_DIR)/lib $$(OBJ_FILES) $$($(1)_SRC) $(LIBS) $$(LIBFILES) -o $$(BIN_DIR)/$$@$$(EXEEXT)
endef

define app_plotclean_template
$(1)_plotclean:
	rm -f -v ./plots/$(1)/*
endef

define app_dataclean_template
$(1)_dataclean: $(1)_plotclean
	rm -f -v ./data/$(1)/*
endef

define app_plots_template
 $(1)_PLT_FILES = $$(wildcard ./apps/$(1)/*.plt)
$(1)_plots: 
	$(foreach plt_file, $$($(1)_PLT_FILES), gnuplot $(plt_file))
	gnuplot -e "app_name = '$(1)'" ./plots/mc_stat.plt
endef

APPS          = $(notdir $(shell find ./apps/* -type d))

$(foreach app, $(APPS), $(eval $(call app_compile_template,$(app))))
$(foreach app, $(APPS), $(eval $(call app_plotclean_template,$(app))))
$(foreach app, $(APPS), $(eval $(call app_dataclean_template,$(app))))
$(foreach app, $(APPS), $(eval $(call app_plots_template,$(app))))

all: $(APPS)

dataclean: $(foreach app, $(APPS), $(app)_dataclean)
plotclean: $(foreach app, $(APPS), $(app)_plotclean)
plots:     $(foreach app, $(APPS), $(app)_plots)




