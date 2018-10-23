
include Makefile.in

MYEXEC   =  $(BIN_DIR)/createStructuredTetMesh.exe
MYFCFUN  =  vtkComMod.f90          \
            vtkTypeParams.f90      \
            vtkXMLParser.f90       \
            main.f90

MYCFUN   =  vtkZpipe.c

MYSRC    =  $(patsubst %,    $(SRC_DIR)/%,  $(MYFCFUN))
MYOBJ    =  $(patsubst %.f90,$(OBJ_DIR)/%.o,$(MYFCFUN))
MYOBJ   +=  $(patsubst %.c,  $(OBJ_DIR)/%.o,$(MYCFUN))

LFLAGS   = $(FCFLAGS) $(Z_LFLAGS)
INCLUDES =

#---- Rules ----#

.PHONY: $(MYEXEC)
$(MYEXEC): $(MYOBJ)
	$(FC) $^ $(LFLAGS) -o $@

$(MYOBJ): | $(OBJ_DIR) $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)

cleanall:
	rm -rf $(BIN_DIR) $(OBJ_DIR)
