# Adjustable parameters ---------------------------------------------------

MAIN = test/test_cag.cpp
#MAIN = test/test_cag_inf.cpp
#MAIN = test/test_perform_inf.cpp
 # ^ Source file containing main()
LIB_NAME = quad1d
 # ^ Bare name of the library (e.g. SIE for libSIE.a)
OTHER_LIBS =
 # ^ Bare names of other nonstandard static libraries to be included in LIB_NAME.
 #   Each bare library name must match the directory name in which the corresponding source files and Makefile are in.
OBJECTS = cr_qag.o cr_qk.o cr_qk15.o cr_qk21.o cr_qk31.o cr_qk41.o cr_qk51.o cr_qk61.o qag_workspace.o \
		  cr_glfixed.o
 # ^ Library files
C_DIRS = cr_qag/ \
		 cr_glfixed/
 # ^ Directories of source files
H_DIRS = quad1d/ ./
 # ^ Directories of included header files
O_DIR = obj/
 # ^ Directory to put object files
B_DIR = bin/
 # ^ Directory to put binary files
D_DIR = dep/
 # ^ Directory to put dependency files
#L_DIRS = /usr/local/atlas/lib
L_DIRS = /usr/local/lib/
 # ^ Directories of other nonstandard libraries linked to
L_NAMES = gsl gslcblas
 # ^ Bare names of other libraries linked to MAIN (e.g. SIE for libSIE.a)
DOC_DIR = doc/
 # ^ Directory to put documentation files generated by Doxygen

#--------------------------------------------------------------------------

# Less adjustables --------------------------------------------------------

PROGRAM = $(basename $(notdir $(MAIN)))
PROGRAM_O = $(O_DIR)$(PROGRAM).o
LIB_FILE = lib$(LIB_NAME).a
OT_LIBS = $(join $(addsuffix /,$(OTHER_LIBS)),$(addprefix lib,$(addsuffix .a,$(OTHER_LIBS))))
vpath %.c $(C_DIRS)
vpath %.cpp $(C_DIRS)
vpath %.h $(H_DIRS)
OBJS = $(addprefix $(O_DIR),$(OBJECTS))
CPPFLAGS = $(addprefix -I,$(H_DIRS))
LDLIBS = -L. $(addprefix -L,$(L_DIRS)) $(addprefix -l,$(LIB_NAME)) \
         $(addprefix -l,$(L_NAMES))
DEPS = $(addprefix $(D_DIR),$(subst .o,.d,$(OBJECTS)))
CFLAGS = -Ofast -march=native -funroll-loops -pedantic -Wall
LDFLAGS = -march=native
CC = gcc -std=c11
CXX = g++ -std=c++11
CXXFLAGS = -Ofast -march=native -funroll-loops -pedantic -Wall

#--------------------------------------------------------------------------

$(LIB_FILE):: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

.PHONY: program
program: $(PROGRAM_O) | $(B_DIR)
ifeq ($(suffix $(MAIN)),.c)
	$(LINK.o) $< $(LDLIBS) -o $(B_DIR)$(PROGRAM)
else
	$(CXX) $(LDFLAGS) $< $(LDLIBS) -o $(B_DIR)$(PROGRAM)
endif

$(B_DIR):
	mkdir $(B_DIR)

.PHONY: $(PROGRAM_O)
$(PROGRAM_O): $(MAIN) $(LIB_FILE) | $(O_DIR)
ifeq ($(suffix $(MAIN)),.c)
	$(COMPILE.c) $(OUTPUT_OPTION) $<
else
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
endif

$(LIB_FILE):: $(OT_LIBS)
	$(foreach ot_lib,$?,mkdir $(notdir $(basename $(ot_lib))); cd $(notdir $(basename $(ot_lib)))/; ar -x ../$(ot_lib))
	$(foreach ot_lib,$?,$(AR) $(ARFLAGS) $@ $(notdir $(basename $(ot_lib)))/*.o)
	$(foreach ot_lib,$?,rm -rf $(notdir $(basename $(ot_lib))))

$(OBJS): | $(O_DIR)

$(O_DIR)%.o: %.c
	$(COMPILE.c) $(OUTPUT_OPTION) $<

$(O_DIR)%.o: %.cpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

$(O_DIR):
	mkdir $(O_DIR)

$(OT_LIBS):
	cd $(dir $@); make

.PHONY: run
run: program
	$(B_DIR)$(PROGRAM)

.PHONY: valg
valg: program
	valgrind --tool=memcheck --leak-check=yes $(B_DIR)$(PROGRAM)

.PHONY: docs
docs: | $(DOC_DIR)
	doxygen Doxyfile

$(DOC_DIR):
	mkdir $(DOC_DIR)

.PHONY: clean
clean:
	rm -f $(LIB_FILE)
	rm -f $(B_DIR)$(PROGRAM)

.PHONY: cleanall
cleanall: clean
	rm -f $(addprefix $(O_DIR),$(OBJECTS))
	rm -f $(DEPS)
	$(foreach ot_lib,$(OT_LIBS),cd $(dir $(ot_lib)); make cleanall)

.PHONY: build_msg
build_msg:
	@printf "#\n# Building $(notdir $(basename $(MAIN)))\n#\n"

.PHONY: help
help:
	$(MAKE) --print-data-base --question |              \
	awk '/^[^.%][-A-Za-z0-9_]*:/                        \
		{ print substr($$1, 1, length($$1)-1) }' |      \
	sort |                                              \
	pr --omit-pagination --width=80 --columns=4

#--------------------------------------------------------------------------

include $(DEPS)

$(D_DIR)%.d: D_FILE = $(notdir $@)
$(D_DIR)%.d: %.c
	$(CC) -M $(CPPFLAGS) $< > $(D_FILE).$$$$;                   \
	sed 's,\($*\)\.o[ :]*,$(O_DIR)\1.o $(D_DIR)$(D_FILE) : ,g' < $(D_FILE).$$$$ > $(D_FILE);  \
	rm -f $(D_FILE).$$$$
	mv $(D_FILE) $@

$(D_DIR)%.d: D_FILE = $(notdir $@)
$(D_DIR)%.d: %.cpp
	$(CXX) -M $(CPPFLAGS) $< > $(D_FILE).$$$$;                   \
	sed 's,\($*\)\.o[ :]*,$(O_DIR)\1.o $(D_DIR)$(D_FILE) : ,g' < $(D_FILE).$$$$ > $(D_FILE);  \
	rm -f $(D_FILE).$$$$
	mv $(D_FILE) $@

$(DEPS): | $(D_DIR)

$(D_DIR):
	mkdir $(D_DIR)
