include Makefile.arch

#------------------------------------------------------------------------------
VPATH    := ../utils
#------------------------------------------------------------------------------

MAINO     = AnaFeedBack.$(ObjSuf)
MAIN      = AnaFeedBack.$(ExeSuf)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    clean

# Build all
all:            $(MAIN)

# Linking step
$(MAIN):        $(MAINO) 
		$(LD) $(MAINO) $(LDFLAGS) $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"
#------------------------------------------------------------------------------


# Explicit rule to compile .cpp files in src/

AnaFeedBack.o: AnaFeedBack.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@
#------------------------------------------------------------------------------

# Clean up
clean:
		@rm -f $(MAIN) *.$(ObjSuf) *.$(DllSuf) core *Dict.* *.*~

.SUFFIXES: .$(SrcSuf)
#------------------------------------------------------------------------------
#.$(SrcSuf).$(ObjSuf):
#	$(CXX) $(CXXFLAGS) -c $< -o $@
#------------------------------------------------------------------------------

