include Makefile.arch
#------------------------------------------------------------------------------

VPATH    := ../utils
#------------------------------------------------------------------------------

MAINO     = src/CalcDisc.$(ObjSuf)
MAIN      = CalcDisc.$(ExeSuf)

#MAINO     = src/EventCheck.$(ObjSuf)
#MAIN	  = EventCheck.$(ExeSuf)

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

src/CalcDisc.o: src/CalcDisc.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@
#------------------------------------------------------------------------------

# Clean up
clean:
		@rm -f $(MAIN) src/*.$(ObjSuf) *.$(DllSuf) core *Dict.* *.*~

.SUFFIXES: .$(SrcSuf)
#------------------------------------------------------------------------------


#.$(SrcSuf).$(ObjSuf):
#	$(CXX) $(CXXFLAGS) -c $< -o $@
#------------------------------------------------------------------------------

