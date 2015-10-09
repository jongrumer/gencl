#
# 	G E N C L
#
# ------------------------------------------------
# 	Modern fortran version of the gencsl(f) 
#
# 	Jon Grumer
# 	jon.grumer@teorfys.lu.se / jongrumer@gmail.com
#
# 	Development version -- use with caution!
# ------------------------------------------------

# ... Compilation setup
FC	       = gfortran
F08FLAGS  = -fcheck=all
F08FLAGS += -fbacktrace
F08FLAGS += -ffpe-trap=zero,overflow,underflow
F08FLAGS += -Wall -Wextra -Wconversion
F08FLAGS += -Wno-unused -Wno-unused-dummy-argument
F08FLAGS += -pedantic
F77FLAGS  = 

# ... Bin directory
BINDIR    = .

# ... Object files
F08_OBJ   = m_parameters.o 	\
				m_globals.o		\
				m_tools.o		\
				m_io.o		\
				m_generate.o	\
				gencl_main.o
F77_OBJ   = gencl77.o

# ... Make all
all: gencl90 gencl77

# ... Make clean (dSYM's are sometime generated in Mac environments)
clean:
	@rm -rf *.o *.mod gencl90.dSYM gencl77.dSYM core

# ... Make cleanall
cleanall:
	@rm -rf *.o *.mod gencl90.dSYM gencl77.dSYM core $(BINDIR)/gencl90 $(BINDIR)/gencl77

# ... Rules

gencl90: $(F08_OBJ)
	$(FC) $(F08FLAGS) $^ -o $@

gencl77: $(F77_OBJ)
	$(FC) $(FFLAGS) $^ -o $@

%.mod: %.o
	$(FC) $(F08FLAGS) -c $<

%.o: %.f08
	$(FC) $(F08FLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $<
