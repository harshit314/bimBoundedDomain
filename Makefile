# generic syntax:
# target: dependency
#(tab)	command

COMPILER = mpif90
FLAGS = -O3 # this is an optimizer which makes loops more efficient
# LINK = -L/home/harshit/packages/lib -llapack -lblas
EXECUTABLE = sccMmatNew.exe
# Object Files for build
OBJS = \
container_mod.o \
rigidBody_mod.o \
solverBoundedGMRES_mod.o \
timeEvolveNew_mod.o \
matrix_mod.o \
global_mod.o

# The commands must be exactly one ‘‘tab’’ over

$(EXECUTABLE) : sccMmatNew.f90 $(OBJS)
	$(COMPILER) $(FLAGS) -o $(EXECUTABLE) sccMmatNew.f90 $(OBJS)

# Generic object dependencies and compilation
%.o : %.f90 
	$(COMPILER) $(FLAGS) -c $<

# $< aliases the first file in the dependency chain

clean: 
	rm -f *.o *.exe *.out 

# Dependency chains:
timeEvolveNew_mod.o : timeEvolveNew_mod.f90 solverBoundedGMRES_mod.o
solverBoundedGMRES_mod.o : solverBoundedGMRES_mod.f90 rigidBody_mod.o
rigidBody_mod.o : rigidBody_mod.f90 container_mod.o
container_mod.o : container_mod.f90 matrix_mod.o
# rigidBodyGMRES_mod.o : rigidBodyGMRES_mod.f90 quadrature_mod.o matrix_mod.o
matrix_mod.o : matrix_mod.f90 global_mod.o
global_mod.o : global_mod.f90
