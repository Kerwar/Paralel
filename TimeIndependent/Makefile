#---------------------- Makefile for parallel -------------------------#
# usage: make compile ; make run or make pbs
#----------------------------------------------------------------------#
#============================= set MPI, compiler ======================#
#If your compiler is NOT on your path (for your shell) then
# you need to insert the full path, e.g. /opt/intel/..../bin/ifort
##-----> set appropriate compiler_wrapper: mpif77 mpif90 mpicc mpic++
COMP = mpifort
##-----> set appropriate extension: f c cpp
EXT = f90
LFLAGs =
#for C: LFLAGs = -lm
##-------------------------- for all:
FLAGs =  -g -Wall -fbacktrace -finit-local-zero -ffpe-trap=invalid,zero,overflow  -fbounds-check -fcheck=all $(MPIlnk) #--showme
# FLAGs = -O3 $(MPIlnk)
MPI = /usr
##---------------------> set path to openMPI on local:
MPIlnk = -I$(MPI)/include -L $(MPI)/lib
#
#=========================== set source code =========================#
##--------------->set names for your PROGram and std I/O files:
PROG = code2D
PARAL = paral
OUTPUT = out
##--------------------> set code components:
#CODE_o = type_m.o param_m.o vel_m.o u_m.o v_m.o p_m.o point_m.o solver_m.o data_m.o main.o
CODEP_o = type_m.o param_m.o interfaces_m.o point_m.o point2_m.o coefun_m.o  paralel_m.o paralel.o

#-==========================create executable: make compile ==========#
#------------> lines below a directive MUST start with TAB <-----------#
$(CODE_o):%.o: %.$(EXT)
	$(COMP) $(FLAGs) -c $< -o $@

$(CODEP_o):%.o: %.$(EXT)
	$(COMP) $(FLAGs) -c $< -o $@
compile:$(CODE_o)
# $(COMP) $(FLAGs) $(CODE_o) -o $(PROG).x $(LFLAGs)
	$(COMP) $(FLAGs) $(CODE_o) -o $(PROG) $(LFLAGs)
	@echo " >>> compiled on ‘hostname -s‘ with $(COMP) <<<"

compilep:$(CODEP_o)
# $(COMP) $(FLAGs) $(CODE_o) -o $(PROG).x $(LFLAGs)
	$(COMP) $(FLAGs) $(CODEP_o) -o $(PARAL) $(LFLAGs)
	@echo " >>> compiled on ‘hostname -s‘ with $(COMP) <<<"
#----------------------- execute: make run --------------------------#
run:
# $(MPI)/bin/mpiexec -n 2 ./$(PROG).x < $(INPUT) > $(OUTPUT)
# $(MPI)/bin/mpirun -np 2 ./$(PROG).x < $(INPUT) > $(OUTPUT)
	$(MPI)/bin/mpirun -np 2 ./$(PROG)

runp:
# $(MPI)/bin/mpiexec -n 2 ./$(PROG).x < $(INPUT) > $(OUTPUT)
# $(MPI)/bin/mpirun -np 2 ./$(PROG).x < $(INPUT) > $(OUTPUT)
	mpirun -np 3 ./$(PARAL) 
#----------------------- execute: make pbs --------------------------#
pbs:
	@ vi PBSscript
	make clean
	qsub PBSscript
#-------------------------- clean -----------------------------------#
clean:
	rm -f o.* DONE *.o watch
#----------------------------------------------------------------------#
