EXE_DIR = .
SRC_DIR = .
INCLUDE_DIRECTORIES = .

# CC = mpicc

CC = icc

EXEC = part

# CFLAGS = -O3 -axSSE4.2 -c -w -traceback -g
CFLAGS = -O1 -c -w -traceback -g

C_CMD	=	$(CC) $(CFLAGS)
C_CMD += $(addprefix -I ,$(INCLUDE_DIRECTORIES))

LINKER_C	= $(CC) $(LFLAGS)

MAIN = main.o

C_OBJS = lag_map.o \
	 handle_stoch_particles.o \
	 $(MAIN)

partMpi:	Makefile $(C_OBJS) 
		$(LINKER_C) *.o -o $(EXEC)

main.o: main.c
	$(C_CMD) main.c
lag_map.o: lag_map.c
	$(C_CMD) lag_map.c
handle_stoch_particles.o: handle_stoch_particles.c
	$(C_CMD) handle_stoch_particles.c
	
clean:
	rm -f Makefile~
	rm -f *.f90~
	rm -f *.mod
	rm -f *.o
	rm -f *.h~
	rm -f *.c~
	rm -f $(EXE_DIR)/$(EXEC)
	rm -f *.vtk
run:
	$(EXE_DIR)/$(EXEC)

