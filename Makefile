#******************************************************
# *         Pragram Name: Makefile                     *
# *         Authore: Junhong Chen                      *
# *         Purpose:My PUI Simulation                  *
# *         Updated time: 1st Oct 2013                 *
# ******************************************************


CC = mpicc -fopenmp
CFLAGS = -std=c99 -O -Wall -g 
GSL = -lgsl -lgslcblas -fopenmp

all: pui
pui: pui_simulation.o error.o hotgas.o menu.o pui_dis.o functions.o e_imp.o fov_int.o print.o set_option.o fitting.o
	$(CC) $^ -o $@ $(GSL)


*.o: commonincludes.h


clean:
	rm -f *.o pui

tar:
	tar cf swics_simulation_Oct232013.tar *.c *.h Makefile
