DEBUG=yes
CC=gcc
FC=gfortran

ifeq ($(DEBUG),yes)
	FFLAGS= -g -Wall
	CFLAGS= -g -I/opt/local/include/ -Wall 
	LDFLAGS= -L/opt/local/lib -lm -ldf -lmfhdf /usr/lib/gcc/i686-apple-darwin11/4.2.1/x86_64/libgfortran.a
else
	FFLAGS= -Wall
	CFLAGS= -I/opt/local/include/ -O3 -ansi 
	LDFLAGS= -L/opt/local/lib -lm -ldf -lmfhdf /usr/lib/gcc/i686-apple-darwin11/4.2.1/x86_64/libgfortran.a
endif

all: edMODISA

edMODISA: edMODISA.o function_read_L3BIN_MODIS_Atmosphere.o subroutine_Ed0moins_at_pixel.o
	$(CC) edMODISA.o function_read_L3BIN_MODIS_Atmosphere.o subroutine_Ed0moins_at_pixel.o $(CFLAGS)  $(LDFLAGS) -o edMODISA

edMODISA.o: edMODISA.c function_read_L3BIN_MODIS_Atmosphere.h
	$(CC) $(CFLAGS) edMODISA.c -c

function_read_L3BIN_MODIS_Atmosphere.o: function_read_L3BIN_MODIS_Atmosphere.c function_read_L3BIN_MODIS_Atmosphere.h
	$(CC) -c $(CFLAGS) function_read_L3BIN_MODIS_Atmosphere.c

subroutine_Ed0moins_at_pixel.o: subroutine_Ed0moins_at_pixel.f
	$(FC) -c $(FFLAGS) subroutine_Ed0moins_at_pixel.f
