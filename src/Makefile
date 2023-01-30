# This makefile compiles the peerMethods routines
# and creates the corrisponding shared library

dlvl = ../.
include $(dlvl)/Makefile.in

OBJECTS = CLab.o utilities.o peerMethods.o

.PHONY: all

all: cleanall lib clean

lib: ${OBJECTS}
	${CC} ${CFLAGS} ${ARCHFLAGS} -o ${CLIB} ${OBJECTS} -Wl,--no-whole-archive

peerMethods.o: peerMethods.c
	${CC} -I../include ${CFLAGS} -c $< ${LIB_FLAGS}

CLab.o: CLab.c
	${CC} -I../include ${CFLAGS} -c $<

utilities.o: utilities.c
	${CC} -I../include ${CFLAGS} -c $<

.PHONY: clean

clean:
	rm -f *.o

rmlib:
	rm -f ${CLIB}

cleanall: clean rmlib
	rm -f *.out