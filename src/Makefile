# This makefile compiles the peerMethods routines
# and creates the corrisponding shared library

dlvl = ../.
include $(dlvl)/Makefile.in

OBJECTS = CLab.o utilities.o peerMethods.o

.PHONY: all

all: cleaner lib clean

lib: ${OBJECTS}
	${CC} ${CFLAGS} ${ARCHFLAGS} -o ${CLIB} ${OBJECTS} -Wl,--no-whole-archive

peerMethods.o: peerMethods.c
	${CC} -I../include ${INCLUDE_DIRS} ${CFLAGS} -c $< ${LINKED_LIB} ${LIB_FLAGS}

CLab.o: CLab.c
	${CC} -I../include ${CFLAGS} -c $<

utilities.o: utilities.c
	${CC} -I../include ${CFLAGS} -c $<

.PHONY: clean

clean:
	rm -f *.o

rmlib:
	rm -f ${CLIB}

cleaner: clean rmlib
	rm -f *.out