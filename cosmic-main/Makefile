CC = g++
CFLAGS = -c -O2 -W -Wall
CCLDR = g++
OBJ_FLAGS = -G
LDFLAGS = -O2
CLIBS = -L./ -l$(U) -lm
OBJS = cosmic.o
SRCS = cosmic.cc

U = cosmo
.cc.o:
	$(CC) $(CFLAGS) $<

all: lib$(U).a cosmic

cosmic: $< cosmic.o lib$(U).a
	$(CCLDR) $(LDFLAGS) -o cosmic $(OBJS) $(CLIBS)

lib$(U).a: $(U).o
	ar -cr lib$(U).a $(U).o

clean:
	rm -f *.o *.l

distclean:
	rm -f *.o *.l libcosmo.a cosmic

cosmo.o: $(U).cc $(U).h
