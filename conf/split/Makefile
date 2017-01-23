CC=gcc

#if you want run it on a machine with opposite endian to that of the configuration:
ENDIAN=-DMACHINE_ENDIAN_IS_OPPOSITE
#ENDIAN=

CFLAGS=-std=c99 $(ENDIAN)


PROGRAM=split_conf
OBJS=split_conf.o

all : $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $< -o $@

#.o: .c
#	$(CC) $(CFLAGS)  -c $<

clean : 
	rm -f $(PROGRAM) $(OBJS)

