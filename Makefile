
TARGET	 = CAP_v6.0.1

OBJS     = main.o input.o calculate.o cluster.o volume_cry.o output.o conformation.o

HEADERS  = crystal.h

CC       = g++
CFLAGS   = -O3 -Wall
LDFLAGS  =
MKFILE   = Makefile

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	@-rm -f $(OBJS)

# dependencies
main.o: main.cpp crystal.h 
input.o: input.cpp crystal.h 
calculate.o: calculate.cpp crystal.h 
cluster.o: cluster.cpp crystal.h 
volume_cry.o: volume_cry.cpp crystal.h 
output.o: output.cpp crystal.h 
conformation.o: conformation.cpp crystal.h 