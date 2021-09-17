CC = gcc
RM = rm -rf
LAPACK_DIR = /opt/libflame/gcc-7.5.0/2.2/
OMP_DIR = /opt/openmpi/gcc-7.5.0/4.1.0/
CFLAGS = -g -O2 -Wall -mcmodel=medium -I$(LAPACK_DIR)/include -I$(OMP_DIR)/include -fopenmp -llapack -lm
LDFLAGS = -mcmodel=medium -L$(LAPACK_DIR)/lib -L$(OMP_DIR)/lib -fopenmp 
LINKS = -llapack -lm
OBJS = srofm.o
TARGET = srofm

.PHONY: all clean dep

all : $(TARGET)

clean :
	$(RM) *.o
	$(RM) $(TARGET)

dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c)

$(TARGET) : $(OBJS)
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ $(OBJS) $(LINKS)

.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<
