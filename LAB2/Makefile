CC   = gcc
LD   = gcc

INC_DIR = /usr/local/include

CFLAGS = -std=gnu99 -O3 -DHAVE_INLINE -Wall -I $(INC_DIR)
LDLIBS = -L/usr/local/lib -lgsl -lgslcblas -lm

TARGETS = phonons

all: $(TARGETS)

#frequencies är länkad till phonons nu
phonons: phonons.o frequencies.o $(COMMON_OBJS)
	$(LD) $^ $(LDLIBS) -o $@

clean:
	rm -f *.o $(TARGETS)