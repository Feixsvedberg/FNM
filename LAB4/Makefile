CC   = gcc
LD   = gcc

INC_DIR = /usr/local/include

CFLAGS = -std=gnu99 -O3 -DHAVE_INLINE -Wall -I $(INC_DIR)
LDLIBS = -L/usr/local/lib -lgsl -lgslcblas -lm

COMMON_OBJS = numerov.o

TARGETS = harm_osc morse morse_rework

all: $(TARGETS)

harm_osc: harm_osc.o $(COMMON_OBJS)
	$(LD) $^ $(LDLIBS) -o $@

morse: morse.o $(COMMON_OBJS)
	$(LD) $^ $(LDLIBS) -o $@
	
morse_rework: morse_rework.o $(COMMON_OBJS)
	$(LD) $^ $(LDLIBS) -o $@

%.o: %.c numerov.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS)