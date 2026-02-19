CC   = gcc
LD   = gcc

INC_DIR = /usr/local/include

CFLAGS = -std=gnu99 -O3 -DHAVE_INLINE -Wall -I $(INC_DIR)
LDLIBS = -L/usr/local/lib -lgsl -lgslcblas -lm

COMMON_OBJS = signal_tools.o

TARGETS = fft_gauss fft_AM fft_decode

all: $(TARGETS)

fft_gauss: fft_gauss.o $(COMMON_OBJS)
	$(LD) $^ $(LDLIBS) -o $@

fft_AM: fft_AM.o $(COMMON_OBJS)
	$(LD) $^ $(LDLIBS) -o $@

fft_decode: fft_decode.o $(COMMON_OBJS)
	$(LD) $^ $(LDLIBS) -o $@

%.o: %.c signal_tools.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS)