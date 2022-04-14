CC = icc
CFLAGS = -O3 -xMIC-AVX512 -DSTREAM_ARRAY_SIZE=80000000 -DNTIMES=20 -qopenmp


all: stream_c.exe


stream_c.exe: stream.c
	$(CC) $(CFLAGS) stream.c -o stream_c.exe

clean:
	rm -f stream_f.exe stream_c.exe *.o
