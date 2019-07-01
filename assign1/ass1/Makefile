CC=g++
CFLAGS=-I.
DEPS = io_bmp.h image_comps.h
OBJ = io_bmp.o filtering_main.o io_comps.o

%.o: %.cpp $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

filtering: $(OBJ)
	$(CC) -g -o $@ $^ $(CFLAGS)
