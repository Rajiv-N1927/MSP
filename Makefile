CC=g++
CFLAGS=-I.
DEPS = io_bmp.h io_custom.h
OBJ = io_bmp.o example_main.o io_custom.o

%.o: %.cpp $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

io_bmp: $(OBJ)
	$(CC) -g -o $@ $^ $(CFLAGS)
