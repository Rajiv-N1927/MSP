CC=g++
CFLAGS=-I.
DEPS = io_bmp.h image_comps.h
OBJ = io_bmp.o motion_main.o

%.o: %.cpp $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

motion: $(OBJ)
	$(CC) -g -o $@ $^ $(CFLAGS)
