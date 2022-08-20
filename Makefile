CC=g++ --std=c++17

WARNINGS= -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

FLAGS= -fsanitize=address

test:prog clean run

prog: main.o Matrix.o
	$(CC) $(FLAGS) $(WARNINGS) -o prog main.o  Matrix.o -lm

main.o:
	$(CC) $(FLAGS) $(WARNINGS) -c main.cpp

Matrix.o:
	$(CC) $(FLAGS) $(WARNINGS) -c Matrix.cpp

.PHONY: clean

clean:
	rm -rf *.o

run:
	./prog