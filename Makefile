

all: libutepnum.a

integers/integer_ops.o: integers/integer_ops.c include/integer_ops.h
	gcc -Iinclude -Wall -O0 -g -c integers/integer_ops.c -o $@

libutepnum.a: integers/integer_ops.o
	ar -rv $@ $^

test: tests/test_integers
	tests/test_integers 1 1 17 42
	tests/test_integers 2 3 99999999999999999999917 888888888888888842

tests/test_integers: libutepnum.a tests/test_integers.o
	gcc -Iinclude -L. -Wall -O0 -g -o $@ tests/test_integers.o libutepnum.a

tests/test_integers.o: tests/test_integers.c include/utepnum.h
	gcc -Iinclude -Wall -O0 -g -c tests/test_integers.c -o $@

clean:
	rm -f libutepnum.a
	rm -f integers/integer_ops.o
	rm -f tests/test_integers.o
	rm -f tests/test_integers


.PHONY: all clean test
