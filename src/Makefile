CC = gcc
UNAME:=$(shell uname -s)
CFLAGS = -Wall -Werror -Wextra -std=c11
GCOVFLAGS = -lgcov --coverage
SOURCES = $(wildcard s21_*.c)
OBJECTS = $(SOURCES:.c=.o)

all := $(shell find ./ -iname "*.c" )

all: s21_matrix.a

rebuild: clean all

s21_matrix.a: $(OBJECTS) s21_matrix.h 
	ar rc $@ $(OBJECTS)
	ranlib $@

lib_cov.a: $(SOURCES) s21_matrix.h
	$(CC) -c *.c -g $(GCOVFLAGS)
	ar rc $@ *.o
	ranlib $@

gcov_report: clean lib_cov.a
	$(CC) test_matrix.c -c -g 
	$(CC) test_matrix.o -g -L. -lcheck lib_cov.a -o test $(GCOVFLAGS)
	-./test 
	lcov -t "t" -o test.info -c -d .
	genhtml -o report test.info 
	open ./report/index.html
	rm -f test

test: clean $(TESTS_SOURSES) s21_matrix.a
	$(CC) -g -o test test_matrix.c s21_matrix.a -lcheck
	./test

clean:
	rm -rf *.a *.o *.test *.gcda *.gcno *.info report *.input *.output */*.exe *.err *.result

clang_format:
	clang-format -i -style=Google *.c *.h
	clang-format -n -style=Google *.c *.h

leak:
ifeq ($(UNAME),Darwin)
	CK_FORK=no leaks -atExit -- ./test
else
	valgrind --leak-check=full -s --track-origins=yes ./test
endif

cppcheck: $(all) 
	cppcheck $^;
