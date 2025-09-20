compile: first second

first:
	clang bigtest.c -Wall -Werror -pedantic -std=c99 -o c-v

second:
	clang++ bigtest.cpp -Wall -Werror -std=c++20 -o cpp-v
