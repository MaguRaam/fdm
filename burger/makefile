CC = g++
CFLAGS = -std=c++17 -Wall -O2

EXEC = burger
SRCS = burger.cpp
OBJ = $(SRCS:.cpp=.o)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean distclean

clean:
	@rm -f $(OBJ) $(EXEC)

distclean: clean
	@rm -f data/* plot/* *.dat
