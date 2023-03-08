CC = g++
CFLAGS = -std=c++17 -Wall -O2
LDFLAGS = -L/home/magu/library/armadillo-11.4.4/install/lib -larmadillo -llapack -lsuperlu
INCLUDE = -I/home/magu/library/armadillo-11.4.4/include

EXEC = poisson
SRCS = poisson.cpp
OBJ = $(SRCS:.cpp=.o)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean distclean

clean:
	rm -f $(OBJ) $(EXEC)

distclean: clean
	rm -f *.dat