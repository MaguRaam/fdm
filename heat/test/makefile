CXX = g++
CXXFLAGS = -O2 -I /home/magu/library/eigen
TARGET = laplace
SOURCES = laplace.cpp
OBJECTS = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean distclean

clean:
	@rm -f $(TARGET) $(OBJECTS)

distclean: clean
	@rm -f *.dat *.vtk