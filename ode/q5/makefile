CXX = g++
CXXFLAGS = -O3 -std=c++17
TARGET = thomas
INCLUDE = -I ../ -I ~/library/eigen

$(TARGET) : $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLUDE)

clean:
	@rm $(TARGET)

distclean:
	@rm $(TARGET) *.dat
