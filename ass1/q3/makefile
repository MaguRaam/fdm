CXX = g++
CXXFLAGS = -O3 -std=c++17
TARGET = ode
INCLUDE_PATH = "../"

$(TARGET) : $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ -I $(INCLUDE_PATH)

clean:
	@rm $(TARGET)

distclean:
	@rm $(TARGET) *.dat
