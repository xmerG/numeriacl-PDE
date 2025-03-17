CXX = g++

CXXFLAGS = -w -std=c++14 -O2  

LDFLAGS = -llapacke -llapack -lblas

TARGET = test

SRCS = test.cpp Circle.cpp TestFunction.cpp

OBJS = $(SRCS:.cpp=.o)


all: $(TARGET)
	./$(TARGET) 

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean