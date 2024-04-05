CXX = g++

CXXFLAGS = -std=c++20 -Wall

TARGET = main

SRCS = main.cpp definitions.cpp

HEADERS = definitions.hpp

OBJS = $(SRCS:.cpp=.o)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)