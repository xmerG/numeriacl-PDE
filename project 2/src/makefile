# 定义编译器
CXX := g++
# 定义编译选项
CXXFLAGS := -O2

# 定义目标可执行文件
TARGET := testall

# 定义所有源文件
SRCS := testall.cpp restriction/Full_weighting.cpp restriction/Injection.cpp \
        prolongation/Linear.cpp prolongation/Quadratic.cpp Vector.cpp \
        Multigrid.cpp Sparse_Matrix.cpp Test.cpp

# 生成对应的对象文件列表
OBJS := $(SRCS:.cpp=.o)

# 默认目标：编译并运行
all: $(TARGET)
	@echo "Running $(TARGET)..."
	@./$(TARGET)

# 链接目标
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# 模式规则编译.cpp文件为.o文件
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 显式声明哪些头文件被哪些源文件依赖
Multigrid.o: Multigrid.cpp Multigrid.h
Test.o: Test.cpp Multigrid.h

# 清理生成的文件
clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean