import numpy as np
from scipy.sparse import lil_matrix

def build_and_print_poisson_matrix():
    N = 5  # 网格尺寸 (5x5=25个点)
    total_points = N * N
    A = lil_matrix((total_points, total_points), dtype=np.float64)
    h = 1.0 / (N - 1)  # 网格间距

    # 构建矩阵
    for i in range(N):
        for j in range(N):
            index = i * N + j
            A[index, index] = 4
            
            # x方向邻居
            if i > 0:
                A[index, index - N] = -1
            if i < N - 1:
                A[index, index + N] = -1 
            
            # y方向邻居
            if j > 0:
                A[index, index - 1] = -1
            if j < N - 1:
                A[index, index + 1] = -1
            
            # Neumann边界调整
            if i == 0 or i == N - 1:
                A[index, index] += -1 
            if j == 0 or j == N - 1:
                A[index, index] += -1 

    # 转换为稠密矩阵并打印
    A_dense = A.toarray()
    for row in A_dense:
        print(' '.join(f"{int(val) if val.is_integer() else val:.1f}" for val in row))

# 执行打印
build_and_print_poisson_matrix()