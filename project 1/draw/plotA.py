import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata


boundaryCondition = {
    0: 'Dirichlet',
    1: 'Neumann',
    2: 'Mixed'
}

domain = {
    0: 'regular',
    1: 'irregular'
}

path = os.path.join('output', 'output.json')
with open(path, 'r') as file:
    data = json.load(file)


for i, j in enumerate(data):
    grids = np.array(j["grids"])
    values = np.array(j["values_on_grids"])
    realValues = np.array(j["real_values"])
    x = grids[:, 0]
    y = grids[:, 1]
    bc = boundaryCondition.get(j["boundary_condition"])
    d = domain.get(j["Domain"])
    z_values = values
    z_real_values = realValues

    fig = plt.figure(figsize=(20, 8))  
    ax1 = fig.add_subplot(121, projection='3d')  
    ax2 = fig.add_subplot(122, projection='3d')  


    ax1.plot_trisurf(x, y, z_real_values, cmap='Blues', alpha=0.7)
    ax1.scatter(x, y, z_real_values, color='red', label='real values')
    ax1.set_title(f"Real Values (domain:{d}, bc:{bc}, n={2**(3 + i % 4)}")
    ax1.set_xlabel("grids x")
    ax1.set_ylabel("grids y")
    ax1.set_zlabel("values")
    ax1.legend()

    ax2.plot_trisurf(x, y, z_values, cmap='YlOrBr', alpha=0.7)
    ax2.scatter(x, y, z_values, color='blue', label='numerical values')
    ax2.set_title(f"Numerical Values (domain:{d}, bc:{bc}, n={2**(3 + i % 4)}")
    ax2.set_xlabel("grids x")
    ax2.set_ylabel("grids y")
    ax2.set_zlabel("values")
    ax2.legend()

    save = os.path.join('figure', f"functype_{1 + i // 24}_{bc}_{d}_{i + 1}.png")
    plt.savefig(save)
    plt.close()
    print(f"figure saved: {save}")
