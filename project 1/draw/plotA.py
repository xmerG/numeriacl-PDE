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

domain={
    0:'regular',
    1:'irregular'
}

path=os.path.join('..', 'output','test.json')
with open(path, 'r') as file:
    data=json.load(file)

def f(x,y):
    return np.exp(y+np.sin(x))

for i,j in enumerate(data):
    grids=np.array(j["grids"])
    values=np.array(j["values_on_grids"])
    realValues=np.array(j["real_values"])
    x=grids[:,0]
    y=grids[:,1]
    bc=boundaryCondition.get(j["boundary_condition"])
    d=domain.get(j["Domain"])
    z_values=values
    z_real_values=realValues
    fig=plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111,projection='3d')
    ax.scatter(x,y,z_real_values,color='red', label='real values')
    tri1=ax.plot_trisurf(x,y,z_real_values, cmap='viridis', alpha=0.7, label='real values')
    ax.scatter(x,y,z_values, color='blue', label='numerical values')
    tri2=ax.plot_trisurf(x,y,z_values, cmap='viridis', alpha=0.7, label='numerical values')
    ax.set_title(f"domain:{d}, bc:{bc}, n={2**(3+i)}")
    ax.set_xlabel("grids x")
    ax.set_ylabel("grids y")
    ax.set_zlabel("values")
    ax.legend()
    save=os.path.join('..', 'figure',f"f_{bc}_{d}_{i+1}.png")
    plt.savefig(save)
    plt.close()

