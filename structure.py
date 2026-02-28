"""
Graphene Tight-Binding (6-site hexagon fragment)

- Defines a lattice fragment as a graph (sites + edges)
- Builds the tight-binding Hamiltonian H (nearest-neighbor hopping)
- Computes the energy spectrum (eigenvalues)
- Prints a clean report and optionally saves outputs (CSV/JSON)
"""


import numpy as np 
import matplotlib.pyplot as plt


#atoms of the hexagon fragment and t-parameter init
N = 6
t = 1.0
#edges (neighbors)
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0)]

# coordinates for drawing

coords = {
    0: (1.0, 0.0),
    1: (0.5, np.sqrt(3)/2),
    2: (-0.5, np.sqrt(3)/2),
    3: (-1.0, 0.0),
    4: (-0.5, -np.sqrt(3)/2),
    5: (0.5, -np.sqrt(3)/2),
}

for i, j in edges:
    x = [coords[i][0], coords[j][0]]
    y = [coords[i][1], coords[j][1]]
    plt.plot(x, y)

# draw atoms
for i in coords:
    plt.scatter(coords[i][0], coords[i][1])
    plt.text(coords[i][0], coords[i][1], str(i))

plt.gca().set_aspect('equal')
plt.title("Graphene Hexagon Fragment")
plt.show()

print("Sites:", list(range(N)))
print("Edges:", edges)