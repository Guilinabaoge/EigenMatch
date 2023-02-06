import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

G_e = np.loadtxt("graphA.txt", int)
G = nx.Graph(G_e.tolist())
#nx.draw(G)
print("Just checking",nx.is_directed(G))
print(G.nodes())
print(G.edges())
pos = nx.spring_layout(G, seed=225)  # Seed for reproducible layout
#nx.draw(G, pos)

Src_e =np.array(G.edges)
Deg = np.diag((np.sum(Src_e, axis=1)))
#G = nx.from_numpy_array(Src_e)
#nx.draw_networkx(G, pos, **options)

options = {
    "font_size": 36,
    "node_size": 3000,
    "node_color": "white",
    "edgecolors": "green",
    "linewidths": 5,
    "width": 4,
}
#nx.draw_random(G,**options)
#nx.draw_networkx(G)
nx.draw_networkx(G,**options)
# Setnx.draw_networkx(G,pos) margins for the axes so that nodes aren't clipped
ax = plt.gca()
ax.margins(0.20)
plt.axis("off")
plt.show()
