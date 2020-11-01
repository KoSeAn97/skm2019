import plotly.graph_objects as go
from plotly.offline import plot
import pandas as pd
import numpy as np

name = 'solution_0-0'
n = 160
m = 160

z = []
for i in range(n):
    z.append([0] * m)

with open(name, 'r') as f:
    for line in f:
        i, j, v = line.split(',')
        z[int(i)][int(j)] = float(v)

x = np.linspace(-1, 2, n)
y = np.linspace(-2, 2, m)

fig = go.Figure(
    data=[
        go.Surface(
            x=x,
            y=y,
            z=z
        )
    ]
)
fig.update_layout(
    scene = dict(
        xaxis_title='x',
        yaxis_title='y',
        zaxis_title='z'
    )
)
plot(fig, filename=name + '.html')
