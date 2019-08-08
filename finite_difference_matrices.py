import numpy as np

n_cuts = 4
n_rows = 4
n_cols = 4

# xx finite difference -------------------------------

x_gradient = []

for h in range(n_cuts):
    for i in range(n_rows):
        for j in range(n_cols):
            x = np.zeros((n_cuts, n_rows, n_cols))
            if 0 < j:
                x[h, i, j-1] = -1.
                x[h, i, j] = 1.
            else:  # If j == 0
                x[h, i, j] = -1.
                x[h, i, j+1] = 1.
            x_gradient.append(x.flatten())

# yy finite difference -------------------------------

y_gradient = []

for h in range(n_cuts):
    for i in range(n_rows):
        for j in range(n_cols):
            y = np.zeros((n_cuts, n_rows, n_cols))
            if 0 < i:
                y[h, i-1, j] = -1.
                y[h, i, j] = 1.
            else:  # If i == 0
                y[h, i, j] = -1.
                y[h, i+1, j] = 1.
            y_gradient.append(y.flatten())

# zz finite difference -------------------------------

z_gradient = []

for h in range(n_cuts):
    for i in range(n_rows):
        for j in range(n_cols):
            z = np.zeros((n_cuts, n_rows, n_cols))
            if 0 < h:
                z[h-1, i, j] = -1.
                z[h, i, j] = 1.
            else:  # If h == 0
                z[h, i, j] = -1.
                z[h+1, i, j] = 1.
            z_gradient.append(z.flatten())
