from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from skimage.draw import ellipse


def tikhonov(signals):
    """Apply the Tikhonov regularization algorithm for a given set of measurements from tomography.

    input:
        signals: array or list of arrays
            should be an array of single measurements from each sensor ordered like in "projections",
            or a list of such arrays
        alpha_x, alpha_y, alpha_z, alpha_norm: float
            regularization weights. DerivativeHorizontal derivative. Vertical derivative. Outside Norm. Inside Norm.
        max_iterations: int
            Maximum number of iterations before algorithm gives up

    output:
        g_list: array or list of arrays
            Reconstructed g vector, or multiple g vectors if multiple signals were provided
        first_g: array
            First iteration g vector. This is the result of the simple Tikhonov regularization
    """

    # -------------------------------------------------------------------------

    fname = 'projections.npy'
    print('Reading:', fname)
    projections = np.load(fname)

    print('projections:', projections.shape, projections.dtype)

    # -------------------------------------------------------------------------

    P = projections.reshape((projections.shape[0], -1))

    print('P:', P.shape, P.dtype)

    # -------------------------------------------------------------------------

    n_cuts = projections.shape[1]
    n_rows = projections.shape[2]
    n_cols = projections.shape[3]

    res_x = 200 / float(n_cols)
    res_y = 200 / float(n_rows)
    res_z = 200 / float(n_cuts)  # x,y,z (mm)

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

    Dx = np.array(x_gradient) / res_x
    Dy = np.array(y_gradient) / res_y
    Dz = np.array(z_gradient) / res_z

    print('Dx:', Dx.shape, Dx.dtype)
    print('Dy:', Dy.shape, Dy.dtype)
    print('Dz:', Dz.shape, Dz.dtype)

    # Masks, negative mask: zeros inside, positive mask: zeros outside -------------------------------

    mask_radius = 85.  # Outside this radius the plasma is expected to have zero emissivity

    ii, jj = ellipse(n_rows / 2., n_cols / 2., mask_radius / res_y, mask_radius / res_x)
    mask_negative = np.ones((n_cuts, n_rows, n_cols))
    mask_negative[:, ii, jj] = 0.
    mask_positive = np.zeros((n_cuts, n_rows, n_cols))
    mask_positive[:, ii, jj] = 1.

    Io = np.eye(n_cuts*n_rows*n_cols) * mask_negative.flatten()

    print('Io:', Io.shape, Io.dtype)

    # -------------------------------------------------------------------------

    Pt = np.transpose(P)
    PtP = np.dot(Pt, P)

    DtDx = np.dot(np.transpose(Dx), Dx)
    DtDy = np.dot(np.transpose(Dy), Dy)
    DtDz = np.dot(np.transpose(Dz), Dz)

    ItIo = np.dot(np.transpose(Io), Io)

    alpha_x = 1e5
    alpha_y = alpha_x
    alpha_z = alpha_x

    alpha_norm = alpha_x * 10

    inv = np.linalg.inv(PtP + alpha_x*DtDx + alpha_y*DtDy + alpha_z*DtDz + alpha_norm*ItIo)

    M = np.dot(inv, Pt)

    #  -------------------------------------------------------------------------------------------------

    f = np.array(signals)

    g = np.dot(M, f).reshape((n_cuts, n_rows, n_cols))

    fig, axes = plt.subplots(1, len(g))

    for ax, cross_section in zip(axes, g):
        ax.imshow(cross_section)

    plt.show()

    return g


signals = np.load("signals.npy")
g = tikhonov(signals)
# fig, axes = plt.subplots(1, len(g))
# for g_cut, ax in zip(g, axes):
#     ax.imshow(g_cut)
