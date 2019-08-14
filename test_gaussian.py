import numpy as np
import matplotlib.pyplot as plt


def gaussian_zz_cone(x, y, z,
                     base_1, base_2, mu_x_1, mu_x_2, mu_y_1, mu_y_2, sigma_x_1, sigma_x_2, sigma_y_1, sigma_y_2):
    """Scalar field of gaussian distributions in the xy plane overlapping in the zz direction. Like a cone.
    The cone is defined by its two bases, `base_1` and `base_2`.
    Base one has center (`mu_x_1`, `mu_y_1`) and standard deviation (`sigma_x_1`, `sigma_y_1`).
    The same goes for base two.
    The center and standard deviation of each distribution change linearly from one base to the other

    Parameters
    ----------
    x, y, z: float
        Point at which the scalar field is computed.
    base_1, base_2: float
        Z-coordinate of base one and two respectively
    mu_x_1, mu_x_2: float
        X-coordinate of the center of base one and two.
    mu_y_1, mu_y_2: float
        Y-coordinate of the center of base one and two.
    sigma_x_1, sigma_x_2: float
        Standard deviation in xx for bases one and two.
    sigma_y_1, sigma_y_2: float
        Standard deviation in yy for bases one and two.
    """

    # mu_x changes linearly between mu_x_1 and mu_x_2
    k = (mu_x_2 - mu_x_1) / (base_2 - base_1)
    b = mu_x_1 - k * base_1
    mu_x = k * z + b

    # mu_y changes linearly between mu_y_1 and mu_y_2
    k = (mu_y_2 - mu_y_1) / (base_2 - base_1)
    b = mu_y_1 - k * base_1
    mu_y = k * z + b

    # sigma_x changes linearly between sigma_x_1 and sigma_x_2
    k = (sigma_x_2 - sigma_x_1) / (base_2 - base_1)
    b = sigma_x_1 - k * base_1
    sigma_x = k * z + b

    # sigma_y changes linearly between sigma_y_1 and sigma_y_2
    k = (sigma_y_2 - sigma_y_1) / (base_2 - base_1)
    b = sigma_y_1 - k * base_1
    sigma_y = k * z + b

    n_x = 1. / np.sqrt( 2 * np.pi * sigma_x ** 2)
    n_y = 1. / np.sqrt( 2 * np.pi * sigma_y ** 2)

    g = n_x * n_y * np.exp( - (x - mu_x) ** 2 / sigma_x ** 2 - (y - mu_y) ** 2 / sigma_y ** 2)

    return g

if __name__ =='__main__':
    # Test plot --------------------------------------------
    # In this example the distribution starts in the center,
    # and then moves to the top right corner while shrinking

    x_min, x_max = (-100, 100)
    y_min, y_max = (-100, 100)
    z_min, z_max = (-148, 148)

    x_points, y_points, z_points = (15, 15, 5)

    x_array = np.linspace(x_min, x_max, x_points)
    y_array = np.linspace(y_max, y_min, y_points)
    z_array = np.linspace(z_min, z_max, z_points)

    scalar_field_coordinates = np.meshgrid(z_array, y_array, x_array, indexing='ij')

    scalar_field_values = gaussian_zz_cone(x=scalar_field_coordinates[2].flatten(),
                                           y=scalar_field_coordinates[1].flatten(),
                                           z=scalar_field_coordinates[0].flatten(),
                                           base_1=-150, base_2=150,
                                           mu_x_1=0, mu_x_2=0,
                                           mu_y_1=0, mu_y_2=0,
                                           sigma_x_1=50, sigma_x_2=15,
                                           sigma_y_1=50, sigma_y_2=15)

    scalar_field_values = scalar_field_values.reshape((z_points, y_points, x_points))

    fig, axes = plt.subplots(1, len(scalar_field_values))

    for ax, cross_section in zip(axes, scalar_field_values):
        ax.imshow(cross_section)

    plt.show()

    # Simulate tomography signals --------------------------------------------------

    fname = 'projections.npy'
    print('Reading:', fname)
    projections = np.load(fname)

    print('projections:', projections.shape, projections.dtype)

    # -------------------------------------------------------------------------

    P = projections.reshape((projections.shape[0], -1))

    print('P:', P.shape, P.dtype)

    signals = np.dot(P, scalar_field_values.flatten())

    print('Signals:', signals.shape, signals.dtype)

    np.save("signals.npy", signals)
