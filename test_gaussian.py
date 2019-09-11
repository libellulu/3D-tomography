import numpy as np
import matplotlib.pyplot as plt


def gaussian_3d(x, y, z, mu_x, mu_y, mu_z, sigma_x, sigma_y, sigma_z, height=None):
    """ Gaussian in 3d space

    Parameters
    ----------
    x, y, z: float or 1darray
        Point at which the scalar field is computed.
    mu_x, mu_y, mu_z: float
        Coordinates of the center of the gaussian
    sigma_x, sigma_y, sigma_z: float
        Standard deviation in each direction
    height: float, optional
        Height of the gaussian function.
        If not specified, the gaussian will default to a normalized version.

    Returns
    -------
        s: float or 1darray
            Scalar field value at coordinates (x, y, z)
    """

    if height is None:  # Normalized gaussian
        n_x = 1. / np.sqrt(2 * np.pi * sigma_x ** 2)
        n_y = 1. / np.sqrt(2 * np.pi * sigma_y ** 2)
        n_z = 1. / np.sqrt(2 * np.pi * sigma_z ** 2)
        s = n_x * n_y * n_z * np.exp(- (x - mu_x) ** 2 / sigma_x ** 2
                                     - (y - mu_y) ** 2 / sigma_y ** 2
                                     - (z - mu_z) ** 2 / sigma_z ** 2)

    else:  # Gaussian with peak height defined
        s = height * np.exp(- (x - mu_x) ** 2 / sigma_x ** 2
                            - (y - mu_y) ** 2 / sigma_y ** 2
                            - (z - mu_z) ** 2 / sigma_z ** 2)
    return s


def flattop(x, y, center_x, center_y, width, height, degree: int = 2):
    """ Define a flattop plasma profile (scalar field) circularly symmetric around `center`.
    The `degree` should be an even number specifying how flat the profile is.

    Parameters
    ----------
        x, y: float or 1darray
            Point at which the scalar field is computed.
        center_x, center_y: float
            x and y coordinates of the center
        width: float
            Total width of the plasma profile. Outside this width the profile is zero
        height: float
            Height of the vertex of the flattop
        degree: even integer
            Degree of the polynomial function used as the profile. The FWHM is a function of the chosen degree.
            FWHM = width / 2**(1/degree) e.g. degree = (   2,    4,    6,    8)
                                                FWHM = (0.70, 0.84, 0.89, 0.92) * width

    Returns
    -------
        s: float or 1darray
            The value of the scalar field at (x, y)
    """
    if (degree % 2 != 0) and (degree < 2):
        raise ValueError("specified degree should be an even integer, value %d given" % degree)

    rho = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    half_width = width / 2

    s = np.zeros_like(x)

    ii = rho < half_width

    s[ii] = - height * (rho[ii] / half_width)**degree + height

    return s


def gaussian_zz_cone(x, y, z,
                     base_1, base_2, mu_x_1, mu_x_2, mu_y_1, mu_y_2, sigma_x_1, sigma_x_2, sigma_y_1, sigma_y_2,height_1=None, height_2=None):
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
    height_1, height_2: float, optional
        Height of the gaussian function. If not specified, the gaussian will
        default to a normalized version.
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

    # height changes linearly between height_1 and height_2
    k = (height_2 - height_1) / (base_2 - base_1)
    b = height_1 - k * base_1
    height = k * z + b

    n_x = 1. / np.sqrt( 2 * np.pi * sigma_x ** 2)
    n_y = 1. / np.sqrt( 2 * np.pi * sigma_y ** 2)

    #if np.sqrt(x**2+y**2)<80.5:

    if (height_1 != None) and (height_2 != None):
        g = height * np.exp( - (x - mu_x) ** 2 / sigma_x ** 2 - (y - mu_y) ** 2 / sigma_y ** 2)
    else:
        g = n_x * n_y * np.exp( - (x - mu_x) ** 2 / sigma_x ** 2 - (y - mu_y) ** 2 / sigma_y ** 2)
    #else:
    #    g=0

    return g


if __name__ == '__main__':

    # Test gaussian_3d and flattop ------------------------------------------------------------

    x_min, x_max = (-100, 100)
    y_min, y_max = (-100, 100)
    z_min, z_max = (-4, 4)

    x_points, y_points, z_points = (30, 30, 3)

    x_array = np.linspace(x_min, x_max, x_points)
    y_array = np.linspace(y_max, y_min, y_points)
    z_array = np.linspace(z_min, z_max, z_points)

    scalar_field_coordinates = np.meshgrid(z_array, y_array, x_array, indexing='ij')

    x_coordinates = scalar_field_coordinates[2].flatten()
    y_coordinates = scalar_field_coordinates[1].flatten()
    z_coordinates = scalar_field_coordinates[0].flatten()

    blob_z = np.linspace(-8, 8, 5)
    phantoms = []

    for z in blob_z:

        phantom = flattop(x=x_coordinates, y=y_coordinates,
                          center_x=0, center_y=0, width=85 * 2, height=1, degree=2)

        phantom += gaussian_3d(x=x_coordinates, y=y_coordinates, z=z_coordinates,
                               mu_x=20, mu_y=20, mu_z=z,
                               sigma_x=5, sigma_y=5, sigma_z=4,
                               height=4)

        phantoms.append(phantom.reshape((z_points, y_points, x_points)))

    fig, axes = plt.subplots(len(phantoms), z_points)

    for sub_axes, phantom in zip(axes, phantoms):
        for ax, cross_section in zip(sub_axes, phantom):
            im = ax.pcolormesh(np.linspace(-100, 100, x_points + 1),
                               np.linspace(100, -100, y_points + 1),
                               cross_section,
                               vmin=0, vmax=2)
            ax.set_aspect('equal', adjustable='box')

    fig.colorbar(im, ax=axes.ravel().tolist())

    plt.show()