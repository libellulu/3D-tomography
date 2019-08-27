import numpy as np
import matplotlib.pyplot as plt


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


if __name__ =='__main__':
    # Test plot --------------------------------------------
    # In this example the distribution starts in the center,
    # and then moves to the top right corner while shrinking

    x_min, x_max = (-100, 100)
    y_min, y_max = (-100, 100)
    z_min, z_max = (-0.2, 0.2)

    x_points, y_points, z_points = (30,30, 3)

    x_array = np.linspace(x_min, x_max, x_points)
    y_array = np.linspace(y_max, y_min, y_points)
    z_array = np.linspace(z_min, z_max, z_points)

    scalar_field_coordinates = np.meshgrid(z_array, y_array, x_array, indexing='ij')



    x_coordinates = scalar_field_coordinates[2].flatten()
    y_coordinates = scalar_field_coordinates[1].flatten()
    z_coordinates = scalar_field_coordinates[0].flatten()

    scalar_field_values = gaussian_zz_cone(x=x_coordinates,
                                           y=y_coordinates,
                                           z=z_coordinates,
                                           base_1=-2, base_2=2,
                                           mu_x_1=0, mu_x_2=0,
                                           mu_y_1=0, mu_y_2=0,
                                           sigma_x_1=50, sigma_x_2=50,
                                           sigma_y_1=50, sigma_y_2=50,
                                           height_1=1, height_2=1)


    scalar_field_values += gaussian_zz_cone(x=x_coordinates,
                                            y=y_coordinates,
                                            z=z_coordinates,
                                            base_1=-2, base_2=2,
                                            mu_x_1=70, mu_x_2=-50,
                                            mu_y_1=80, mu_y_2=-30,
                                            sigma_x_1=10, sigma_x_2=10,
                                            sigma_y_1=10, sigma_y_2=10,
                                            height_1=1, height_2=1)

    scalar_field_values[x_coordinates**2 + y_coordinates**2 > 85**2] = 0.0

    scalar_field_values = scalar_field_values.reshape((z_points, y_points, x_points))

    fig, axes = plt.subplots(1, len(scalar_field_values))
    # print('lentgh',len(scalar_field_values[0]))
    # for i in range(0,len(scalar_field_values[0])):
    #     print('premier', scalar_field_values[1][i])
    #     if np.sqrt(scalar_field_values[1][i]**2+scalar_field_values[2][i]**2)>80.5:
    #         scalar_field_values[1][i]=0
    #         scalar_field_values[2][i]=0

    for ax, cross_section in zip(axes, scalar_field_values):
        #ax.imshow(cross_section)
        ax.pcolormesh(np.linspace(-100,100,x_points+1),np.linspace(100,-100,y_points+1),cross_section)

        ax.set_aspect('equal', adjustable='box')

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
    np.save("phantom.npy", scalar_field_values)
