import numpy as np
from test_gaussian import gaussian_zz_cone
from algo_tomo_new import final_function
from algo_tomo_new import max_z_among_all_CCD
from reconstructions import tikhonov
import matplotlib.pyplot as plt

def comparison(printplot,nb_cell_x,nb_cell_z,spacing_x,spacing_z,nb_voxel_x,nb_voxel_y,nb_voxel_z,radius_tokamak,pinhole_to_CCD):

    real_projection,maxz=final_function(nb_cell_x,nb_cell_z,spacing_x,spacing_z,nb_voxel_x,nb_voxel_y,nb_voxel_z,radius_tokamak,pinhole_to_CCD)

    print(maxz)

    # Test plot --------------------------------------------
    # In this example the distribution starts in the center,
    # and then moves to the top right corner while shrinking


    x_min, x_max = (-radius_tokamak, radius_tokamak)
    y_min, y_max = (-radius_tokamak, radius_tokamak)
    z_min, z_max = (-maxz,maxz)

    x_points, y_points, z_points = (nb_voxel_x,nb_voxel_y,nb_voxel_z)

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

    if printplot==1:
        fig, axes = plt.subplots(1, len(scalar_field_values))

        for ax, cross_section in zip(axes, scalar_field_values):
            ax.imshow(cross_section)


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

    signals = np.load("signals.npy")

    g = tikhonov(signals)


    fig, axes = plt.subplots(1, len(g))
    for g_cut, ax in zip(g, axes):
        ax.imshow(g_cut)
    plt.show()
    return real_projection

comparison(1,12,12,7,7,20,20,4,100,20)
# for pinhole_to_CCD in (4,7):
#     print('pinhole',pinhole_to_CCD)
#     comparison(1,5,5,2,2,15,15,5,100,pinhole_to_CCD)

# tuning of the x and z number of cells
#for n in (2,4,6,8,10,12):
#    comparison(1,12,12,7,7,20,20,4,100,20)
