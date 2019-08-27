import numpy as np
from test_gaussian import gaussian_zz_cone
from algo_tomo_new import final_function
from algo_tomo_new import max_z_among_all_CCD
from reconstructions import tikhonov
import matplotlib.pyplot as plt
from collections import defaultdict

def comparison(plotphantom,plotreconstruction,alpha):
    """
    Function that is faster than trial when there is just need to chznge the
    alpha. The signal and the projections have to be done before by runing first
    algo_tomo_new.py with the fixed parameter then test_gaussian_then
    reconstruction

    Parameters
    ----------
    plotphantom: 0 if we do not want to print the phantom, 1 if we want
    plotreconstruction : 0 if we do not want to plot the reconstruction, 1
    if we want .
    alpha : the parameter to twick for a better accuracy

    Returns
    -------
    g : the value of the reconstruction
    scalar_field_values : the value of the phantom


    """

    scalar_field_values= np.load("phantom.npy")

    if plotphantom==1:

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


    g = tikhonov(signals,alpha)

    if plotreconstruction==1:

        fig, axes = plt.subplots(1, len(g))
        for g_cut, ax in zip(g, axes):
            ax.imshow(g_cut)
        plt.show()
    return g,scalar_field_values

def unitest():
    """Function for when we want to try a specific combination of parameter and
    we just have to tune the alpha
    Will print the best alpha

    Parameters
    ----------
    None, the alpha and the parameters are changed in the declaration of the
    function and we call it without parameters

    Returns
    -------
    accuracy : the best accuracy
    index_max : the index in the lists where we have the parameters for this
    best accuracy
    """
    accuracy_list= []

    alpha_list=[]

    for alpha in [1295,1305]:
        g,plasma=comparison(0,0,alpha)
        accuracy=1-(np.sum(np.abs(g-plasma))/(np.sum(plasma)))
        accuracy_list.append(accuracy)

        alpha_list.append(alpha)
        print('alpha', alpha)
        print('mylist',accuracy_list)

    maximum=max(accuracy_list)
    print('max array',maximum)
    index_max=accuracy_list.index(maximum)
    print('index',index_max)
    print('optimum parameters','alpha=',alpha_list[index_max],'accurac=',accuracy_list[index_max])
    return maximum,index_max


    return accuracy

unitest()
