# 3D-tomography

This project simulates a 3D tomography set up in the tokamak of ISTTOK. By changing parameters, it is possible to simulate 3D tomography in any other tokamak.

This project help you to find the best configuration to implement a 3D tomography system inside a tokamak.


## Getting started

This README is going to explain you how the code can be used an adapted to each personnal case (for other tokamaks and set up than ISTTOK). Then inside of each files, the functions are described (Goal, Parameters, Output).

The parameters that are passed in the function when you first upload them are consistent with ISTTOK set up.

### Prerequisites

You need to have the packages for python 3 installed.

## First stage, simulate the tokamak of you choice 

By editing algo_tomo_new.py you can change inside of the file the parameters to adapt. 
For example, the minor radius of the tokamak and the place of the camera inside of it.
Running it will print a 3D plot of the tokamak and the lines of sight. 
![Lines of sight of 3D set up](https://github.com/libellulu/3D-tomography/blob/master/images/LOS_cylinderclose.png)

You can also check the coverage of the space by the lines of sight with the file space_coverage.py . The only thing to do is to pass the same parameter in final function that you did in algo_tomo_new for the result to be consitent.
The two following pictures show different examples. The first one shows the coverage of every lines of sight, and the second picture shows the coverage for the top camera only. To verify the coverage more specifically like on the second picture, you need to know how many LOS each camera has. For example if each camera have 30 LOS, then in the code to see the full coverage, you have to put in linspace from 0 to 89 (90 LOS in total). To see only the top from 0 to 29, the second camera is from 30 to 59... 

![Coverage for every camera](https://github.com/libellulu/3D-tomography/blob/master/images/20%2C11%2C065%2C004%2C50%2C50%2C3.png)

![Coverage for the top camera](https://github.com/libellulu/3D-tomography/blob/master/images/dist_pinhole%3D9.png)

You can simulate quite a lot with space_coverage.py. You can also decide to simulate only one row of sensors. If you want to check the central row of the second CCD for example, following the previous example, it would be in linspace from 40 to 49.

The black circle represents the limiter in ISTTOK, radius 85. It enables to see where the plasma is in comparison with the lines of sight. You can decide not to print the limiter or change its value.

## Second stage, try to find the best reconstruction with the phantom evolving in space

My advice is to begin by using the most simple plasma phantom of test_gaussian_V1.py. This phantom shows the plasma at the same time in 3 different places.The second version test_gaussian_V2 shows the 3 places at 5 different times.

When running the program, if you chose the parameter printplot=1 (see explanation in the code) to print the phantom and the reconstruction, here is the phantom you should get :

![Phantom of V1](https://github.com/libellulu/3D-tomography/blob/master/images/phantom2_zoom.png)

Use the file trial_function_V1 to determine the parameter that give the best reconstruction. The different functions inside of the file are each looping over one parameter which you define the range , and return the configuration that gives the best accuracy of reconstruction. 

Once the best set up is found, use the function unitest (still inside of test_gaussian_V1.py) to run only over the alpha. The alpha parameter can take any value above zero and is very sensitive. So a lot of tests have to be done....

Notice that this program charges the projection.npy and the signals.npy files every time. Once you decide to look for the best alpha, you can comment this part for the program to be faster.

## Third stage, using the space and time phantom

This phantom is used in trial_function_V2. This function only has unitest because normally, one should just try to find the best alpha.

This part is very tricky because the alpha are different for each time. The goal is to find a compromise. 

## Tips 
The trial function can be very long to run if you try a lot of configurations. I solved this problem by avoiding the print and using a server, more powerful than my PC. If you are just running over one final parameter, you can avoid charging projection.npy and signal.npy, it saves computing power.  

To run over the alpha, if you have no idea, try first by digit 0.1,1,10,100,1000, it should give you a clearer idea of the range you should look into.

