# 3D-tomography
This project simulates a 3D tomography set up in the tokamak of ISTTOK. By changing parameters, it is possible to simulate 3D tomography in any other tokamak.
This project help you to find the best configuration to implement a 3D tomography system inside of a tokamak.


##Getting started
This README is going to explain you how the code can be used an adapted to each personnal case for other tokamak and set up than ISTTOK

### Prerequisites
You need to have the packages for python 3 installed.

##First stage, simulate the tokamak of you choice 
By editing algo_tomo_new.py you can change inside of file the parameter to adapt. 
For example, the minor radius of the tokamak and the place of the camera inside of it.
Running it will printa 3D plot of the tokamak and the lines of sight. 

You can also check the coverage by the lines of sight with the file space_coverage.py . The only thing to do is to pass the same parameter in final function that you did in algo_tomo_new for the result to be consitent.

##Second stage, try to find the best reconstruction with the phantom evolving in space
My advice, begin by using the most simple plasma phantom of test_gaussian_V1. This phantom shows the plasma at the same time in 3 different place.
The second version test_gaussian_V2, shows the 3 places at 5 different times.
Use the file trial_function_V1 to determine the parameter that give the best reconstruction. The different functions inside of the file are each looping aver one parameter which you define the range , and return the configuration that gives the best accuracy of reconstruction. 
Once the best set up is found, use unitest to run only over the alpha. The alpha parameter can take any values and is very sensitive. So a lot of test have to be done....
Notice that this programm charge the projection.npy file every time. Once you decide to look for the best alpha, you can comment this part for the program to be faster.

###Third stage, using the space and time phantom
This phantom is used in trial_function_V2. This function only has unitest because normally, one should just try to find the best alpha.
This part is very tricky because the alpha are different for each time. The goal is to find a compromise. 

