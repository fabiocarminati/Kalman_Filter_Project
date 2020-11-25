### Kalman Filter-Localization Navigation and Smart Mobility Project:

In this project we handle a localization problem by formulating the model, implementing a positioning method and assessing the related performance in **Matlab** environment.

### Getting Started

Clone or download the repository somewhere on your PC.

In the GR35 folder are placed the measurements for the various tasks.

The are 8 Access Points and 1 user location.

Each task could be useful for the next ones thus their order of execution is relevant. 

### Prerequisites

**Matlab** is required.

### Tasks

##### 1) Learning the AP positions and the measurement model:

First we need to find the exact location of the APs (that won't change during the tasks) using measurements with no noise.

<img src="C:\Users\carmi\AppData\Roaming\Typora\typora-user-images\image-20201125173722812.png" alt="image-20201125173722812" style="zoom: 80%;" />

Then the covariance matrix of TOA measurements is computed.

##### 2) Learning the motion model dataset for motion models

 We observe the trajectory dataset and since
1)Ux,Ux,Vx,Vy are given
2)There is a random acceleration process whose mean is 0 while the variance is not 0  we came to the conclusion that this is a M3 motion model: **Random Force model.**

<img src="C:\Users\carmi\AppData\Roaming\Typora\typora-user-images\image-20201125173806414.png" alt="image-20201125173806414" style="zoom:80%;" />

(first 10 M3 user trajectories)

##### 3) Kalman filter testing on known trajectory

Implementation of an **Extended Kalman Filter** over a known user trajectory and using the previous results obtained by tasks 1,2.

<img src="C:\Users\carmi\AppData\Roaming\Typora\typora-user-images\image-20201125173850402.png" alt="image-20201125173850402" style="zoom: 67%;" />

(red is the real trajectory; green is the predicted one by the **Extended Kalman Filter**)

##### 4) Kalman filter testing on unknown trajectory

Implement an **Extended Kalman Filter** by using the variable “rhoUEAP.mat” in “Task4_rhoUEAP_GRXX.mat”. An Iterative NLS procedure is used in order to cope with the unknown trajectory scenario.

The reconstructed trajectory should emulate a car, which at the beginning is moving at 50 km/h (almost 14 m/s ) along the x-axis.

<img src="C:\Users\carmi\AppData\Roaming\Typora\typora-user-images\image-20201125174429103.png" alt="image-20201125174429103" style="zoom:67%;" />

(predicted velocities in m/s)

##### 5) Kalman filter testing on unknown trajectory

Implement an **Extended Kalman Filter** by using the variable “rhoUEAP.mat” in “Task5_rhoUEAP_GRXX.mat”.An Iterative NLS procedure is used in order to cope with the unknown trajectory scenario.

##### 6) Kalman filter testing on unknown trajectory and with missing measurements

In this task we deal with missing measurements from some/all the 8 APs that make the estimated location less accurate w.r.t. previous tasks.

### Authors

- **Fabio Carminati [Github](https://github.com/fabiocarminati)** *[fabio3.carminati@mail.polimi.it](mailto:fabio3.carminati@mail.polimi.it)*
- **Emanuele Gallone [Github](https://github.com/EmanueleGallone/)** *[emanuele.gallone@mail.polimi.it](mailto:emanuele.gallone@mail.polimi.it)*

### Reference:

In the plots we use a pre-existing function the *plot_dir* 

Kangwon Lee (2020). Plot With Direction (https://www.mathworks.com/matlabcentral/fileexchange/1676-plot-with-direction), MATLAB Central File Exchange.

### License

This project is licensed under the MIT License - see the [LIccccCENSE.md](https://github.com/FabioCarminati/LNSM_Project/blob/master/LICENSE.md) file for details