# adaptive-koopman-mpc
This repository contains supplementary material to the paper "Adaptive Koopman Model Predictive Control of Simple Serial Robots" (submitted to IROS 2025). 

### Real system experiments  
To test the adaptive Koopman Model Predictive Control (KMPC) algorithm, trajectory tracking experiments were conducted on a 1R and a belt-driven 2R robot. The following video shows an experiment on the 2R system, which demonstrates the controllers ability to reject force disturbances. 

![](2R_experiments.gif)


#### Controller parameters
The table below shows the MPC weights used in real system experiments for the 1R and 2R robot systems.  

| Controller | System | $Q = diag(.)$               |$Q_f = diag(.)$               | $R  = diag(.)$| Buffer size  |
| :----------|:------ |:----------------------------|:---------------------------- |:--------------|:-------------|
| lin.   MPC | 2R     | $[2.5, 1.1, 0.01, 0.01]$    | $[5.5, 2.5, 0.01, 0.01]$     | $[0.30, 0.30]$| -            |
| adapt. KMPC| 2R     | $[3.5, 3.0, 0, \ldots, 0]$  | $[3.5, 3.0, 0, \ldots, 0]$   | $[0.02, 0.02]$| 1500         | 
| stat.  KMPC| 2R     | $[10.5, 8.45, 0, \ldots, 0]$| $[10.5, 8.45, 0, \ldots, 0]$ | $[0.15, 0.15]$| 1500         |
| lin.   MPC | 1R     | $[1.1, 0.01]$               | $[5.5, 0.01]$                | $[0.1]$       | -            |
| adapt. KMPC| 1R     | $[3.0, 0, \ldots, 0]$       | $[3.0, 0, \ldots, 0]$        | $[1.0]$       | 1000         |
| stat.  KMPC| 1R     | $[8.45, 0, \ldots, 0]$      | $[8.45, 0, \ldots, 0]$       | $[0.4]$       | 1000         |

These weights determine how much emphasis is placed on tracking of intermediate states ($Q$) and final state ($Q_f$), and applied controls ($R$). They correspond to the unaugmented optimization problem, as denoted in equation (6a) in the paper. After state augmentation and rewriting the optimization problem in condensed form, the overall weight matrices are: 

$$\begin{align*}
\mathbf{Q} =   
    \begin{bmatrix} 
    Q & 0 & \cdots & 0\\ 
    0 & Q & \ddots  & 0\\ 
    \vdots & \ddots & \ddots & \vdots\\     
    0 & 0 & \cdots & Q_f\\         
\end{bmatrix},
\mathbf{R} =   
    \begin{bmatrix} 
    R & 0 & \cdots & 0\\ 
    0 & R & \ddots  & 0\\ 
    \vdots & \ddots & \ddots & \vdots\\     
    0 & 0 & \cdots & R\\         
\end{bmatrix}.
\end{align*}$$

The buffer size relates to the number of data points stored within the circular buffer, which is used for the extended dynamic mode decomposition in the Koopman based controllers. 


### Simulation code 
In the paper, experimental comparison between adaptive KMPC, linearization MPC and a further Koopman based controller, denoted static KMPC, is done. This repository contains code implementations of adaptive KMPC and comparative controllers, which allows to simulate tracking control and mimic the hardware experiments under simplified conditions. Variable names in the source files are consistent with the notation used in the paper, such that equations can be easily identified. 
Exemplary usage is shown in two Jupyter notebooks, which simulate reference tracking control of [1R system](src/reference_tracking_single_pendulum.ipynb) and [2R system](src/reference_tracking_double_pendulum.ipynb), along with a comparison of the controllers.   


