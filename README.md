# adaptive-koopman-mpc
This repository contains supplementary material to the paper "Adaptive Koopman Model Predictive Control of Simple Serial Robots" (submitted to IROS 2025). 

### Real system experiments  
To test the adaptive Koopman Model Predictive Control (KMPC) algorithm, trajectory tracking experiments were conducted on a 1R and a belt-driven 2R robot. The following video shows an experiment, where force disturbance rejection of the controller was tested on the 2R system.
![](2R_experiments.gif)



#### Controller parameters
The following table shows MPC weights used in real system experiments with 1R and 2R robot systems. 

| Controller | System | $Q = diag(.)$               |$Q_f = diag(.)$               | $R  = diag(.)$| Buffer size|
| :----------|:------ |:----------------------------|:---------------------------- |:--------------|:-----------|
| lin.   MPC | 2R     | $[2.5, 1.1, 0.01, 0.01]$    | $[5.5, 2.5, 0.01, 0.01]$     | $[0.30, 0.30]$|
| adapt. KMPC| 2R     | $[3.5, 3.0, 0, \ldots, 0]$  | $[3.5, 3.0, 0, \ldots, 0]$   | $[0.02, 0.02]$| 
| stat.  KMPC| 2R     | $[10.5, 8.45, 0, \ldots, 0]$| $[10.5, 8.45, 0, \ldots, 0]$ | $[0.15, 0.15]$|
| lin.   MPC | 1R     | $[1.1, 0.01]$               | $[5.5, 0.01]$                | $[0.1]$       |
| adapt. KMPC| 1R     | $[3.0, 0, \ldots, 0]$       | $[3.0, 0, \ldots, 0]$        | $[1.0]$       | 
| stat.  KMPC| 1R     | $[8.45, 0, \ldots, 0]$      | $[8.45, 0, \ldots, 0]$       | $[0.4]$       |

The weights determine how much emphasis is placed in tracking of intermediate states $Q$ and final state $Q_f$, and applied controls (R). They correspond to the unaugmented optimization problem, as denoted in equation (6a) in the paper. After state augmentation and rewriting the optimization problem in condensed form, the overall weight matrices are 

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

