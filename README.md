# EZ-NMA
Normal Mode Analysis with Elastic Network Model


### Principles of Normal Mode Analysis
_(please refresh the page if equations are not shown properly)_
In elastic network model, the potential is expressed as:

<a href="https://www.codecogs.com/eqnedit.php?latex=E=\frac{1}{2}\sum_{d_{ij}<R_c}K(d_{ij}-d_{ij}^0)^2\quad(1)" target="_blank"><img src="https://latex.codecogs.com/png.latex?E=\frac{1}{2}\sum_{d_{ij}<R_c}K(d_{ij}-d_{ij}^0)^2\quad(1)" title="E=\frac{1}{2}\sum_{d_{ij}<R_c}K(d_{ij}-d_{ij}^0)^2\quad(1)" /></a>

Moreover, the potential energy can be expanded as:

<a href="https://www.codecogs.com/eqnedit.php?latex=V(\mathbf{r})=V(\mathbf{r}^0)&plus;\sum_i&space;\big(&space;\frac{\partial&space;V}{\partial&space;r_i}&space;\big)^0(r_i-r_i^0)&space;&plus;\frac{1}{2}\sum_i\sum_j&space;\big(&space;\frac{\partial^2&space;V}{\partial&space;r_i\partial&space;r_j}&space;\big)^0(r_i-r_i^0)(r_j-r_j^0)&plus;\dots\quad(2)" target="_blank"><img src="https://latex.codecogs.com/png.latex?V(\mathbf{r})=V(\mathbf{r}^0)&plus;\sum_i&space;\big(&space;\frac{\partial&space;V}{\partial&space;r_i}&space;\big)^0(r_i-r_i^0)&space;&plus;\frac{1}{2}\sum_i\sum_j&space;\big(&space;\frac{\partial^2&space;V}{\partial&space;r_i\partial&space;r_j}&space;\big)^0(r_i-r_i^0)(r_j-r_j^0)&plus;\dots\quad(2)" title="V(\mathbf{r})=V(\mathbf{r}^0)+\sum_i \big( \frac{\partial V}{\partial r_i} \big)^0(r_i-r_i^0) +\frac{1}{2}\sum_i\sum_j \big( \frac{\partial^2 V}{\partial r_i\partial r_j} \big)^0(r_i-r_i^0)(r_j-r_j^0)+\dots\quad(2)" /></a>

Here the system is assumed to be near equilibrium, thus the first order derivative is zero. The ground state energy can also be assigned as zero, thus only the second order derivative is left:

<a href="https://www.codecogs.com/eqnedit.php?latex=V(\mathbf{q})=\frac{1}{2}\mathbf{\Delta&space;q^TH\Delta&space;q},~~~~\mathrm{where}~H_{ij}=\left(\frac{\partial^2V}{\partial&space;q_i\partial&space;q_j}&space;\right)^0\quad(3)" target="_blank"><img src="https://latex.codecogs.com/png.latex?V(\mathbf{q})=\frac{1}{2}\mathbf{\Delta&space;q^TH\Delta&space;q},~~~~\mathrm{where}~H_{ij}=\left(\frac{\partial^2V}{\partial&space;q_i\partial&space;q_j}&space;\right)^0\quad(3)" title="V(\mathbf{q})=\frac{1}{2}\mathbf{\Delta q^TH\Delta q},~~~~\mathrm{where}~H_{ij}=\left(\frac{\partial^2V}{\partial q_i\partial q_j} \right)^0\quad(3)" /></a>

According to Newton's Law II, we have:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{M}\frac{\mathrm{d^2}\mathbf{\Delta&space;q}}{\mathrm{d}t^2}&plus;\mathbf{H\Delta&space;q}=0\quad(4)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\mathbf{M}\frac{\mathrm{d^2}\mathbf{\Delta&space;q}}{\mathrm{d}t^2}&plus;\mathbf{H\Delta&space;q}=0\quad(4)" title="\mathbf{M}\frac{\mathrm{d^2}\mathbf{\Delta q}}{\mathrm{d}t^2}+\mathbf{H\Delta q}=0\quad(4)" /></a>

The general solution of this differential equations has the following form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{u}_k(t)=\mathbf{a}_ke^{-i\omega_kt}\quad(5)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\mathbf{u}_k(t)=\mathbf{a}_ke^{-i\omega_kt}\quad(5)" title="\mathbf{u}_k(t)=\mathbf{a}_ke^{-i\omega_kt}\quad(5)" /></a>

Plug this solution back to the equation, we have an eigenvalue problem. By matrix diagonalization, we can obtain the eigenvalue and eigenvectors (normal modes) of the sysmte:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Hu}_k=\omega^2_k\mathbf{Mu}_k\xrightarrow{mass-weighted}&space;\mathbf{\tilde{H}\tilde{U}=\tilde{U}\Lambda}&space;\Longleftrightarrow\color{red}\mathbf{\tilde{U}^T\tilde{H}\tilde{U}=\Lambda}\color{black}\quad(6)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\mathbf{Hu}_k=\omega^2_k\mathbf{Mu}_k\xrightarrow{mass-weighted}&space;\mathbf{\tilde{H}\tilde{U}=\tilde{U}\Lambda}&space;\Longleftrightarrow\color{red}\mathbf{\tilde{U}^T\tilde{H}\tilde{U}=\Lambda}\color{black}\quad(6)" title="\mathbf{Hu}_k=\omega^2_k\mathbf{Mu}_k\xrightarrow{mass-weighted} \mathbf{\tilde{H}\tilde{U}=\tilde{U}\Lambda} \Longleftrightarrow\color{red}\mathbf{\tilde{U}^T\tilde{H}\tilde{U}=\Lambda}\color{black}\quad(6)" /></a>

### Notes
This is the CPU version. Parallel version with new OpenMP features and GPU version are under testing and debugging...
A short demo of spastin is provided:

<img src="https://github.com/wangqi1990uc/ez-nma/blob/master/nma-demo.gif" width="40%" height="40%" />

