# EZ-NMA
Normal Mode Analysis with Elastic Network Model


### Principles of Normal Mode Analysis
In elastic network model, the potential is expressed as:

![Elastic network model](<img src="http://www.sciweavers.org/tex2img.php?eq=E%3D%5Cfrac%7B1%7D%7B2%7D%5Csum_%7Bd_%7Bij%7D%3CR_c%7DK%28d_%7Bij%7D-d_%7Bij%7D%5E0%29%5E2%5Cquad%281%29&bc=White&fc=Black&im=png&fs=12&ff=mathpazo&edit=0" align="center" border="0" alt="E=\frac{1}{2}\sum_{d_{ij}<R_c}K(d_{ij}-d_{ij}^0)^2\quad(1)" width="219" height="49" ></R_c>)

Moreover, the potential energy can be expanded as:

![Taylor expansion](<img src="http://www.sciweavers.org/tex2img.php?eq=V%28%5Cmathbf%7Br%7D%29%3DV%28%5Cmathbf%7Br%7D%5E0%29%2B%5Csum_i%20%5Cbig%28%20%5Cfrac%7B%5Cpartial%20V%7D%7B%5Cpartial%20r_i%7D%20%5Cbig%29%5E0%28r_i-r_i%5E0%29%20%2B%5Cfrac%7B1%7D%7B2%7D%5Csum_i%5Csum_j%20%5Cbig%28%20%5Cfrac%7B%5Cpartial%5E2%20V%7D%7B%5Cpartial%20r_i%5Cpartial%20r_j%7D%20%5Cbig%29%5E0%28r_i-r_i%5E0%29%28r_j-r_j%5E0%29%2B%5Cdots%5Cquad%282%29&bc=White&fc=Black&im=png&fs=12&ff=mathpazo&edit=0" align="center" border="0" alt="V(\mathbf{r})=V(\mathbf{r}^0)+\sum_i \big( \frac{\partial V}{\partial r_i} \big)^0(r_i-r_i^0) +\frac{1}{2}\sum_i\sum_j \big( \frac{\partial^2 V}{\partial r_i\partial r_j} \big)^0(r_i-r_i^0)(r_j-r_j^0)+\dots\quad(2)" width="593" height="50" />)

Here the system is assumed to be near equilibrium, thus the first order derivative is zero. The ground state energy can also be assigned as zero, thus only the second order derivative is left:

![Proximity to equilibrium assumption
](<img src="http://www.sciweavers.org/tex2img.php?eq=V%28%5Cmathbf%7Bq%7D%29%3D%5Cfrac%7B1%7D%7B2%7D%5Cmathbf%7B%5CDelta%20q%5ETH%5CDelta%20q%7D%2C%7E%7E%7E%7E%5Cmathrm%7Bwhere%7D%7EH_%7Bij%7D%3D%5Cleft%28%5Cfrac%7B%5Cpartial%5E2V%7D%7B%5Cpartial%20q_i%5Cpartial%20q_j%7D%20%5Cright%29%5E0%5Cquad%283%29&bc=White&fc=Black&im=png&fs=12&ff=mathpazo&edit=0" align="center" border="0" alt="V(\mathbf{q})=\frac{1}{2}\mathbf{\Delta q^TH\Delta q},~~~~\mathrm{where}~H_{ij}=\left(\frac{\partial^2V}{\partial q_i\partial q_j} \right)^0\quad(3)" width="383" height="56" />)

According to Newton's Law II, we have:

![Newton’s law II](<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cmathbf%7BM%7D%5Cfrac%7B%5Cmathrm%7Bd%5E2%7D%5Cmathbf%7B%5CDelta%20q%7D%7D%7B%5Cmathrm%7Bd%7Dt%5E2%7D%2B%5Cmathbf%7BH%5CDelta%20q%7D%3D0%5Cquad%283%29&bc=White&fc=Black&im=png&fs=12&ff=mathpazo&edit=0" align="center" border="0" alt="\mathbf{M}\frac{\mathrm{d^2}\mathbf{\Delta q}}{\mathrm{d}t^2}+\mathbf{H\Delta q}=0\quad(4)" width="189" height="40" />)

The general solution of this differential equations has the following form:
![solution of Newton's law II](<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cmathbf%7Bu%7D_k%28t%29%3D%5Cmathbf%7Ba%7D_ke%5E%7B-i%5Comega_kt%7D%5Cquad%285%29&bc=White&fc=Black&im=png&fs=12&ff=mathpazo&edit=0" align="center" border="0" alt="\mathbf{u}_k(t)=\mathbf{a}_ke^{-i\omega_kt}\quad(5)" width="158" height="22" />)

Plug this solution back to the equation, we have an eigenvalue problem. By matrix diagonalization, we can obtain the eigenvalue and eigenvectors (normal modes) of the sysmte:

![Diagonalization of Hessian](<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cmathbf%7BHu%7D_k%3D%5Comega%5E2_k%5Cmathbf%7BMu%7D_k%5Cxrightarrow%7Bmass-weighted%7D%20%5Cmathbf%7B%5Ctilde%7BH%7D%5Ctilde%7BU%7D%3D%5Ctilde%7BU%7D%5CLambda%7D%20%5CLongleftrightarrow%5Ccolor%7Bred%7D%5Cmathbf%7B%5Ctilde%7BU%7D%5ET%5Ctilde%7BH%7D%5Ctilde%7BU%7D%3D%5CLambda%7D%5Ccolor%7Bblack%7D%5Cquad%286%29&bc=White&fc=Black&im=png&fs=12&ff=mathpazo&edit=0" align="center" border="0" alt="\mathbf{Hu}_k=\omega^2_k\mathbf{Mu}_k\xrightarrow{mass-weighted} \mathbf{\tilde{H}\tilde{U}=\tilde{U}\Lambda} \Longleftrightarrow\color{red}\mathbf{\tilde{U}^T\tilde{H}\tilde{U}=\Lambda}\color{black}\quad(6)" width="460" height="28" />)

### Notes
This is the CPU version. Parallel version with new OpenMP features and GPU version are under testing and debugging...
A short demo of spastin is provided:

<img src="https://github.com/wangqi1990uc/ez-nma/blob/master/nma-demo.gif" width="40%" height="40%" />

