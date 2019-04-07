# EZ-NMA
Normal Mode Analysis with Elastic Network Model


### Principles of Normal Mode Analysis
In elastic network model, the potential is expressed as:

![equation 1](demo/eqn/eqn1.png)

Moreover, the potential energy can be expanded as:

![equation 2](demo/eqn/eqn2.png)

Here the system is assumed to be near equilibrium, thus the first order derivative is zero. The ground state energy can also be assigned as zero, thus only the second order derivative is left:

![equation 3](demo/eqn/eqn3.png)

According to Newton's Law II, we have:

![equation 4](demo/eqn/eqn4.png)

The general solution of this differential equations has the following form:

![equation 5](demo/eqn/eqn5.png)

Plug this solution back to the equation, we have an eigenvalue problem. By matrix diagonalization, we can obtain the eigenvalue and eigenvectors (normal modes) of the system:

![equation 6](demo/eqn/eqn6.png)

### Notes
- Both CPU and GPU implementations are preliminary and under testing and debugging. 

- A short demo of spastin is provided:

![nma demo](demo/eqn/nma-demo.gif)

