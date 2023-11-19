# Exciting particle in infinite potential well from level n=1 to n=3.

![Exciting a quantum system from level 1 to level 3 in infinite potential well](schrodinger12T-6.gif).

# Exciting H-atom from n=1 to n=2.

![Exciting H atom from n=1 to n=2](h-atom.gif).

# Time-Dependent 1D Schrödinger Equation

The time-dependent 1D Schrödinger equation is given by:

```math
i\hbar \frac{\partial \Psi(t, x)}{\partial t} = -\frac{\hbar^2}{2m} \frac{\partial^2 \Psi(t, x)}{\partial x^2} + V(t, x) \Psi(t, x)
```

where:
- Ψ(t, x) is the wave function of the particle.
- ħ is the reduced Planck constant.
- m is the mass of the particle.
- V(t, x) is the potential energy as a function of position x and time t.

In a one-dimensional infinite potential well (also known as a particle in a box), the eigenfunctions and eigenvalues are fundamental concepts that describe the quantized states of a particle confined in a well with infinitely high walls. The potential well is typically defined as having zero potential energy inside the well and infinite potential energy outside.

# Infinite potential well
```math
V(x) = 
    \begin{cases}
        V_0 \sin(\frac{\pi x}{L}) .* cos(\omega t), & 0 < x < L \\
        \infty, & \text{otherwise}
    \end{cases}
```

[Details about Infinite potential well](./Readme_inf_potential_well.md).


# Coulomb potential
```math
V(x) = \frac{V_0}{r} \cos(\omega t), 0 < t < 50 \frac{h}{\Delta E_{12}}
```
[Details about H-atom](./Readme_h_atom.md).

