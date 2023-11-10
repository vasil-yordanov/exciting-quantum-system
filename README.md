# Exciting 1D quantum system in infinite potential well from level n=1 to n=3.

![Exciting a quantum system from level 1 to level 3 in infinite potential well](schrodinger12T-6.gif).

# Time-Dependent 1D Schrödinger Equation

The time-dependent 1D Schrödinger equation is given by:

$$
i\hbar \frac{\partial \Psi(x, t)}{\partial t} = -\frac{\hbar^2}{2m} \frac{\partial^2 \Psi(x, t)}{\partial x^2} + V(x, t) \Psi(x, t)
$$

where:
- Ψ(x, t) is the wave function of the particle.
- ħ is the reduced Planck constant.
- m is the mass of the particle.
- V(x, t) is the potential energy as a function of position x and time t.

In a one-dimensional infinite potential well (also known as a particle in a box), the eigenfunctions and eigenvalues are fundamental concepts that describe the quantized states of a particle confined in a well with infinitely high walls. The potential well is typically defined as having zero potential energy inside the well and infinite potential energy outside.

1. **Eigenfunctions**: 
   The eigenfunctions of a particle in a one-dimensional infinite potential well of width `L` are given by the sine functions that satisfy the boundary conditions (the wavefunction must be zero at the walls of the well). The eigenfunctions are:

   ```math
   \psi_n(x) = \sqrt{\frac{2}{L}} \sin\left(\frac{n \pi x}{L}\right)

Here:

n is a positive integer (1, 2, 3, ...), representing the quantum number of the state.
L is the width of the potential well.
x is the position within the well (ranging from 0 to L).

1. **Eigenvalues (Energy Levels):
   The eigenvalues correspond to the quantized energy levels of the particle. They are given by the formula:

   ```math
E_n = \frac{n^2 \pi^2 \hbar^2}{2mL^2}

Here:

E_n is the energy of the n-th state.
\hbar is the reduced Planck constant.
m is the mass of the particle.
L is the width of the potential well.

These eigenfunctions and eigenvalues are derived from solving the time-independent Schrödinger equation for a particle in a one-dimensional box with the given boundary conditions. The solutions demonstrate the principle of quantization in quantum mechanics, where only certain discrete energy levels are allowed for the particle, and the spatial distribution of the particle (described by the wavefunctions) shows characteristic standing wave patterns.

To calculate the matrix element for the transition from the first excited state (n=1) to the third excited state (n=3) in a one-dimensional infinite potential well, we start by defining the eigenfunctions in the well and then proceed to compute the integral that represents the matrix element.

1. **Eigenfunctions**:
   The eigenfunctions for a particle in a one-dimensional infinite potential well are given by:

   ```math
   \psi_n(x) = \sqrt{\frac{2}{L}} \sin\left(\frac{n \pi x}{L}\right)
   ```
   
where n is the quantum number, L is the width of the well, and x is the position within the well.

** Matrix Element Calculation:
The matrix element V 13 for the transition from the first to the third excited state is defined as:

   ```math
   V_{13} = \int_0^L \psi_3^*(x) V(x) \psi_1(x) \, dx
   ```

Assuming a perturbation potential V(x) and substituting the eigenfunctions:

```math
V_{13} = \int_0^L \sqrt{\frac{2}{L}} \sin\left(\frac{3 \pi x}{L}\right) V(x) \sqrt{\frac{2}{L}} \sin\left(\frac{\pi x}{L}\right) \, dx
```

Simplifying, we get:

```math
V_{13} = \frac{2}{L} \int_0^L \sin\left(\frac{3 \pi x}{L}\right) V(x) \sin\left(\frac{\pi x}{L}\right) \, dx
```

If V(x) is:

```math
V(x)=V_0 sin(\pi x/L)
```

Then, we get:

```math
V_{13} = \frac{2}{L} \int_0^L \sin\left(\frac{3 \pi x}{L}\right) V_0 sin(\pi x/L) \sin\left(\frac{\pi x}{L}\right) \, dx = - \frac{8}{15} \frac{V_0}{L}
```
