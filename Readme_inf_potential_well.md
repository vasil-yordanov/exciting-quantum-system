# Particle in Infinite potential well

## Eigenfunctions
The eigenfunctions of a particle in a one-dimensional infinite potential well of width L are given by the sine functions that satisfy the boundary conditions (the wavefunction must be zero at the walls of the well). The eigenfunctions are:

```math
\psi_n(x) = \sqrt{\frac{2}{L}} \sin\left(\frac{n \pi x}{L}\right)
```

Here:

- n is a positive integer (1, 2, 3, ...), representing the quantum number of the state.
- L is the width of the potential well.
- x is the position within the well (ranging from 0 to L).

## Energy Levels
The eigenvalues correspond to the quantized energy levels of the particle. They are given by the formula:

```math
E_n = \frac{n^2 \pi^2 \hbar^2}{2mL^2}
```

Here:

- E_n is the energy of the n-th state.
- ħ is the reduced Planck constant.
- m is the mass of the particle.
- L is the width of the potential well.

These eigenfunctions and eigenvalues are derived from solving the time-independent Schrödinger equation for a particle in a one-dimensional box with the given boundary conditions. The solutions demonstrate the principle of quantization in quantum mechanics, where only certain discrete energy levels are allowed for the particle, and the spatial distribution of the particle (described by the wavefunctions) shows characteristic standing wave patterns.

To calculate the matrix element for the transition from the first excited state (n=1) to the third excited state (n=3) in a one-dimensional infinite potential well, we start by defining the eigenfunctions in the well and then proceed to compute the integral that represents the matrix element.

## Eigenfunctions
The eigenfunctions for a particle in a one-dimensional infinite potential well are given by:

```math
\psi_n(x) = \sqrt{\frac{2}{L}} \sin\left(\frac{n \pi x}{L}\right)
```

where n is the quantum number, L is the width of the well, and x is the position within the well.

## Matrix Element Calculation
The matrix element for the transition from the first to the third excited state is defined as:

```math
V_{13}(t) = \int_0^L \psi_3^*(x) V(t,x) \psi_1(x) \, dx
```

Assuming a perturbation potential V(x) and substituting the eigenfunctions:

```math
V_{13}(t) = \int_0^L \sqrt{\frac{2}{L}} \sin\left(\frac{3 \pi x}{L}\right) V(t,x) \sqrt{\frac{2}{L}} \sin\left(\frac{\pi x}{L}\right) \, dx
```

Simplifying, we get:

```math
V_{13}(t) = \frac{2}{L} \int_0^L \sin\left(\frac{3 \pi x}{L}\right) V(t,x) \sin\left(\frac{\pi x}{L}\right) \, dx
```

If V(x) is:

```math
V(t,x)=V_0 \sin\left(\frac{\pi x}{L}\right) \cos\left(\omega t\right)
```

where
```math
\omega = \frac{\Delta E_{13}}{\hbar}
```

Then, we get:

```math
V_{13}(t) = \frac{2}{L} V_0 \cos\left(\omega t\right) \int_0^L \sin\left(\frac{3 \pi x}{L}\right) \sin\left(\frac{\pi x}{L}\right) \sin\left(\frac{\pi x}{L}\right) \, dx = - \frac{8}{15} \frac{V_0}{L}
```
## Rabbi frequency
```math
\Omega = \frac{e E d}{\hbar},
```
where
- E is the electric field amplitude of the applied field
- d is the dipole moment associated with the transition between the two states

The electric dipole moment for a transition between two states in a quantum system is given by the matrix element of the position operator between these states.
If we denote the wavefunctions of the initial and final states as ψ_i and ψ_f, respectively, the dipole moment d is calculated as:

```math
d = \frac{2}{L} \int_0^L \sin\left(\frac{3 \pi x}{L}\right) x \sin\left(\frac{\pi x}{L}\right) \, dx = 0
```
