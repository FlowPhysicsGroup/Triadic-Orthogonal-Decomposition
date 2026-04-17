# Triadic orthogonal decomposition in MATLAB
![Graphical abstract](aux/TOD_widelogo.png "TOD modes")
TOD() is a MATLAB implementation of triadic orthogonal decomposition (TOD), a method tailored to triadic interactions among three wave components that enable interscale energy transfer, developed by Yeung, Chu & Schmidt [[1](https://doi.org/10.1017/jfm.2026.11183)]. For applications in fluid dynamics, TOD distinguishes three components - a momentum recipient, donor, and catalyst - and recovers laws governing triad conservation. It identifies coherent flow structures optimally capturing spectral momentum transfer, quantifies their coupling and energy exchange in an energy-budget bispectrum, and reveals the regions where they interact.

The direct numerical simulation data provided along with this example is a $Re=100$ cylinder wake calculated using the polyharmonic splines with polynomial augmentation radial basis function-finite differences implementation of the fractional-step, staggered-grid incompressible Navier-Stokes solver developed by Chu & Schmidt [[2](https://doi.org/10.1016/j.jcp.2022.111756)]. A discussion of the results in these examples can be found in [[1](https://doi.org/10.1017/jfm.2026.11183)].

## Files
| File        |     Description     |
| ------------- |:-------------|
| tod.m | Triadic orthogonal decomposition in MATLAB |
| example_1.m | Load data, build differential operators, and calculate and plot TOD |
| example_2.m | Apply TOD to 3D data |
| example_3.m | Apply TOD to datasets that do not include velocities |
| aux/bluered.m | Blue-white-red colormap |
| aux/Dmats_SBP.m | Finite-difference operators implemented by Aaron Towne |
| aux/drawWavyLine.m | Plot a wavy line |
| aux/findnearest.m | Find nearest element in a list |
| aux/plotBispectrum.m | Plot TOD mode bispectrum or modal energy budget |
| aux/plotModes2D.m | Plot 2D TOD modes and visualize them in Feynman-type diagrams |
| aux/plotModes3D.m | Plot 3D TOD modes as isosurfaces and visualize them in Feynman-type diagrams |
| aux/plotGeneralObservable.m | Plot TOD modes computed from non-velocity data |
| aux/cylinder_Re100.mat | $Re=100$ cylinder wake database |
| aux/graphical_abstract.jpg | Graphical abstract |
| LICENSE.txt | License |

## Usage
```
           fn
            ^
   fk_      |________
    |\     /|       |
       \ /  |       |
       / \  |       |
     /     \|       |
----+-------+-------+-> fl
    |       |\     / 
    |       |  \ /   
    |       |  / \  
    |_______|/     \
            |
```
Figure. Sketch of the $f_l$ - $f_n$ plane. Triads are expressed as frequency
        triplets, $\{f_k=f_{n-l},f_l,f_n\}$, or frequency index triplets,
        $(k=n-l,l,n)$. The principal region consists of $f_n\ge0$ for real
        data, all triads otherwise.

`[L,P,F] = TOD(X)` returns the triadic orthogonal decomposition of the
data matrix X. The first dimension of X must be time, the second
dimension variable indices. X can have any number of additional spatial
dimensions. The mode bispectrum is returned in L and the modes in cell
array P. The spatial dimensions of the modes are identical to those of
X. The convective modes are stored in P{l,n}(1,...) and the recipient
modes in P{l,n}(2,...). The cell indices {l,n} of P are the donor and
recipient indices. F is the frequency vector. If DT is not specified,
the frequency index is returned in F. Although TOD(X) automatically
chooses default spectral estimation parameters, it is recommended to
manually specify problem-dependent parameters on a case-to-case basis.

`[L,P,F] = TOD(X,WINDOW)` uses a temporal window. If WINDOW is a
vector, X is divided into segments of the same length as WINDOW. Each
segment is then weighted (pointwise multiplied) by WINDOW. If WINDOW is
a scalar, a Hamming window of length WINDOW is used. If WINDOW is
omitted or empty, a Hamming window is used.

`[L,P,F] = TOD(X,WINDOW,WEIGHT)` uses a spatial inner product weight,
usually quadrature weights. WEIGHT must have the same spatial dimensions
as X.

`[L,P,F] = TOD(X,WINDOW,WEIGHT,NOVERLAP)` increases the number of
segments by overlapping consecutive blocks by NOVERLAP snapshots.
NOVERLAP defaults to 50% of the length of WINDOW if not specified.

`[L,P,F] = TOD(X,WINDOW,WEIGHT,NOVERLAP,DT)` uses the time step DT
between consecutive snapshots to determine a physical frequency F.

`[L,P,F] = TOD(X,WINDOW,WEIGHT,NOVERLAP,DT,OPTS)` specifies options:  
- `OPTS.Q`           function that takes inputs (q1,q2) and returns the
                 nonlinear term of the equation [function |
                 {@(q1,q2) permute(q1.\*q2,[2 1 3])}]  
- `OPTS.LHS`         function that takes input (q) and returns the left hand
                 side of the equation [function | {@(q) q}]  
- `OPTS.nmode`       number of modes (ranks) per triad to store [integer |
                 {nBlks}]  
- `OPTS.threads`     number of threads if MATLAB supports parallel
                 functionality; if not, or if OPTS.threads<=1, the
                 code runs serially [integer | {maxNumCompThreads}]  
- `OPTS.isreal`      reality of data, which determines regions of the
                 bispectrum to compute, see figure caption above
                 [logical | {isreal(X)}]  
- `OPTS.precision`   compute TOD in single precision ['single' | 'double' |
                 {class(X)}]  
- `OPTS.mean`        provide a mean that is subtracted from each
                 snapshot [array of size X | {0}]  
- `OPTS.nfreq`       restrict computation to |l|,|k|,|n|<=OPTS.nfreq
                 [integer | {all}]  

`[L,P,F,T] = TOD(...)` returns the modal energy budget in T. T has the
same dimensions as L.

`[L,P,F,T,A] = TOD(...)` returns the expansion coefficients in A. The
matrices of coefficients associated with the convective and recipient
modes are stored respectively in A(1,l,n,:,:) and A(2,l,n,:,:).

`[L,P,F,T,A,Xi] = TOD(...)` returns the donor modes in Xi{l,n}(1,...) and
catalyst modes in Xi{l,n}(2,...), where l and n are the same indices that
index into L, P, T, and A.

## References
[[1](https://doi.org/10.1017/jfm.2026.11183)] Yeung, B., Chu, T., and Schmidt, O. T., Triadic orthogonal decomposition reveals nonlinearity in fluid flows, J. Fluid Mech. 1031, A34, 2026  
DOI 10.1017/jfm.2026.11183

[[2](https://doi.org/10.1016/j.jcp.2022.111756)] Chu, T. and Schmidt, O. T., RBF-FD discretization of the Navier-Stokes equations on scattered but staggered nodes, J. Comp. Phys. 474, 111756, 2023  
DOI 10.1016/j.jcp.2022.111756
