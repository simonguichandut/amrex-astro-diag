# Convective gradients

This tool constructs the various convective gradients.  It is built
off of the AMReX fgradient tool.  Note: unlike fgradient, we don't
allow different BCs per variable, but instead just use hoextrap.

The three gradients that are evaluated are:
- The true temperature gradient in the model:
$$\nabla=\left(\frac{d\ln T}{d\ln P}\right)$$
- The adiabatic gradient:
$$\nabla_{\rm ad}=\left(\frac{d\ln T}{d\ln P}\right)_S=\frac{\chi_T}{\Gamma_1c_V}\frac{P}{\rho T}$$
- The Ledoux gradient:
$$\nabla_{\rm L}=\nabla_{\rm ad}+B$$
For the composition term $B$, we use MESA's formulation (Paxton et al. 2013, Equation 8). At some grid point $k$ with pressure $P_k=P(\rho_k,T_k,\bm{X}_k)$:
$$B_k=-\frac{1}{\chi_T}\frac{\ln P(\rho_k,T_k,\bm{X}_{k+1})-\ln P_k}{\ln P_{k+1}-\ln P_k}$$
(The top left term is the local pressure evaluated with the composition of the above grid point.)

To build, do:

```
make DIM=2
```

changing the `DIM` line to match the dimension of your plotfile.

It is also important that the network you build with matches
the one used for generating the plotfile.  This is set via
the `NETWORK_DIR` parameter in the `GNUmakefile`.

Runtime parameters are managed by AMReX's ParmParse. To run,
you specify the plotfile via `diag.plotfile`, either in an inputs
file or on the command line, e.g.:

```
./fconvgrad.gnu.ex diag.plotfile=plt00000
```




