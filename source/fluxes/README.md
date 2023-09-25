# EOS demo

This tool evaluates the following fluxes:
- The true convective heat flux 
$$F_{\rm conv}=\rho c_p v \delta T$$
- The *mixing-length theory* convective heat flux
$$F_{\rm conv,MLT}=\frac{\rho c_pT}{QgH_P}v^3$$
- The kinetic energy flux
$$F_{\rm kin}=\rho v^3$$
- The radiative flux
$$F_{\rm rad}=-\frac{4ac T^3}{3\kappa\rho}\frac{dT}{dr}$$
- The chemical fluxes of species $i$
$$F_{\rm X_i}=\rho v \delta X_i$$

The first three fluxes have units of $\rm{erg/s/cm}^2$, while the composition fluxes have units of $\rm{g/s/cm}^2$

For the first flux, we use $\delta T=T-\bar{T}$, which exists as a code variable in MAESTROeX. We use the same definition for $\delta X_i$, however we have to manually compute the horizontal average.

In the second flux, note that $Q=\left(\frac{d\ln\rho}{d\ln T}\right)_P=\frac{\chi_T}{\chi_\rho}$ (neglecting composition gradients) comes from the EOS. Also note that it is assumed that the mixing-length $\ell$ is equal to the pressure scale height $H_P=P/\rho g$.


To build, do:

```
make DIM=2
```

changing the `DIM` line to match the dimension of your plotfile.

It is also important that the network you build with matches
the one used for generating the plotfile.  This is set via
the `NETWORK_DIR` parameter in the `GNUmakefile`.

Runtime parameters are managed by AMReX's ParmParse.  To run,
you specify the plotfile via `diag.plotfile`, either in an inputs
file or on the command line, e.g.:

```
./fluxes.gnu.ex diag.plotfile=plt00000
```



