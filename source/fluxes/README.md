# EOS demo

This tool evaluates the following fluxes:

- The true convective heat flux 
$$F_{\rm conv}=\rho c_p v \delta T$$
This is currently set up to use $\delta T=T-\bar{T}$, which exists as a code variable in MAESTROeX.

- The *mixing-length theory* (MLT) convective heat flux
$$F_{\rm MLT}=\rho c_p\sqrt{g\delta H_p}(\nabla-\nabla_{\rm ad})$$
where $\delta=\left(\frac{d\ln\rho}{d\ln T}\right)_P=\frac{\chi_T}{\chi_\rho}$ comes from the EOS. This expression assumes efficient convection ($\nabla_{\rm element}\approx\nabla_{\rm ad}$), ignores composition gradients, assumes the mixing length $\ell$ is equal to the pressure scale height $H_p=p/\rho g$, and ignores other order unity parameters from various MLT formulations.

- The MLT flux can also be expressed as a function of the velocity, under the same assumptions as above:
$$F_{\rm MLT,v}=\frac{\rho c_pT}{\delta gH_p}|v|^3$$


- The kinetic energy flux
$$F_{\rm kin}=\rho v^3$$
This should be the same as the MLT flux, up to some constant scaling, except for not using the absolute value of the velocity.

- The radiative flux
$$F_{\rm rad}=-\frac{4ac T^3}{3\kappa\rho}\frac{dT}{dr}$$
This requires that you include `CONDUCTIVITY_DIR` in the makefile

- The chemical fluxes of species $i$
$$F_{\rm X_i}=\rho v X_i$$
Currently, only the hydrogen flux is included in the code, but it is straightforward to include other species.

The first three fluxes have units of $\rm{erg/s/cm}^2$, while the composition fluxes have units of $\rm{g/s/cm}^2$

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



