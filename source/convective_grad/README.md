# Convective gradients

This tool constructs the various convective gradients.  It is built
off of the AMReX fgradient tool.  Note: unlike fgradient, we don't
allow different BCs per variable, but instead just use hoextrap.

To build, do:

```
make DIM=2
```

changing the `DIM` line to match the dimension of your plotfile.

It is also important that the network you build with matches
the one used for generating the plotfile.  This is set via
the `NETWORK_DIR` parameter in the `GNUmakefile`.

Runtime parameters are managed by AMReX's ParmParse.



