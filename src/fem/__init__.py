"""Lobule-scale finite elements for LiPNet.

The network model in src/model/ represents a whole lobule as one edge of conductance
w = r^4 / l. This package opens that edge up: anisotropic Darcy perfusion of a hexagonal
lobule, transport of oxygen (with its haemoglobin-bound pool), ethanol, acetaldehyde and a
pericentrally cleared drug, a nitric-oxide wall sub-model, and the two quantities that go
back to the graph, the geometric shape factor kappa and the flow-dependent resistance load
z_R(z_F).

    params     physical constants, in SI, with the source of each value
    mesh       structured triangular mesh of the hexagon minus the central vein
    fem_core   P1 assembly, conservative SUPG advection-diffusion-reaction, Newton
    darcy      perfusion, lobule resistance, transit time, wall shear stress
    transport  O2, ethanol, acetaldehyde, drug, metabolic zonation
    two_phase  sinusoid and hepatocyte as two continua with a finite membrane
    no_wall    NO around a single sinusoid (it is a boundary layer, not a lobule field)
    coupling   shape factor, shear-NO autoregulation, patch into a LiPNet graph
    figures    one function per figure, results/fem_*.pdf and .png
    run_fem    command line

See docs/FEM.md.
"""
