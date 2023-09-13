# PESMD-FES-RATES Workflow: Single Force Molecular Dynamics
Workflow runs single force Molecular Dynamics simulations of a two-particle system having their interaction described by a simple potential energy function.
The potential function can be 1D with two wells representins a bound and unbound state or 2D with a third well representing an itermediate state.
In the 1D case a force is applied to the interparticle distance as to pull particles apart and the same distance is biased. 
In the 2D case both components of distance can be biased while pulling on one of the components.
The biasing is done with the metadynamics method and FE surfaces and kinetic rates as a function of pulling force can be obtained for each case.
