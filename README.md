# PESMD-FES-RATES Workflow: Single Force Molecular Dynamics
Workflow runs single force Molecular Dynamics simulations of a two-particle system having their interaction described by a simple potential energy function.
The potential function can be 1D with two wells representins a bound and unbound state or 2D with a third well representing an itermediate state.
In the 1D case a force is applied to the interparticle distance as to pull particles apart and the same distance is biased. 
In the 2D case both components of distance can be biased while pulling on one of the components.
The biasing is done with the metadynamics method and FE surfaces and kinetic rates as a function of pulling force can be obtained for each case.

## Contents

+ **`./thumb`:** contains thumbnail for workflow.
+ **`./parsl_utils`:** contains scripts necessay to set up local and remote environments.
+ **`./plumed_inputs`:** contains a general input files  which will be edited and copied to corresponding running directories.
+ **`./utils`:** contains analysis code.
+ **`./requirements`:** contains yaml files specifying packages needed by this workflow.  These dependencies are installed automatically by `parsl_utils`.
+ **`parsl_wrapper.sh`:** bash script to run scripts in parsl_utils.
+ **`main.py`:** is the the workflow code that submits and manages jobs. Analysis is done after all jobs are finished
+ **`workflow.xml`:** is the file that takes user input in xml format for the workflow