# CHiMaD3_dendritic_growth
https://pages.nist.gov/pfhub/benchmarks/benchmark3.ipynb/

https://github.com/eastglow-zz/moose

Solution of ChiMaD benchmark problem #3; dendritic growth using MOOSE-framework

The treatment of anisotropy function is generalized and simplified by newly written kernels associated with Function Parsed system of MOOSE-framework.

--Kernel descriptions--

[AnisotropicGradEnergy] (to be merged with the master branch of MOOSE repository soon)

This kernel calculates the residual and Jacobian of the strong form;

 -L * div( (d/dgrad_aeta)(0.5 * kappa(grad_aeta) * grad_aeta^2) ), where aeta is a scalar field.

 Inputs: L, kappa(grad_aeta), the components of grad_aeta, and the threshold magnitude of gradient of aeta

 where, L is assumed to be a constant, kappa(grad_aeta) is an anisotropy function, can be provided by DerivativeParsedMaterial as a functions of gradient components of aeta, and each component of grad_aeta can be calculated by using VariableGradientComponent AuxKernel from aeta declared as a Variable.
 This kernel does not calculate anisotropy terms where the gradient magnitude of aeta is less than the threshold value.

 To use VariableGradientComponent, each component of grad_aeta, namely, in 2D, dpx and dpy, should be declared as an elemental AuxVariable (e.g. MONOMIAL family), respectively.

 CHiMaD3_DK_w_postprocess.i might be useful as an example for usage of this kernel.


[AnisotropicTimeDerivative] (to be merged with the master branch of MOOSE repository soon)

This kernel calculates the residual and Jacobian of the strong form;

tau(grad_aeta) * (d/dt) aeta, where aeta is a scalar field.

 Inputs: L, kappa(grad_aeta), and the components of grad_aeta,

 where, L is assumed to be a constant, tau(grad_aeta) is an anisotropy function, can be provided by DerivativeParsedMaterial as a functions of gradient components of aeta, and each component of grad_aeta can be calculated by using VariableGradientComponent AuxKernel from aeta declared as a Variable.
 This kernel does not calculate anisotropy terms where the gradient magnitude of aeta is less than the threshold value.

 To use VariableGradientComponent, each component of grad_aeta, namely, in 2D, dpx and dpy, should be declared as an elemental AuxVariable (e.g. MONOMIAL family), respectively.

 CHiMaD3_DK_w_postprocess.i might be useful as an example for usage of this kernel.
 
 Last update: July 18th 2018
