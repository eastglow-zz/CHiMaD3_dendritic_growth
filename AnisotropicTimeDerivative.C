//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AnisotropicTimeDerivative.h"

registerMooseObject("PhaseFieldApp", AnisotropicTimeDerivative);

template <>
InputParameters
validParams<AnisotropicTimeDerivative>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addClassDescription("Interfacial-orientation-dependent time derivative Kernel, tau(grad_aeta) * (d/dt) (aeta), where tau(grad_aeta) is anisotropy function in terms of gradient components of aeta");
  params.addParam<MaterialPropertyName>("tau_name","tau_op","Orientation dependent anisotropy function in terms of gradient components of aeta, which can be provided by DerivativeParsedMaterial with derivative_order 2");
  params.addCoupledVar("gradient_component_names", "Name vector of gradient components of aeta, arguments of the tau function, in x y z order, e.g.) in 2D, dpx dpy, in 3D, dpx dpy dpz, where dpx, dpy, and dpz are defined as AuxVariables(FIRST order, MONOMIAL family) and calculated by VariableGradientComponent AuxKernel (execute_on = LINEAR)");
  params.addParam<Real>("gradmag_threshold",1e-7,"Threshold value to turn on anisotropy term calculations; grad_mag > thres ? anisotropic calc. : isotropic calc.");
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  return params;
}

AnisotropicTimeDerivative::AnisotropicTimeDerivative(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<TimeKernel>>(parameters),
  _nvar(_coupled_moose_vars.size()),
  _tau(getMaterialProperty<Real>("tau_name")),
  _gradmag_threshold(getParam<Real>("gradmag_threshold")),
  _dtau_darg(_nvar)
{
  /// Get derivative data
  for (unsigned int i = 0; i < _nvar; ++i)
  {
    MooseVariable * ivar = _coupled_standard_moose_vars[i];
    const VariableName iname = ivar->name();
    if (iname == _var.name())
      paramError("gradient_component_names",\
                 "The kernel variable should not be specified in the coupled `gradient_components` parameter.");

    /// The 1st derivatives
    _dtau_darg[i] = &getMaterialPropertyDerivative<Real>("tau_name", iname);
  }
}

void
AnisotropicTimeDerivative::initialSetup()
{
  validateCoupling<Real>("tau_name");
}

Real
AnisotropicTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * _u_dot[_qp] * _tau[_qp];
}

Real
AnisotropicTimeDerivative::computeQpJacobian()
{
  if (_nvar > 0 && _grad_u[_qp] * _grad_u[_qp] > _gradmag_threshold * _gradmag_threshold)
  {
    Real dtau_du = 0;
    for (unsigned int i = 0; i < _nvar; i++)
    {
      dtau_du += (*_dtau_darg[i])[_qp] * _grad_phi[_j][_qp](i);
    }
    return (dtau_du * _u_dot[_qp] + _du_dot_du[_qp] * _tau[_qp]) * _test[_i][_qp];
  }else{
    return (_du_dot_du[_qp] * _tau[_qp]) * _test[_i][_qp];
  }
}

Real
AnisotropicTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (0)
  {
    const unsigned int cvar = mapJvarToCvar(jvar);
    return _test[_i][_qp] * _u_dot[_qp] * (*_dtau_darg[cvar])[_qp] * _phi[_j][_qp];
  }else{
    return 0;
  }
}

void
AnisotropicTimeDerivative::computeJacobian()
{
  if (_lumping)
  {
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());

    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
          ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
  }
  else
    TimeKernel::computeJacobian();
}
