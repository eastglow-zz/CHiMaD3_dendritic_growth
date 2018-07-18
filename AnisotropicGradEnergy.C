//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html



#include "AnisotropicGradEnergy.h"

registerMooseObject("PhaseFieldApp", AnisotropicGradEnergy);

template <>
InputParameters
validParams<AnisotropicGradEnergy>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates a general form of -L * div((d/dgrad_aeta) (0.5 * kappa(grad_aeta) * (grad_aeta)^2))");
  params.addParam<MaterialPropertyName>("mob_name","L","The mobility used with the kernel, assumed as a constant");
  params.addParam<MaterialPropertyName>("kappa_name","kappa_op","Orientation dependent anisotropy function in terms of gradient components of aeta, which can be provided by DerivativeParsedMaterial with derivative_order 2");
  params.addCoupledVar("gradient_component_names","Name vector of gradient components of aeta, arguments of the kappa function, in x y z order, e.g.) in 2D, dpx dpy, in 3D, dpx dpy dpz, where dpx, dpy, and dpz are defined as AuxVariables(FIRST order, MONOMIAL family) and calculated by VariableGradientComponent AuxKernel (execute_on = LINEAR)");
  params.addParam<Real>("gradmag_threshold",1e-7,"Threshold value to turn on anisotropy term calculations; grad_mag > thres ? anisotropic calc. : isotropic calc.");
  return params;
}

AnisotropicGradEnergy::AnisotropicGradEnergy(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
  _nvar(_coupled_moose_vars.size()),
  _L(getMaterialProperty<Real>("mob_name")),
  _kappa(getMaterialProperty<Real>("kappa_name")),
  _gradmag_threshold(getParam<Real>("gradmag_threshold")),
  _dkappa_darg(_nvar),
  _d2kappa_darg2(_nvar)
{
  /// Get derivative data
  for (unsigned int i = 0; i < _nvar; ++i)
  {
    MooseVariable * ivar = _coupled_standard_moose_vars[i];
    const VariableName iname = ivar->name();
    if (iname == _var.name())
      paramError("gradient_component_names",\
                 "The kernel variable should not be specified in the coupled `gradient_component_names` parameter.");

    /// The 1st derivatives
    _dkappa_darg[i] = &getMaterialPropertyDerivative<Real>("kappa_name", iname);

    /// The 2nd derivatives
    _d2kappa_darg2[i].resize(_nvar);
    for (unsigned int j = 0; j < _nvar; ++j)
    {
      const VariableName jname = _coupled_moose_vars[j]->name();
      if (jname == _var.name())
        paramError("gradient_component_names",\
                   "The kernel variable should not be specified in the coupled `gradient_component_names` parameter.");
      _d2kappa_darg2[i][j] = &getMaterialPropertyDerivative<Real>("kappa_name", iname, jname);
    }
  }
}

void
AnisotropicGradEnergy::initialSetup()
{
  validateCoupling<Real>("kappa_name");
}

RealGradient
AnisotropicGradEnergy::get_dkappa_darg(unsigned int qp) // This function must be called in computeQp* functions
{
  RealGradient v0(0.0, 0.0, 0.0);
  switch (_nvar) {
    case 1:
      {
        RealGradient v1((*_dkappa_darg[0])[qp], 0.0, 0.0);
        return v1;
      }
      break;
    case 2:
      {
        RealGradient v2((*_dkappa_darg[0])[qp], (*_dkappa_darg[1])[qp], 0.0);
        return v2;
      }
      break;
    case 3:
      {
        RealGradient v3((*_dkappa_darg[0])[qp], (*_dkappa_darg[1])[qp], (*_dkappa_darg[2])[qp]);
        return v3;
      }
      break;
    default:
      return v0;
  }
}

RealGradient
AnisotropicGradEnergy::get_d2kappa_darg2(unsigned int i, unsigned int qp)  // This function must be called in computeQp* functions
{
  RealGradient v0(0.0, 0.0, 0.0);
  switch (_nvar) {
    case 1:
      {
        RealGradient v1((*_d2kappa_darg2[i][0])[qp], 0.0, 0.0);
        return v1;
      }
      break;
    case 2:
      {
        RealGradient v2((*_d2kappa_darg2[i][0])[qp], (*_d2kappa_darg2[i][1])[qp], 0.0);
        return v2;
        break;
      }
      break;
    case 3:
      {
        RealGradient v3((*_d2kappa_darg2[i][0])[qp], (*_d2kappa_darg2[i][1])[qp], (*_d2kappa_darg2[i][2])[qp]);
        return v3;
      }
      break;
    default:
      return v0;
  }
}

RealGradient
AnisotropicGradEnergy::get_dargv_darg(unsigned int i)
{
  RealGradient v0(0.0, 0.0, 0.0);
  switch (i) {
    case 0:
      {
        RealGradient v1(1.0, 0.0, 0.0);
        return v1;
      }
      break;
    case 1:
      {
        RealGradient v2(0.0, 1.0, 0.0);
        return v2;
        break;
      }
      break;
    case 2:
      {
        RealGradient v3(0.0, 0.0, 1.0);
        return v3;
      }
      break;
    default:
      return v0;
  }
}

Real
AnisotropicGradEnergy::computeQpResidual()
{
  Real grad_u_sq = _grad_u[_qp] * _grad_u[_qp];
  Real grad_u_dot_grad_test = _grad_u[_qp] * _grad_test[_i][_qp];
  if (_nvar > 0 && grad_u_sq > _gradmag_threshold * _gradmag_threshold)
  {
    RealGradient dkappa_dgradaeta = get_dkappa_darg(_qp);
    Real dkappa_dgradaeta_dot_grad_test = dkappa_dgradaeta * _grad_test[_i][_qp];
    return _L[_qp] * (_kappa[_qp] * grad_u_dot_grad_test\
                      + 0.5 * grad_u_sq * dkappa_dgradaeta_dot_grad_test);
  }else{
    return _L[_qp] * _kappa[_qp] * grad_u_dot_grad_test;
  }
}

Real
AnisotropicGradEnergy::computeQpJacobian()
{

  Real grad_u_sq = _grad_u[_qp] * _grad_u[_qp];
  Real kappa_gradphi_dot_grad_test = _kappa[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  if (_nvar > 0 && grad_u_sq > _gradmag_threshold * _gradmag_threshold)
  {
    RealGradient dkappa_dgradaeta = get_dkappa_darg(_qp);
    Real dkappa_dgradaeta_dot_grad_phi = dkappa_dgradaeta * _grad_phi[_j][_qp];
    Real grad_u_dot_grad_test = _grad_u[_qp] * _grad_test[_i][_qp];
    Real grad_u_dot_grad_phi = _grad_u[_qp] * _grad_phi[_j][_qp];
    Real dkappa_dgradaeta_dot_grad_test = dkappa_dgradaeta * _grad_test[_i][_qp];
    Real d2kappa_dgradaeta2_dot_grad_phi_dot_grad_test = 0.0;
    for (unsigned int i = 0; i < _nvar; ++i)
    {
      for (unsigned int j = 0; j < _nvar; ++j)
      {
        d2kappa_dgradaeta2_dot_grad_phi_dot_grad_test += \
          _grad_test[_i][_qp](i) * (*_d2kappa_darg2[i][j])[_qp] * _grad_phi[_j][_qp](j);
      }
    }
    return _L[_qp] * (kappa_gradphi_dot_grad_test \
                     + dkappa_dgradaeta_dot_grad_phi * grad_u_dot_grad_test \
                     + grad_u_dot_grad_phi * dkappa_dgradaeta_dot_grad_test \
                     + 0.5 * grad_u_sq * d2kappa_dgradaeta2_dot_grad_phi_dot_grad_test);
  }else{
    return _L[_qp] * kappa_gradphi_dot_grad_test;
  }
}


Real
AnisotropicGradEnergy::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (0)
  {
    // get the coupled variable jvar is referring to
    const unsigned int cvar = mapJvarToCvar(jvar);
    Real grad_u_sq = _grad_u[_qp] * _grad_u[_qp];
    Real xcvar_dot_grad_test = get_dargv_darg(cvar) * _grad_test[_i][_qp];
    Real grad_u_dot_grad_test = _grad_u[_qp] * _grad_test[_i][_qp];
    Real dkappa_dgradaeta_dot_grad_test = get_dkappa_darg(_qp) * _grad_test[_i][_qp];
    Real d2kappa_dgradaeta_dcvar_dot_grad_test = get_d2kappa_darg2(cvar, _qp) * _grad_test[_i][_qp];
    return _L[_qp] * _phi[_j][_qp] * ((*_dkappa_darg[cvar])[_qp] * grad_u_dot_grad_test \
                      + _kappa[_qp] * xcvar_dot_grad_test \
                      + _grad_u[_qp](cvar) * dkappa_dgradaeta_dot_grad_test \
                      + 0.5 * grad_u_sq * d2kappa_dgradaeta_dcvar_dot_grad_test);
  }else{
    //return _L[_qp] * _phi[_j][_qp] * _kappa[_qp] * xcvar_dot_grad_test;
    return 0;
  }
}
