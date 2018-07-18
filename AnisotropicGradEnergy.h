//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ANISOTROPICGRADENERGY_H
#define ANISOTROPICGRADENERGY_H

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

class AnisotropicGradEnergy;

template <>
InputParameters validParams<AnisotropicGradEnergy>();

class AnisotropicGradEnergy : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  AnisotropicGradEnergy(const InputParameters & parameters);
  virtual void initialSetup();

protected:
  RealGradient get_dkappa_darg(unsigned int qp);
  RealGradient get_d2kappa_darg2(unsigned int cvar, unsigned int qp);
  RealGradient get_dargv_darg(unsigned int cvar);
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const unsigned int _nvar;
  /// Phase-field mobility (assumed to be constant)
  const MaterialProperty<Real> & _L;
  /// Gradient energy coefficient; 0th order derivative from Material data
  const MaterialProperty<Real> & _kappa;
  const Real _gradmag_threshold;
  /// Gradient energy coefficient; 1st order derivatives from Material data
  std::vector<const MaterialProperty<Real> *> _dkappa_darg;
  /// Gradient energy coefficient; 2nd order derivatives from Material data
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2kappa_darg2;
};

#endif // ANISOTROPICGRADENERGY_H
