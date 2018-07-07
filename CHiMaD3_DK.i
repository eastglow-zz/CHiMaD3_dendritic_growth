
#
# CHiMaD benchmark problem #3, dendritic growth
# https://pages.nist.gov/pfhub/benchmarks/benchmark3.ipynb/
# Dong-Uk Kim
# 06/20/2018

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 960
  ny = 960
  xmax = 960
  ymax = 960
[]

[Variables]
  [./p] #phase-field variable
    order = FIRST
    family = LAGRANGE
  [../]
  [./u] #dimensionless temperature variable
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./dpx]
    order = FIRST
    family = MONOMIAL
  [../]
  [./dpy]
    order = FIRST
    family = MONOMIAL
  [../]

#  [./FEdensity]
#    order = FIRST
#    family = MONOMIAL
#  [../]
[]

[ICs]
  [./p_circleIC]
    type = SmoothCircleIC
    variable = p
    x1 = 0
    y1 = 0
    radius = 8
    int_width = 1
    invalue = 1.0
    outvalue = -1.0
  [../]
  [./u_constantIC]
    type = ConstantIC
    variable = u
    value = -0.3
  [../]
[]

[BCs]
  #Do nothing means no-flux boundary condition
[]

[AuxKernels]
  [./get_dpx]
    type = VariableGradientComponent
    variable = dpx
    gradient_variable = p
    component = x
    execute_on = LINEAR
  [../]
  [./get_dpy]
    type = VariableGradientComponent
    variable = dpy
    gradient_variable = p
    component = y
    execute_on = LINEAR
  [../]
[]

[Kernels]
  #Phase-field
  [./time_derivative_p]
    type = AnisotropicTimeDerivative
    variable = p
    tau_name = tau_aniso
    gradient_component_names = 'dpx dpy'
  [../]
  [./InterfacialE]
    type = AnisotropicGradEnergy
    variable = p
    mob_name = L
    kappa_name = Wsq_aniso
    gradient_component_names = 'dpx dpy'
  [../]
  [./localFE]
    type = AllenCahn
    variable = p
    f_name = f_doublewell
    mob_name = L
  [../]
  [./DrivingForce]
    type = ACSwitching
    variable = p
    Fj_names = 'lambda_U'
    hj_names = 'h'
    args = 'p u'
    mob_name = L
  [../]

  #Dimensionless temperature field
  [./time_derivative_u]
    type = TimeDerivative
    variable = u
  [../]
  [./div_grad_mu]
    type = ACInterface
    variable = u
    mob_name = L
    kappa_name = D
    variable_L = false
  [../]
  [./Partition_term]
    type = CoupledSwitchingTimeDerivative
    variable = u
    v = p
    Fj_names = 'oneovertwo'
    hj_names = 'switching_U'
    args =     'p'
  [../]
[]

[Materials]
  [./f_doublewell]
    type = DerivativeParsedMaterial
    f_name = f_doublewell
    args = 'p'
    function = '-0.5*p^2+0.25*p^4'
    derivative_order = 2
    #outputs = exodus
  [../]

  [./lambda_U]
    type = DerivativeParsedMaterial
    f_name = lambda_U
    material_property_names = 'lambda'
    args = 'u'
    function = 'lambda*u'
  [../]

  [./h]
    type = DerivativeParsedMaterial
    f_name = h
    args = 'p'
    function = 'p*(1-2*p^2/3+p^4/5)'
    derivative_order = 2
    #outputs = exodus
  [../]

  [./switching_U]
    type = DerivativeParsedMaterial
    f_name = switching_U
    args = 'p'
    function = 'p'
    derivative_order = 1
    #outputs = exodus
  [../]

  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'tau0 W0 D  eps0 L oneovertwo'
    prop_values = '1    1  10 0.05 1 -0.5'
  [../]

  [./lambda]
    type = ParsedMaterial
    f_name = lambda
    material_property_names = 'D tau0 W0'
    function = 'D*tau0/0.6267/W0/W0'
  [../]

  [./Anisotropic_Wsquare]
    type = DerivativeParsedMaterial
    f_name = Wsq_aniso
    material_property_names = 'W0 eps4'
    args = 'dpx dpy'
    function = 'if(eps4 > 0, W0^2 * (1 + eps4 * (dpx^4 + dpy^4 - 6*dpx^2*dpy^2)/(dpx^2 + dpy^2)^2)^2, W0^2)'
    derivative_order = 2
    #outputs = exodus
  [../]

  [./Anisotropic_tau]
    type = DerivativeParsedMaterial
    f_name = tau_aniso
    material_property_names = 'tau0 eps4'
    args = 'dpx dpy'
    function = 'if(eps4 > 0, tau0 * (1 + eps4 * (dpx^4 + dpy^4 - 6*dpx^2*dpy^2)/(dpx^2 + dpy^2)^2)^2, tau0)'
    derivative_order = 2
    #outputs = exodus
  [../]

  #To prevent divided by zero
  [./eps4]
    type = ParsedMaterial
    f_name = eps4
    material_property_names = 'eps0'
    args = 'p'
    function = 'if((0.9-p)*(0.9+p) >= 0, eps0, 0)'
    #outputs = exodus
  [../]

  # For the postprocess
#  [./solid_volume]
#    type = ParsedMaterial
#    f_name = solid_volume_per_element
#    args = 'p'
#    function = '(1 + p)/2'
#  [../]
#  [./element_volume]
#    type = ParsedMaterial
#    f_name = volume_per_element
#    function = '1'
#  [../]
#  [./FEdensity_material]
#    type = ParsedMaterial
#    f_name = f_density
#    material_property_names = 'Wsq_aniso f_doublewell lambda_U h'
#    args = 'p u dpx dpy'
#    function = '0.5*Wsq_aniso*(dpx^2+dpy^2 + f_doublewell + lambda_U*h)'
#    #outputs = exodus
#  [../]
[]

[Preconditioning]
  #active = ''
  [./cw_coupling]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  scheme = bdf2

  #petsc_options_iname = '-pc_type -sub_pc_type'
  #petsc_options_value = 'asm      lu'
  petsc_options_iname = '-pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      1'

  l_max_its = 20
  l_tol = 1e-4
  nl_max_its = 15
  nl_rel_tol = 1e-8

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 2
    cutback_factor = 0.8
    optimal_iterations = 20
    iteration_window = 10
  [../]
  end_time = 1500
[]

[Outputs]
  exodus = true
  #csv = true
  print_perf_log = true
[]

#[Debug]
#  show_var_residual_norms = true
#[]

#[Postprocessors]
#  [./Total_solid_volume]
#    type = ElementIntegralMaterialProperty
#    mat_prop = solid_volume_per_element
#  [../]
#  #[./Total_volume]
#  #  type = ElementIntegralMaterialProperty
#  #  mat_prop = volume_per_element
#  #[../]
#  [./Total_FE]
#    type = ElementIntegralMaterialProperty
#    mat_prop = f_density
#  [../]
#  [./Interface_location_along_x_axis]
#    type = FindValueOnLine
#    start_point = '0 0 0'
#    end_point = '960 0 0'
#    target = 0
#    depth = 13
#    tol = 1e-1
#    v = p
#  [../]
#[]
