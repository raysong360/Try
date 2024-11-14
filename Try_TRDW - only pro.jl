# only process side

using ModelingToolkit, MethodOfLines, DifferentialEquations, DomainSets

# Define the desiccant wheel model PDESystem

# One-dimensional honeycombed desiccant wheel (air: 1-dimensional, desiccant: 0-dimensional)

### 0. Independent and state Variables ###
@parameters t z                                    # Independent variables: Time and spatial variables  
@variables ωa(..) Ta(..) ωd(..) Td(..) W(..)    # State variables      # W: moisture content of the desiccant material (kg/kg) 
Dt = Differential(t)
Dz = Differential(z)
@variables Da(..) Ky(..) φw(..) Pws(..) Qst(..)

# @parameters h_const, ρa, A, P, Cpa, Cpv, Cpd, Cpl, Cpm, hv, fd, fm, Patm, p, p0, T0, u, Sh_const    # Dependent variables

Ka = 0.0321           # λ in Zhang 2023; Thermal conductivity of dry air (W/(m·K))
Nu = 2.45
a = 0.0015            # Half duct height (m)
b = 0.0015            # Half duct width (m)
r0 = 0.2              # Radius of desiccant wheel (m)
A = 2 * a * b         # Cross-sectional area (m^2)
P = 2b + 2 * sqrt(b^2 + (a * π)^2) * (3 + (2b / (a * π))^2) / (4 + (2b / (a * π))^2)  # Perimeter (m)
h_const = (Nu * Ka * P) / (4 * A)  # heat transfer coefficient

h_const = h_const                 # Calculated heat transfer coefficient
Sh_const = 2.45                   # Sherwood number as constant
u = 3.0                          # Air velocity on both the process and regeneration side (m/s)
ρa = 1.2                          # Air density (kg/m3)
A = A                        # Cross-sectional area
P = P                      # Perimeter
Cpa = 1009                        # Specific heat of dry air (J/(kg·K))
Cpv = 2028                        # Specific heat of water vapor (J/(kg·K))
Cpd = 921                         # Specific heat of desiccant (J/(kg·K))
Cpl = 4179                       # Specific heat of liquid water (J/(kg·K))
Cpm = 900                         # Specific heat of matrix materials (J/(kg·K))
hv = 2358.0e3                     # Evaporation latent heat of water (J/(kg))
fd = 0.005                        # Specific desiccant mass (kg/m)
fm = 0.003                       # Specific matrix mass (kg/m)
Patm = 1.013e5                    # Atmospheric pressure (Pa)
p = 1.013e5                       # Air pressure (Pa)
p0 = 0.98e5                       # Reference pressure (Pa)
T0 = 256.0


### 1. Governing Equations ###

# Eq 1: Moisture conservation in air
moist_cons_air = Dt(ωa(t,z)) + u * Dz(ωa(t,z)) ~ (Ky(t,z) * P / (ρa * A)) * (ωd(t,z) - ωa(t,z))



# Eq 2: Energy conservation in air
energy_cons_air = Dt(Ta(t,z)) + u * Dz(Ta(t,z)) ~ (h_const * P / (ρa * A * (Cpa + ωa(t,z) * Cpv))) * (Td(t,z) - Ta(t,z)) + (Ky(t,z) * Cpv * P / (ρa * A * (Cpa + ωa(t,z) * Cpv))) * (ωd(t,z) - ωa(t,z)) * (Ta(t,z) - Td(t,z))

# Eq 3: Moisture conservation in desiccant     # ∂W/∂t
moist_cons_des = Dt(W(t,z)) ~ (2 * Ky(t,z) * P / fd) * (ωa(t,z) - ωd(t,z))


# Eq 4: Energy conservation in desiccant

energy_cons_des = Dt(Td(t,z)) ~ (2 * h_const * P / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (Ta(t,z) - Td(t,z)) +
                              (2 * Ky(t,z) * P * Qst(t,z) / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (ωa(t,z) - ωd(t,z)) +
                              (2 * Ky(t,z) * P * Cpv / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (ωa(t,z) - ωd(t,z)) * (Td(t,z) - Ta(t,z))


### 2. Auxiliary Equations ###

# Define Da as a function of Ta, Ky as a function of Da
# @variables Da(t, z) Ky(t, z) 
Da_eq = Da(t,z) ~ 2.302e-5 * p0 / p * (Ta(t,z) / T0)^1.81   # Mass diffusion coefficient (m2 s-1)
Ky_eq = Ky(t,z) ~ ρa * Sh_const * Da(t,z) * P / (4 * A)

# The equilibrium relative humidity on the surface of RD silica gel
φw_W_eq = φw(t,z) ~ 0.0078 - 0.0579 * W(t,z) + 24.1655 * W(t,z)^2 - 124.478 * W(t,z)^3 + 204.226 * W(t,z)^4

# ωd: Humidity ratio of air in equilibrium with the desiccant (kg kg-1)
Pws_Td_eq = Pws(t,z) ~ exp(23.1964 - 3816.44 / max(Td(t,z) - 46.13, 1e-3))
wd_φw_eq = ωd(t,z) ~ 0.622 * φw(t,z) / (Patm / Pws(t,z) - φw(t,z))

# Qst: heat of adsorption [J/kg H2O]
Qst_eq = Qst(t,z) ~ hv * (1.0 + 0.2843 * exp(-10.28 * W(t,z)))


### 3. Construct the PDESystem without parameters, initial, or boundary conditions ###

# Collect all equations
eqs = [moist_cons_air, energy_cons_air, moist_cons_des, energy_cons_des, wd_φw_eq, Pws_Td_eq, φw_W_eq, Da_eq, Ky_eq, Qst_eq]

Ta_pro_in = 273.15 + 34.3
ωa_pro_in = 0.02


# domain
L = 0.2             # Thickness of desiccant wheel (m)
time_cycle = 512.0
time_pro_per = 0.75
time_reg_per = 1 - time_pro_per

domains = [t ∈ IntervalDomain(0.0, time_cycle * time_pro_per),
           z ∈ IntervalDomain(0.0, L)]


# Initial and boundary conditions 
## process
bcs = [Ta(t,0) ~ Ta_pro_in,             # bc: fixed inlet air temperature
       ωa(t,0) ~ ωa_pro_in,                      # bc: fixed inlet humidity ratio  0.01925 at 56.2%
       
       Dz(Td(t,0)) ~ 0.,
       Dz(Td(t,L)) ~ 0.,
       Dz(W(t,0)) ~ 0.,
       Dz(W(t,L)) ~ 0.,

       ωa(0,z) ~ ωa_pro_in,            
       Ta(0,z) ~ Ta_pro_in,   
       W(0,z) ~ 0.001,                   
       ωd(0,z) ~ 0.01e-3,                 
       Td(0,z) ~ 273.15 + 30, 
]


vars = [ωa(t,z), Ta(t,z), W(t,z), Td(t,z), ωd(t,z), φw(t,z), Pws(t,z), Da(t,z), Ky(t,z), Qst(t,z)]




# Define the PDE system with only equations, time, and spatial variables
@named DW_model_pdesys = PDESystem(eqs, bcs, domains, [t,z], vars) # parameters

dz = L/25
discretization = MOLFiniteDifference([z=>dz], t, approx_order=2)


odeprob = discretize(DW_model_pdesys, discretization)


sol = solve(odeprob, Rodas5(); reltol=1e-6, abstol=1e-8, maxiters=10000)

#=
prob = discretize(DW_model_pdesys, discretization)
sys = structural_simplify(DW_model_pdesys)

sol = solve(simplified_prob, Rodas5())  # ; reltol=1e-6, abstol=1e-8, maxiters=10000
=#


sol[ωa(t,z)] .*1000
sol[Td(t,z)] .-273.15
sol[ωd(t,z)] .*1000
sol[W(t,z)]

sol[φw(t,z)]
sol[Pws(t,z)]
sol[Da(t,z)]
sol[Ky(t,z)]
sol[Qst(t,z)]

sol[Ta(t,z)] .-273.15

