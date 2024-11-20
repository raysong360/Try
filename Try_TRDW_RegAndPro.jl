using ModelingToolkit, MethodOfLines, DifferentialEquations, DomainSets, Interpolations
using Plots

# Define the desiccant wheel model PDESystem

# One-dimensional honeycombed desiccant wheel (air: 1-dimensional, desiccant: 0-dimensional)

### 0. Independent and state Variables ###
@parameters t z                                    # Independent variables: Time and spatial variables  
@variables ωa(..) Ta(..) ωd(..) Td(..) W(..)    # State variables      # W: moisture content of the desiccant material (kg/kg) 
Dt = Differential(t)
Dz = Differential(z)
@variables Da(..) Ky(..) φw(..) Pws(..) Qst(..)


Ka = 0.0321           # λ in Zhang 2023; Thermal conductivity of dry air (W/(m·K))
Nu = 2.45
a = 0.0015            # Half duct height (m)
b = 0.0015            # Half duct width (m)
r0 = 0.2              # Radius of desiccant wheel (m)
A = 2 * a * b         # Cross-sectional area (m^2)
P = 2b + 2 * sqrt(b^2 + (a * π)^2) * (3 + (2b / (a * π))^2) / (4 + (2b / (a * π))^2)  # Perimeter (m)
h = (Nu * Ka * P) / (4 * A)  # heat transfer coefficient
Sh = 2.45                   # Sherwood number as constant
ρa = 1.2                          # Air density (kg/m3)
Cpa = 1009                        # Specific heat of dry air (J/(kg·K))
Cpv = 2028                        # Specific heat of water vapor (J/(kg·K))
Cpd = 921                         # Specific heat of desiccant (J/(kg·K))
Cpl = 4179                       # Specific heat of liquid water (J/(kg·K))
Cpm = 900                         # Specific heat of matrix materials (J/(kg·K))
Lv = 2358.0e3                     # Evaporation latent heat of water (J/(kg))
fd = 0.005                        # Specific desiccant mass (kg/m)
fm = 0.003                       # Specific matrix mass (kg/m)
Patm = 1.013e5                    # Atmospheric pressure (Pa)
p = 1.013e5                       # Air pressure (Pa)
p0 = 0.98e5                       # Reference pressure (Pa)
T0 = 256.0


@parameters u    # Dependent variables
pars = [u]
pars_value = Dict(u => 3.0)



### 1. Governing Equations ###

# Eq 1: Moisture conservation in air
moist_cons_air = Dt(ωa(t,z)) + u * Dz(ωa(t,z)) ~ (Ky(t,z) * P / (ρa * A)) * (ωd(t,z) - ωa(t,z))



# Eq 2: Energy conservation in air
energy_cons_air = Dt(Ta(t,z)) + u * Dz(Ta(t,z)) ~ (h * P / (ρa * A * (Cpa + ωa(t,z) * Cpv))) * (Td(t,z) - Ta(t,z)) + (Ky(t,z) * Cpv * P / (ρa * A * (Cpa + ωa(t,z) * Cpv))) * (ωd(t,z) - ωa(t,z)) * (Ta(t,z) - Td(t,z))

# Eq 3: Moisture conservation in desiccant     # ∂W/∂t
moist_cons_des = Dt(W(t,z)) ~ (2 * Ky(t,z) * P / fd) * (ωa(t,z) - ωd(t,z))


# Eq 4: Energy conservation in desiccant

energy_cons_des = Dt(Td(t,z)) ~ (2 * h * P / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (Ta(t,z) - Td(t,z)) +
                              (2 * Ky(t,z) * P * Qst(t,z) / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (ωa(t,z) - ωd(t,z)) +
                              (2 * Ky(t,z) * P * Cpv / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (ωa(t,z) - ωd(t,z)) * (Td(t,z) - Ta(t,z))


### 2. Auxiliary Equations ###

# Define Da as a function of Ta, Ky as a function of Da
# @variables Da(t, z) Ky(t, z) 
Da_eq = Da(t,z) ~ 2.302e-5 * p0 / p * (Ta(t,z) / T0)^1.81   # Mass diffusion coefficient (m2 s-1)
Ky_eq = Ky(t,z) ~ ρa * Sh * Da(t,z) * P / (4 * A)

# The equilibrium relative humidity on the surface of RD silica gel
φw_W_eq = φw(t,z) ~ 0.0078 - 0.0579 * W(t,z) + 24.1655 * W(t,z)^2 - 124.478 * W(t,z)^3 + 204.226 * W(t,z)^4

# ωd: Humidity ratio of air in equilibrium with the desiccant (kg kg-1)
Pws_Td_eq = Pws(t,z) ~ exp(23.1964 - 3816.44 / max(Td(t,z) - 46.13, 1e-3))
wd_φw_eq = ωd(t,z) ~ 0.622 * φw(t,z) / (Patm / Pws(t,z) - φw(t,z))

# Qst: heat of adsorption [J/kg H2O]
Qst_eq = Qst(t,z) ~ Lv * (1.0 + 0.2843 * exp(-10.28 * W(t,z)))


### 3. Construct the PDESystem without parameters, initial, or boundary conditions ###

# Collect all equations
eqs = [moist_cons_air, energy_cons_air, moist_cons_des, energy_cons_des, wd_φw_eq, Pws_Td_eq, φw_W_eq, Da_eq, Ky_eq, Qst_eq]


Ta_pro_in = 273.15 + 34.3
ωa_pro_in = 0.02
Ta_reg_in = 273.15 + 100
ωa_reg_in = 0.005


L = 0.2             # Thickness of desiccant wheel (m)
time_cycle = 512.0
time_pro_per = 0.75
time_reg_per = 1 - time_pro_per

n_seg = 20
dz = L/n_seg

z_points = range(0.0, stop=L, length=n_seg+1)


function generate_bcs(time_cycle = 512.0, time_reg_per = 0.25; W_ini, ωa_ini, Ta_ini, ωd_ini, Td_ini, n_seg, L, whichphase)

    dz = L/n_seg
    z_points = range(0.0, stop=L, length=n_seg+1)

    if length(W_ini) != length(z_points) || length(ωa_ini) != length(z_points) ||
        length(Ta_ini) != length(z_points) || length(ωd_ini) != length(z_points) ||
        length(Td_ini) != length(z_points)
         throw(ArgumentError("The length of initial value arrays must match the length of z_points."))
     end
    
    initial_values = Dict(
        :W => W_ini,
        :ωa => ωa_ini,
        :Ta => Ta_ini,
        :ωd => ωd_ini,
        :Td => Td_ini
    )

    interpolators = Dict()
    functions = Dict()

    for (key, values) in initial_values
        let key = key, values = values
            # Print debugging information
            println("Processing key: $key")
            println("z_points: ", z_points)
            println("Values: ", values)
            println("Boundary check: ", values[1], " (z=0.0), ", values[end], " (z=L)")

            # Linear interpolation
            local itp
            itp = interpolate((z_points,), values, Gridded(Linear()))

            local func_name
            func_name = Symbol("sitp_$(key)_function")
            
            @eval begin
                $(func_name)(z) = $itp(z)
                @register_symbolic $(func_name)(z)
            end

            interpolators[key] = itp
            functions[key] = eval(func_name)
        end
    end

    if whichphase == "Reg"
        return [
            Ta(t,0) ~ Ta_reg_in,
            ωa(t,0) ~ ωa_reg_in,
            Dz(Td(t,0)) ~ 0.0,
            Dz(Td(t,L)) ~ 0.0,
            Dz(W(t,0)) ~ 0.0,
            Dz(W(t,L)) ~ 0.0,

            ωa(0,z) ~ functions[:ωa](z),
            Ta(0,z) ~ functions[:Ta](z),
            W(0,z) ~ functions[:W](z),
            ωd(0,z) ~ functions[:ωd](z),
            Td(0,z) ~ functions[:Td](z)
        ]
    elseif whichphase == "Pro"
        return [
            Ta(t,0) ~ Ta_pro_in,
            ωa(t,0) ~ ωa_pro_in,
            Dz(Td(t,0)) ~ 0.0,
            Dz(Td(t,L)) ~ 0.0,
            Dz(W(t,0)) ~ 0.0,
            Dz(W(t,L)) ~ 0.0,

            ωa(time_cycle * time_reg_per, z) ~ functions[:ωa](z),
            Ta(time_cycle * time_reg_per,z) ~ functions[:Ta](z),
            W(time_cycle * time_reg_per,z) ~ functions[:W](z),
            ωd(time_cycle * time_reg_per,z) ~ functions[:ωd](z), 
            Td(time_cycle * time_reg_per,z) ~ functions[:Td](z)
        ]
    else
        throw(ArgumentError("Invalid phase type: $whichphase. Expected 'Reg' or 'Pro'."))
    end
end



function ini_vect_0(n_seg)
    (
        W_ini_vect=[0.1 for _ in 1:n_seg+1], 
        ωa_ini_vect=[0.02 for _ in 1:n_seg+1], 
        Ta_ini_vect=[273.15 + 40. for _ in 1:n_seg+1], 
        ωd_ini_vect=[0.015 for _ in 1:n_seg+1], 
        Td_ini_vect=[273.15 + 40. for _ in 1:n_seg+1]
    )
end


# 1 pde
bcs_reg = generate_bcs(
    W_ini = ini_vect_0(20).W_ini_vect, 
    ωa_ini = ini_vect_0(20).ωa_ini_vect, 
    Ta_ini = ini_vect_0(20).Ta_ini_vect, 
    ωd_ini = ini_vect_0(20).ωd_ini_vect, 
    Td_ini = ini_vect_0(20).Td_ini_vect, 
    n_seg=20, L=0.2, whichphase="Reg"
    )

vars = [ωa(t,z), Ta(t,z), W(t,z), Td(t,z), ωd(t,z), φw(t,z), Pws(t,z), Da(t,z), Ky(t,z), Qst(t,z)]


domains_reg = [t ∈ IntervalDomain(0.0, 128.0),
           z ∈ IntervalDomain(0.0, L)]

@named reg_pdesys = PDESystem(eqs, bcs_reg, domains_reg, [t,z], vars, pars; defaults=pars_value) 
discretization = MOLFiniteDifference([z=>L/n_seg], t, approx_order=2)
odeprob_reg = discretize(reg_pdesys, discretization)

# 注意pars的顺序：[u, Ka, Nu, a, b, r0, A, P, h, Sh, ρa, Cpa, Cpv, Cpd, Cpl, Cpm, Lv, fd, fm, Patm, p, p0, T0]
# odeprob = remake(odeprob, p = [3])
sol_reg = solve(odeprob_reg, Rodas5())


sol_reg[ωa(t,z)] .*1000
sol_reg[Td(t,z)] .-273.15
sol_reg[ωd(t,z)] .*1000
sol_reg[W(t,z)]

sol_reg[φw(t,z)]
sol_reg[Pws(t,z)]
sol_reg[Da(t,z)]
sol_reg[Ky(t,z)]
sol_reg[Qst(t,z)]

sol_reg[Ta(t,z)] .-273.15

sol_reg[W(t,z)][end, :]


bcs_pro = generate_bcs(
    W_ini = reverse(sol_reg[W(t,z)][end, :]), 
    ωa_ini = reverse(sol_reg[ωa(t,z)][end, :]), 
    Ta_ini = reverse(sol_reg[Ta(t,z)][end, :]), 
    ωd_ini = reverse(sol_reg[ωd(t,z)][end, :]), 
    Td_ini = reverse(sol_reg[Td(t,z)][end, :]), 
    n_seg=20, L=0.2, whichphase="Pro"
    )

domains_pro = [t ∈ IntervalDomain(time_cycle * time_reg_per, time_cycle),
            z ∈ IntervalDomain(0.0, L)]
@named pro_pdesys = PDESystem(eqs, bcs_pro, domains_pro, [t,z], vars, pars; defaults=pars_value) 
discretization = MOLFiniteDifference([z=>L/n_seg], t, approx_order=2)
odeprob_pro = discretize(pro_pdesys, discretization)
sol_pro = solve(odeprob_pro, Rodas5())

sol_pro[Ta(t,z)] .-273.15

t_values_reg = sol_reg.t
t_values_pro = sol_pro.t

Ta_reg_out = sol_reg[Ta(t,z)][:, end] .-273.15
Ta_pro_out = sol_pro[Ta(t,z)][:, end] .-273.15

plot(t_values_reg, Ta_reg_out, label="Reg Outlet T",
     xlabel="Time (s)", ylabel="Temperature (°C)", color=:red)
hline!([Ta_reg_in - 273.15], label="Reg Inlet T", linestyle=:dash, color=:red)
plot!(t_values_pro, Ta_pro_out, label="Pro Outlet T",
     xlabel="Time (s)", color=:purple)
hline!([Ta_pro_in - 273.15], label="Pro Inlet T", linestyle=:dash, color=:purple)




ωa_reg_out = sol_reg[ωa(t,z)][:, end] .*1000
ωa_pro_out = sol_pro[ωa(t,z)][:, end] .*1000

plot(t_values_reg, ωa_reg_out, label="Reg Outlet ω",
     xlabel="Time (s)", ylabel="Humidity ratio (g/kg)", color=:red)
hline!([ωa_reg_in] .*1000, label="Reg Inlet ω", linestyle=:dash, color=:red)
plot!(t_values_pro, ωa_pro_out, label="Pro Outlet ω",
     xlabel="Time (s)", color=:purple)
hline!([ωa_pro_in] .*1000, label="Pro Inlet ω", linestyle=:dash, color=:purple)
