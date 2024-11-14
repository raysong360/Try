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


# Governing Equations
moist_cons_air = Dt(ωa(t,z)) + u * Dz(ωa(t,z)) ~ (Ky(t,z) * P / (ρa * A)) * (ωd(t,z) - ωa(t,z))
energy_cons_air = Dt(Ta(t,z)) + u * Dz(Ta(t,z)) ~ (h_const * P / (ρa * A * (Cpa + ωa(t,z) * Cpv))) * (Td(t,z) - Ta(t,z)) +
                                             (Ky(t,z) * Cpv * P / (ρa * A * (Cpa + ωa(t,z) * Cpv))) * (ωd(t,z) - ωa(t,z)) * (Ta(t,z) - Td(t,z))
moist_cons_des = Dt(W(t,z)) ~ (2 * Ky(t,z) * P / fd) * (ωa(t,z) - ωd(t,z))
energy_cons_des = Dt(Td(t,z)) ~ (2 * h_const * P / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (Ta(t,z) - Td(t,z)) +
                             (2 * Ky(t,z) * P * Qst(t,z) / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (ωa(t,z) - ωd(t,z)) +
                             (2 * Ky(t,z) * P * Cpv / (fd * (Cpd + W(t,z) * Cpl) + fm * Cpm)) * (ωa(t,z) - ωd(t,z)) * (Td(t,z) - Ta(t,z))

# Auxiliary Equations
Da_eq = Da(t,z) ~ 2.302e-5 * p0 / p * (Ta(t,z) / T0)^1.81
Ky_eq = Ky(t,z) ~ ρa * Sh_const * Da(t,z) * P / (4 * A)
φw_W_eq = φw(t,z) ~ 0.0078 - 0.0579 * W(t,z) + 24.1655 * W(t,z)^2 - 124.478 * W(t,z)^3 + 204.226 * W(t,z)^4
Pws_Td_eq = Pws(t,z) ~ exp(23.1964 - 3816.44 / max(Td(t,z) - 46.13, 1e-3))
wd_φw_eq = ωd(t,z) ~ 0.622 * φw(t,z) / (Patm / Pws(t,z) - φw(t,z))
Qst_eq = Qst(t,z) ~ hv * (1.0 + 0.2843 * exp(-10.28 * W(t,z)))

# Equation system
eqs = [moist_cons_air, energy_cons_air, moist_cons_des, energy_cons_des, wd_φw_eq, Pws_Td_eq, φw_W_eq, Da_eq, Ky_eq, Qst_eq]

# Boundary conditions for process and regen phases
Ta_in_pro, ωa_in_pro = 273.15 + 34.3, 0.02
Ta_in_reg, ωa_in_reg = 273.15 + 100.0, 0.005

function bcs_process(; Ta_in, ωa_in, Ta_init, ωa_init, W_init, ωd_init, Td_init)
    return [
        Ta(t, 0) ~ Ta_in, ωa(t, 0) ~ ωa_in,
        Dz(Td(t, 0)) ~ 0.0, Dz(W(t, 0)) ~ 0.0, # Dz(Td(t, L)) ~ 0.0,  Dz(W(t, L)) ~ 0.0, 
        Ta(0, z) ~ Ta_init, ωa(0, z) ~ ωa_init, W(0, z) ~ W_init, ωd(0, z) ~ ωd_init, Td(0, z) ~ Td_init
    ]
end

function bcs_regen(; Ta_in, ωa_in, Ta_init, ωa_init, W_init, ωd_init, Td_init)
    return [
        Ta(t, L) ~ Ta_in, ωa(t, L) ~ ωa_in,
        Dz(Td(t, 0)) ~ 0.0, Dz(W(t, 0)) ~ 0.0, # Dz(Td(t, 0)) ~ 0.0, Dz(W(t, 0)) ~ 0.0, 
        Ta(0, z) ~ Ta_init, ωa(0, z) ~ ωa_init, W(0, z) ~ W_init, ωd(0, z) ~ ωd_init, Td(0, z) ~ Td_init
    ]
end

# Simulation function
function run_cycles(total_cycles)
    results = []
    # Initial values
    Ta_init, ωa_init, W_init, ωd_init, Td_init = 273.15 + 60., 0.02, 0.1, 0.015, 273.15 + 60.

    dz = L/25

    vars = [ωa(t,z), Ta(t,z), W(t,z), Td(t,z), ωd(t,z), φw(t,z), Pws(t,z), Da(t,z), Ky(t,z), Qst(t,z)]

    for cycle in 1:total_cycles
        # Regeneration phase
        bcs_regen_phase = bcs_regen(Ta_in=Ta_in_reg, ωa_in=ωa_in_reg, Ta_init=Ta_init, ωa_init=ωa_init, W_init=W_init, ωd_init=ωd_init, Td_init=Td_init)
        @named pdesys_regen = PDESystem(eqs, bcs_regen_phase, [t ∈ IntervalDomain(0.0, 0.25*time_cycle), z ∈ IntervalDomain(0.0, L)], [t, z], vars)
        sol_regen = solve(discretize(pdesys_regen, MOLFiniteDifference([z=>dz], t)), Rodas5(); reltol=1e-6, abstol=1e-8, maxiters=10000)

        # Update for process phase
        Ta_init, ωa_init, W_init, ωd_init, Td_init = sol_regen[Ta(t,z)][end], sol_regen[ωa(t,z)][end], sol_regen[W(t,z)][end], sol_regen[ωd(t,z)][end], sol_regen[Td(t,z)][end]

        # Process phase
        bcs_process_phase = bcs_process(Ta_in=Ta_in_pro, ωa_in=ωa_in_pro, Ta_init=Ta_init, ωa_init=ωa_init, W_init=W_init, ωd_init=ωd_init, Td_init=Td_init)
        @named pdesys_process = PDESystem(eqs, bcs_process_phase, [t ∈ IntervalDomain(0.0, 0.75*time_cycle), z ∈ IntervalDomain(0.0, L)], [t, z], vars)
        sol_process = solve(discretize(pdesys_process, MOLFiniteDifference([z=>dz], t)), Rodas5(); reltol=1e-6, abstol=1e-8, maxiters=10000)

        # Update for next regeneration phase
        Ta_init, ωa_init, W_init, ωd_init, Td_init = sol_process[Ta(t,z)][end], sol_process[ωa(t,z)][end], sol_process[W(t,z)][end], sol_process[ωd(t,z)][end], sol_process[Td(t,z)][end]
        
        push!(results, (sol_regen, sol_process))
    end

    return results
end

# Run simulation for multiple cycles
cycles_number = 1
results = run_cycles(cycles_number)



#= Results
# reg
results[1][1][Ta(t, z)] .-273.15
results[1][1][Td(t, z)] .-273.15

results[1][1][ωa(t, z)] .*1000
results[1][1][ωd(t, z)] .*1000 

# pro
results[1][2][Ta(t, z)] .-273.15
results[1][2][Td(t, z)] .-273.15

results[1][2][ωa(t, z)] .*1000
results[1][2][ωd(t, z)] .*1000


# plot(results[1][1].t, results[1][1][Ta(t, z)][:, end] .-273.15, label="Air Outlet Temperature (°C)", xlabel="Time (s)", ylabel="Temperature (°C)", title="Outlet Air Temperature at z=L during Regeneration")
=#


using Plots


sol_reg = results[cycles_number][1]
sol_pro = results[cycles_number][2]

sol_reg[Ta(t, z)].-273.15
sol_pro[Ta(t, z)].-273.15

# outlet temperature
# reg_last second
t_values_reg = sol_reg.t
Ta_outlet_reg = sol_reg[Ta(t, z)][:, 1] .- 273.15  

# pro_last second
t_values_pro = sol_pro.t
Ta_outlet_pro = sol_pro[Ta(t, z)][:, end] .- 273.15 

plot(t_values_reg, Ta_outlet_reg, label="Reg Outlet T at z=L (°C)",
     xlabel="Time (s)", ylabel="Temperature (°C)", title="Outlet Air Temperature at z=L", color=:red)
hline!([Ta_in_reg - 273.15], label="Reg Inlet T (°C)", linestyle=:dash, color=:red)
plot!(t_values_pro .+ t_values_reg[end], Ta_outlet_pro, label="Pro Outlet T at z=L (°C)", color=:blue)
hline!([Ta_in_pro - 273.15], label="Pro Inlet T at z=0 (°C)", linestyle=:dash, legend=:right, color=:blue)


#= outlet humidity ratio
ωa_outlet_reg = sol_reg[ωa(t, z)][:, end] .*1000   
ωa_outlet_pro = sol_pro[ωa(t, z)][:, end] .*1000 
plot(t_values_reg, ωa_outlet_reg, label="Reg Inlet ω at z=L (g/kg)",
     xlabel="Time (s)", ylabel="humidity ratio (g/kg)", title="Outlet Air Humidity ratio at z=L", color=:red)
hline!([ωa_in_reg].*1000, label="Reg Inlet ω", linestyle=:dash, color=:red)
plot!(t_values_pro .+ t_values_reg[end], ωa_outlet_pro, label="Pro Outlet ω at z=L", color=:blue)
hline!([ωa_in_pro].*1000, label="Pro Inlet ω", linestyle=:dash, color=:blue)
=#

#= humidity ratio
ωa_outlet_regen = sol_regen[ωa(t, z)][:, end] .*1000   
ωa_outlet_process = sol_process[ωa(t, z)][:, end] .*1000 
plot(t_values_regen, ωa_outlet_regen, label="Reg Inlet ω at z=L (g/kg)",
     xlabel="Time (s)", ylabel="humidity ratio (g/kg)", title="Outlet Air Humidity ratio at z=L", color=:red)
hline!([ωa_in_reg].*1000, label="Reg Inlet ω", linestyle=:dash, color=:red)
plot!(t_values_process .+ t_values_regen[end], ωa_outlet_process, label="Pro Outlet ω at z=L", color=:blue)
hline!([ωa_in_pro].*1000, label="Pro Inlet ω", linestyle=:dash, color=:blue)

# pro阶段最后时刻 t=512 Ta随z的变化
z_values = range(0, stop=L, length=size(sol_process[Ta(t, z)], 2))
Ta_z_final = sol_process[Ta(t, z)][end, :] .- 273.15 
plot(z_values, Ta_z_final, label="Air Temperature at Final Time (°C)", xlabel="z (m)", ylabel="Temperature (°C)", title="Temperature Distribution along z at Final Time (Process Stage)")

# pro阶段最后时刻 t=512 ωa随z的变化
ωa_z_final = sol_process[ωa(t, z)][end, :] .*1000
plot(z_values, ωa_z_final, label="Air Humidity ratio at Final Time (°C)", xlabel="z (m)", ylabel="Humidity ratio (g/kg)", title="Humidity Ratio Distribution along z at Final Time (Process Stage)")
=#