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

# domain
L = 0.2             # Thickness of desiccant wheel (m)
time_cycle = 512.0
time_pro_per = 0.75
time_reg_per = 1 - time_pro_per

# 定义不同的时间域
domain_regeneration = [t ∈ IntervalDomain(0.0, time_cycle * time_reg_per), z ∈ IntervalDomain(0.0, L)]
domain_process = [t ∈ IntervalDomain(0.0, time_cycle * time_pro_per), z ∈ IntervalDomain(0.0, L)]

Ta_pro_in = 273.15 + 34.3
ωa_pro_in = 0.02
Ta_reg_in = 273.15 + 100.0  # 再生空气入口温度（假设）
ωa_reg_in = 0.005          # 再生空气入口含湿量（假设）

# 定义再生阶段的边界条件模板
function bcs_regeneration(; Ta_in_bc, ωa_in_bc, Ta_ic, ωa_ic, W_ic, ωd_ic, Td_ic)
       return [
           Ta(t, 0) ~ Ta_in_bc,
           ωa(t, 0) ~ ωa_in_bc,
           Dz(Td(t, 0)) ~ 0.0,
           Dz(Td(t, L)) ~ 0.0,
           Dz(W(t, 0)) ~ 0.0,
           Dz(W(t, L)) ~ 0.0,
           ωa(0, z) ~ ωa_ic,
           Ta(0, z) ~ Ta_ic,
           W(0, z) ~ W_ic,
           ωd(0, z) ~ ωd_ic,
           Td(0, z) ~ Td_ic
       ]
end
   
   # 定义处理阶段的边界条件模板
function bcs_process(; Ta_in_bc, ωa_in_bc, Ta_ic, ωa_ic, W_ic, ωd_ic, Td_ic)
       return [
           Ta(t, 0) ~ Ta_in_bc,
           ωa(t, 0) ~ ωa_in_bc,
           Dz(Td(t, 0)) ~ 0.0,
           Dz(Td(t, L)) ~ 0.0,
           Dz(W(t, 0)) ~ 0.0,
           Dz(W(t, L)) ~ 0.0,
           ωa(0, z) ~ ωa_ic,
           Ta(0, z) ~ Ta_ic,
           W(0, z) ~ W_ic,
           ωd(0, z) ~ ωd_ic,
           Td(0, z) ~ Td_ic
       ]
end
   
# 主循环函数
function run_cycles(total_cycles)
       results = []  # 存储每个周期的结果
       # t=0的初始条件
       Ta_reg_ic, ωa_reg_ic, W_reg_ic, ωd_reg_ic, Td_reg_ic = Ta_reg_in, ωa_reg_in, 0.01, 0.02, Ta_reg_in  # 初始条件
       
       for cycle in 1:total_cycles
           # 再生阶段
           bcs_regen = bcs_regeneration(
              Ta_in_bc = Ta_reg_in, 
              ωa_in_bc = ωa_reg_in, 
              Ta_ic = Ta_reg_ic, 
              ωa_ic = ωa_reg_ic, 
              W_ic = W_reg_ic,  
              ωd_ic = ωd_reg_ic, 
              Td_ic = Td_reg_ic)
           @named DW_regen_pdesys = PDESystem(eqs, bcs_regen, domain_regeneration, [t, z], vars)
           odeprob_regen = discretize(DW_regen_pdesys, discretization)
           sol_regen = solve(odeprob_regen, Rodas5(); reltol=1e-6, abstol=1e-8, maxiters=10000)
   
           # 更新初始条件为再生阶段的终点
           Ta_pro_ic = sol_regen[Ta(t,z)][end]
           ωa_pro_ic = sol_regen[ωa(t,z)][end]
           W_pro_ic = sol_regen[W(t,z)][end]
           ωd_pro_ic = sol_regen[ωd(t,z)][end]
           Td_pro_ic = sol_regen[Td(t,z)][end]
   
           # 处理阶段
           bcs_pro = bcs_process(
              Ta_in_bc = Ta_pro_in, 
              ωa_in_bc = ωa_pro_in, 
              Ta_ic = Ta_pro_ic, 
              ωa_ic = ωa_pro_ic, 
              W_ic = W_pro_ic, 
              ωd_ic = ωd_pro_ic,
              Td_ic = Td_pro_ic)
           @named DW_process_pdesys = PDESystem(eqs, bcs_pro, domain_process, [t, z], vars)
           odeprob_process = discretize(DW_process_pdesys, discretization)
           sol_process = solve(odeprob_process, Rodas5(); reltol=1e-6, abstol=1e-8, maxiters=10000)
   
           # 更新初始条件为处理阶段的终点
           Ta_reg_ic = sol_process[Ta(t,z)][end]
           ωa_reg_ic = sol_process[ωa(t,z)][end]
           W_reg_ic = sol_process[W(t,z)][end]
           ωd_reg_ic = sol_process[ωd(t,z)][end]
           Td_reg_ic = sol_process[Td(t,z)][end]
   
           # 存储当前周期的结果
           push!(results, (sol_regen, sol_process))
       end
   
       return results  # 返回每个周期的结果
end

total_cycles = 1     # 总循环次数，定义运行多少个“再生-处理”周期
# 运行多个周期
results = run_cycles(total_cycles)


using Plots

## Results
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

# 第一个周期的再生和处理解
sol_regen = results[1][1]
sol_process = results[1][2]

# 再生阶段的时间和出口温度（最后一列）
t_values_regen = sol_regen.t
Ta_outlet_regen = sol_regen[Ta(t, z)][:, end] .- 273.15  

# 处理阶段的时间和出口温度（最后一列）
t_values_process = sol_process.t
Ta_outlet_process = sol_process[Ta(t, z)][:, end] .- 273.15 

# reg出口温度
plot(t_values_regen, Ta_outlet_regen, label="Reg Inlet T at z=L (°C)",
     xlabel="Time (s)", ylabel="Temperature (°C)", title="Outlet Air Temperature at z=L", color=:red)
# reg进口温度
hline!([Ta_reg_in - 273.15], label="Reg Inlet T (°C)", linestyle=:dash, color=:red)
# pro出口温度
plot!(t_values_process .+ t_values_regen[end], Ta_outlet_process, label="Pro Outlet T at z=L (°C)", color=:blue)
# pro进口温度
hline!([Ta_pro_in - 273.15], label="Pro Inlet T at z=0 (°C)", linestyle=:dash, legend=:bottomright, color=:blue)

# humidity ratio
ωa_outlet_regen = sol_regen[ωa(t, z)][:, end] .*1000   
ωa_outlet_process = sol_process[ωa(t, z)][:, end] .*1000 
plot(t_values_regen, ωa_outlet_regen, label="Reg Inlet ω at z=L (g/kg)",
     xlabel="Time (s)", ylabel="humidity ratio (g/kg)", title="Outlet Air Humidity ratio at z=L", color=:red)
hline!([ωa_reg_in].*1000, label="Reg Inlet ω", linestyle=:dash, color=:red)
plot!(t_values_process .+ t_values_regen[end], ωa_outlet_process, label="Pro Outlet ω at z=L", color=:blue)
hline!([ωa_pro_in].*1000, label="Pro Inlet ω", linestyle=:dash, color=:blue)

# pro阶段最后时刻 t=512 Ta随z的变化
z_values = range(0, stop=L, length=size(sol_process[Ta(t, z)], 2))
Ta_z_final = sol_process[Ta(t, z)][end, :] .- 273.15 
plot(z_values, Ta_z_final, label="Air Temperature at Final Time (°C)", xlabel="z (m)", ylabel="Temperature (°C)", title="Temperature Distribution along z at Final Time (Process Stage)")

# pro阶段最后时刻 t=512 ωa随z的变化
ωa_z_final = sol_process[ωa(t, z)][end, :] .*1000
plot(z_values, ωa_z_final, label="Air Humidity ratio at Final Time (°C)", xlabel="z (m)", ylabel="Humidity ratio (g/kg)", title="Humidity Ratio Distribution along z at Final Time (Process Stage)")