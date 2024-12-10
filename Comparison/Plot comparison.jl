cd(@__DIR__) 

using CSV, DataFrames, Plots, Interpolations

df_MTK = CSV.read("MOL_results.csv", DataFrame)
df_lit = CSV.read("Literature_data.csv", DataFrame)

# 分段数据 (t <= 128 为 reg，再生阶段；t > 128 为 pro，处理阶段)
df_MTK_reg = filter(row -> row.time < 128, df_MTK)
df_MTK_pro = filter(row -> row.time >= 128, df_MTK)

df_lit_reg1 = filter(row -> !ismissing(row.Time1) && row.Time1 < 128, df_lit)
df_lit_pro1 = filter(row -> !ismissing(row.Time1) && row.Time1 >= 128, df_lit)

df_lit_reg2 = filter(row -> !ismissing(row.Time2) && row.Time2 < 128, df_lit)
df_lit_pro2 = filter(row -> !ismissing(row.Time2) && row.Time2 >= 128, df_lit)

# 绘制原始数据曲线
# T
plot(df_MTK_reg.time, df_MTK_reg.temperature, label="Julia: T reg", xlabel="Time (s)", ylabel="Temperature (°C)", color=:red)
plot!(df_MTK_pro.time, df_MTK_pro.temperature, label="Julia: T pro", color=:red)

plot!(df_lit_reg1.Time1, df_lit_reg1.temperature, label="Literature: T reg", color=:gray)
plot!(df_lit_pro1.Time1, df_lit_pro1.temperature, label="Literature: T pro", color=:gray)

# ω
plot(df_MTK_reg.time, df_MTK_reg.humidity_ratio, label="Julia: ω reg", xlabel="Time (s)", ylabel="Humidity ratio (g/kg)", color=:red)
plot!(df_MTK_pro.time, df_MTK_pro.humidity_ratio, label="Julia: ω pro", color=:red)

plot!(df_lit_reg2.Time2, df_lit_reg2.humidity_ratio .* 1000, label="Literature: ω reg", color=:gray)
plot!(df_lit_pro2.Time2, df_lit_pro2.humidity_ratio .* 1000, label="Literature: ω pro", color=:gray)


# Error of T
# 创建插值函数
interp_lit_reg1 = LinearInterpolation(df_lit_reg1.Time1, df_lit_reg1.temperature, extrapolation_bc=NaN)
interp_lit_pro1 = LinearInterpolation(df_lit_pro1.Time1, df_lit_pro1.temperature, extrapolation_bc=NaN)

# 对齐 Literature 数据到 MTK 数据的时间
lit_temp_reg1 = [interp_lit_reg1(t) for t in df_MTK_reg.time]
lit_temp_pro1 = [interp_lit_pro1(t) for t in df_MTK_pro.time]

# 计算误差
error_reg_T = abs.((df_MTK_reg.temperature .- lit_temp_reg1)./ lit_temp_reg1)  
error_pro_T = abs.((df_MTK_pro.temperature .- lit_temp_pro1) ./ lit_temp_pro1)

# 绘制再生阶段 (reg) 和处理阶段 (pro) 的误差曲线
plot(df_MTK_reg.time, error_reg_T .* 100, label="Error (Regeneration)", xlabel="Time (s)", ylabel="Error of temperature (%)", color=:black, linestyle=:dash)
plot!(df_MTK_pro.time, error_pro_T .* 100, label="Error (Process)", color=:black, linestyle=:dash)


# Error of ω
interp_lit_reg2 = LinearInterpolation(df_lit_reg2.Time2, df_lit_reg2.humidity_ratio, extrapolation_bc=NaN)
interp_lit_pro2 = LinearInterpolation(df_lit_pro2.Time2, df_lit_pro2.humidity_ratio, extrapolation_bc=NaN)

# 对齐 Literature 数据到 MTK 数据的时间
lit_temp_reg2 = [interp_lit_reg2(t) for t in df_MTK_reg.time]
lit_temp_pro2 = [interp_lit_pro2(t) for t in df_MTK_pro.time]

# 计算误差
error_reg_ω = abs.((df_MTK_reg.humidity_ratio ./ 1000 .- lit_temp_reg2)./ lit_temp_reg2)  
error_pro_ω = abs.((df_MTK_pro.humidity_ratio ./ 1000 .- lit_temp_pro2) ./ lit_temp_pro2)

# 绘制再生阶段 (reg) 和处理阶段 (pro) 的误差曲线
plot(df_MTK_reg.time, error_reg_ω .* 100, label="Error (Regeneration)", xlabel="Time (s)", ylabel="Error of humidity ratio (%)", color=:black, linestyle=:dash)
plot!(df_MTK_pro.time, error_pro_ω .* 100, label="Error (Process)", color=:black, linestyle=:dash)






#=
# Use thershold to create NaN for discontinuity
threshold = 25
for i in 1:(length(df_MTK.temperature)-1)
    if abs(df_MTK.temperature[i+1] - df_MTK.temperature[i]) > threshold
        x = insert!(df_MTK.time, i+1, NaN)
        y = insert!(df_MTK.temperature, i+1, NaN)
    end

    if abs(df_lit.T[i+1] - df_lit.T[i]) > threshold
        x = insert!(df_lit.Time1, i+1, NaN)
        y = insert!(df_lit.T, i+1, NaN)
    end
end


plot(df_MTK.time, df_MTK.temperature, label="Julia T",
     xlabel="Time (s)", ylabel="Temperature (°C)", color=:red)
plot!(df_lit.Time1, df_lit.T, label="Literature T",
     xlabel="Time (s)", ylabel="Temperature (°C)", color=:gray)
=#