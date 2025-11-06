#!/usr/bin/env julia
"""
Pellet Ablation Simulation
Author: wenjun xiong
Description:
    Demonstrate pellet ablation inside fusion plasma using PAM.jl.
    Outputs:
      - pellet_trajectory.png
      - pellet_radius.png
      - ablation_rate.png
"""

import PAM
import IMAS
using Plots

# ========== 设置绘图为无显示模式 ==========
ENV["GKSwstype"] = "100"  # 禁止 GUI 弹窗 (headless 环境)

# ========== 加载输入数据 ==========
example_dir = joinpath(@__DIR__, "examples")
json_file = joinpath(example_dir, "/gemini/code/PAM.jl/examples/template_D3D_1layer_2species.json")

if !isfile(json_file)
    error("❌ 示例文件未找到: $json_file")
end

println("🔹 Loading IMAS data from JSON...")
dd = IMAS.json2imas(json_file)
###通过RL训练的时候引入agent控制这几个参数, 接口到时候可以直接通过接口改动速度
dd.pellets.time_slice[].pellet[1].velocity_initial = 300.0
dd.pellets.time_slice[].pellet[1].layer[1].thickness = 0.5
dd.pellets.time_slice[].pellet[1].shape.size[1] = 0.5
dd.pellets.time_slice[].pellet[1].path_geometry.second_point.phi =0.0
dd.pellets.time_slice[].pellet[1].path_geometry.first_point.phi = 0.0
dd.pellets.time_slice[].pellet[1].path_geometry.first_point.z = 1.0
"""
表示为从装置中心向外的距离为r
dd.pellets.time_slice[].pellet[1].path_geometry.first_point.r = 45.0"
表示为从装置中心的垂直距离为多少
dd.pellets.time_slice[].pellet[1].path_geometry.first_point.z = 45.0"
"""
# ========== 模型参数 ==========
params = (
    t_start = 0.0,
    t_finish = 0.45,
    time_step = 0.001,
    drift_model = :HPI2,        # 可选 :Parks / :none
    Bt_dependance = true,
    update_plasma = false,
)

println("🚀 Running Pellet Ablation Model...")
pellet = PAM.run_PAM(dd; params...)

println("✅ 模拟完成！")

# ========== 绘图 ==========
# 1. 轨迹图 (r-z)
plot(pellet.r .+ pellet.R_drift, pellet.z,
     xlabel="r (m)", ylabel="z (m)", title="Pellet Trajectory (with drift)", legend=false)
savefig("pellet_trajectory_drift.png")
println("📈 Saved: pellet_trajectory.png")

# 2. 半径随时间变化
plot(pellet.time, pellet.radius,
    xlabel="Time (s)", ylabel="Radius (m)",
    title="Pellet Radius Evolution", lw=2)
savefig("pellet_radius.png")
println("📉 Saved: pellet_radius.png")

# 3. 烧蚀速率随时间变化
plot(pellet.time, pellet.ablation_rate,
    xlabel="Time (s)", ylabel="Ablation Rate (particles/s)",
    title="Pellet Ablation Rate", lw=2)
savefig("ablation_rate.png")
println("🔥 Saved: ablation_rate.png")

# 4. 密度沉积随时间
deposition = sum(pellet.density_source; dims=2)
plot(pellet.time, deposition,
    xlabel="Time (s)", ylabel="Integrated Source",
    title="Pellet Density Deposition", lw=2)
savefig("density_source.png")
println("🌫️ Saved: density_source.png")

println("\n🎯 All results saved in current directory:")
println("  - pellet_trajectory.png")
println("  - pellet_radius.png")
println("  - ablation_rate.png")
println("  - density_source.png")