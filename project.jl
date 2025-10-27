# plot_pam.jl
#
# 该脚本演示如何使用 RL-in-nuclear 仓库中的 PAM 模型生成示例图。
# 运行前请确保已在仓库目录下执行:
#   using Pkg
#   Pkg.activate(".")
#   Pkg.instantiate()

using PAM
using IMAS
using Plots

# -------------------------------------------------------------------------
# 配置部分
# -------------------------------------------------------------------------
# 模板文件路径，可使用 json 或 h5 模板
template_path = joinpath(@__DIR__, "examples", "template_D3D_1layer_2species.json")

# 仿真时间窗口（秒）
t_start  = 0.0
t_finish = 0.0045
time_step = 0.0001

# 漂移模型选择，可为 :none、:HPI2 或 :Parks
drift_model_choice = :HPI2

# 是否考虑磁场依赖及更新背景等离子体剖面
include_B_field_dependence = true
update_plasma_profiles     = true

# -------------------------------------------------------------------------
# 加载输入数据并运行 PAM 模型
# -------------------------------------------------------------------------
println("加载 IMAS 模板: ", template_path)
dd = IMAS.json2imas(template_path)

println("开始运行 PAM 模拟…")
pellet = PAM.run_PAM(
    dd;
    t_start  = t_start,
    t_finish = t_finish,
    time_step = time_step,
    drift_model = drift_model_choice,
    Bt_dependance = include_B_field_dependence,
    update_plasma = update_plasma_profiles
)

println("模拟完成，开始绘图…")

# -------------------------------------------------------------------------
# 图1: 颗粒半径随时间变化
# -------------------------------------------------------------------------
plt_radius = plot(
    pellet.time,
    pellet.radius,
    xlabel = "Time (s)",
    ylabel = "Pellet Radius (m)",
    title  = "Pellet Radius vs Time",
    legend = false
)
savefig(plt_radius, "pellet_radius_vs_time.png")
println("已保存图像: pellet_radius_vs_time.png")

# -------------------------------------------------------------------------
# 图2: 烧蚀率随时间变化
# -------------------------------------------------------------------------
plt_ablation = plot(
    pellet.time,
    pellet.ablation_rate,
    xlabel = "Time (s)",
    ylabel = "Ablation Rate (kg/s)",
    title  = "Ablation Rate vs Time",
    legend = false
)
savefig(plt_ablation, "ablation_rate_vs_time.png")
println("已保存图像: ablation_rate_vs_time.png")

# -------------------------------------------------------------------------
# 图3: (r,z) 平面沉积分布
# -------------------------------------------------------------------------
#
# -------------------------------------------------------------------------
# 图3: 径向烧蚀与沉积分布对比
# -------------------------------------------------------------------------

angles = [0, 45, 60]
models = [:PAM, :HPI2]

colors = Dict(
    :PAM  => [:red, :blue, :green],
    :HPI2 => [:red, :blue, :green]
)

linestyle = Dict(
    :ablation   => :dash,
    :deposition => :solid
)

# 新建空白图
plt_profile = plot(
    xlabel = "ρ",
    ylabel = "δn_D (m⁻³)",
    title  = "LFS Z=0.4 injection",
    legend = :topright
)

for model in models
    for (i, angle) in enumerate(angles)

        println("Running model=$(model), InjAngle=$(angle)° ...")
        pellet = PAM.run_PAM(
            dd;
            t_start  = t_start,
            t_finish = t_finish,
            time_step = time_step,
            drift_model = model == :PAM ? :Parks : :HPI2,
            inj_angle = angle,
            Bt_dependance = include_B_field_dependence,
            update_plasma = update_plasma_profiles
        )

        rho = pellet.rho
        ablation = pellet.ablation_profile
        deposition = pellet.deposition_profile

        # 颜色按角度区分，线型按 ablation/deposition 区分
        plot!(
            rho, ablation,
            color = colors[model][i],
            linestyle = linestyle[:ablation],
            label = "$(model), ablation, InjAngle=$(angle)°"
        )

        plot!(
            rho, deposition,
            color = colors[model][i],
            linestyle = linestyle[:deposition],
            label = "$(model), deposition, InjAngle=$(angle)°"
        )
    end
end

savefig(plt_profile, "pellet_density_profile.png")
println("✅ 已保存图像: pellet_density_profile.png")
println("所有图像生成完毕。")