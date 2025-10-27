# =========================================================================
#  plot_pam_final_fixed.jl
#  兼容当前 PAM.jl 源码的版本：基于 density_source 和 ρ 计算 ablation/deposition 剖面
# =========================================================================

using PAM
using IMAS
using Plots
using Statistics

# -------------------------------------------------------------------------
# 基本配置
# -------------------------------------------------------------------------
template_path = joinpath(@__DIR__, "examples", "template_D3D_1layer_2species.json")

t_start  = 0.0
t_finish = 0.0045
time_step = 0.0001
include_B_field_dependence = true
update_plasma_profiles     = true

# -------------------------------------------------------------------------
# 载入 IMAS 模板
# -------------------------------------------------------------------------
println("加载 IMAS 模板: ", template_path)
dd = IMAS.json2imas(template_path)

# -------------------------------------------------------------------------
# 辅助函数：积分沉积分布、分箱烧蚀分布
# -------------------------------------------------------------------------

# (A) 从 pellet.density_source 计算沉积剖面 (沿时间积分)
function deposition_profile_vs_rho(dd, pellet)
    depot_ts = pellet.density_source
    Nt, Nsurf = size(depot_ts)

    t = pellet.time
    dt = [t[1]; diff(t)]
    depot_int = vec(sum(depot_ts .* dt, dims=1))   # 时间积分

    # 获取每个磁通面的代表性 ρ
    eqt = dd.equilibrium.time_slice[]
    fw  = IMAS.first_wall(dd.wall)
    surfaces = IMAS.trace_surfaces(eqt, fw...)
    _, _, RHO = IMAS.ρ_interpolant(eqt)

    ρ_surf = similar(depot_int)
    for (i, s) in enumerate(surfaces)
        r̄, z̄ = mean(s.r), mean(s.z)
        ρ_surf[i] = RHO(r̄, z̄)
    end

    ord = sortperm(ρ_surf)
    return ρ_surf[ord], depot_int[ord]
end

# (B) 从 ablation_rate 和 ρ(t) 得到烧蚀分布
function ablation_profile_vs_rho(pellet; nbins::Int=200)
    ρt = pellet.ρ
    t  = pellet.time
    dt = [t[1]; diff(t)]
    G  = pellet.ablation_rate .* dt
    mask = (ρt .>= 0.0) .& (ρt .<= 1.0)
    ρt, G = ρt[mask], G[mask]

    edges = range(0.0, 1.0; length=nbins+1) |> collect
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    accum = zeros(nbins)
    for (ρ, g) in zip(ρt, G)
        idx = clamp(searchsortedlast(edges, ρ), 1, nbins)
        accum[idx] += g
    end
    return centers, accum
end

# -------------------------------------------------------------------------
# 绘图主程序
# -------------------------------------------------------------------------
models = [:PAM, :HPI2]
angles = [0, 45, 60]
colors = Dict(:PAM => [:red, :blue, :green], :HPI2 => [:red, :blue, :green])

plt = plot(
    xlabel = "ρ",
    ylabel = "normalized δnᴰ (arb.)",
    title  = "LFS Z=0.4 injection",
    legend = :topright,
    lw = 2,
    grid = false
)

for model in models
    for (i, angle) in enumerate(angles)
        local pellet
        println("🚀 模型=$(model), 注入角度=$(angle)° ...")

        # 注意：PAM.run_PAM 不接受 inj_angle 参数
        # 如果要扫角度，请修改模板文件中的几何参数 path_geometry
        pellet = PAM.run_PAM(
            dd;
            t_start  = t_start,
            t_finish = t_finish,
            time_step = time_step,
            drift_model = model == :PAM ? :Parks : :HPI2,
            Bt_dependance = include_B_field_dependence,
            update_plasma = update_plasma_profiles
        )

        # --- 计算剖面 ---
        ρ_depo, depo = deposition_profile_vs_rho(dd, pellet)
        ρ_abl , abl  = ablation_profile_vs_rho(pellet)

        scale = maximum(depo) > 0 ? maximum(depo) : 1.0
        depo_n = depo ./ scale
        abl_n  = (maximum(abl) > 0 ? abl ./ maximum(abl) : abl)

        # --- 绘图 ---
        plot!(
            ρ_depo, depo_n,
            color = colors[model][i],
            linestyle = :solid,
            label = "$(model), deposition, $(angle)°"
        )
        plot!(
            ρ_abl, abl_n,
            color = colors[model][i],
            linestyle = :dash,
            label = "$(model), ablation, $(angle)°"
        )
    end
end

# 保存图像
savefig(plt, "pellet_density_profile.png")
println("✅ 已生成图像: pellet_density_profile.png")