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
# 加载模板
# -------------------------------------------------------------------------
println("加载 IMAS 模板: ", template_path)
dd = IMAS.json2imas(template_path)

# -------------------------------------------------------------------------
# 工具函数：计算 ablation 与 deposition 剖面
# -------------------------------------------------------------------------
function deposition_profile_vs_rho(dd, pellet)
    depot_ts = pellet.density_source
    t = pellet.time
    dt = [t[1]; diff(t)]
    depot_int = vec(sum(depot_ts .* dt, dims=1))

    eqt = dd.equilibrium.time_slice[]
    fw  = IMAS.first_wall(dd.wall)
    surfaces = IMAS.trace_surfaces(eqt, fw...)
    _, _, RHO = IMAS.ρ_interpolant(eqt)

    ρ_surf = similar(depot_int)
    for (i, s) in enumerate(surfaces)
        ρ_surf[i] = mean([RHO(r, z ) for (r, z) in zip(s.r, s.z)])
    end

    ord = sortperm(ρ_surf)
    return ρ_surf[ord], depot_int[ord]
end

function ablation_profile_vs_rho(pellet; nbins=200)
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
markers = Dict(:PAM => :utriangle, :HPI2 => :dtriangle)

plt = plot(
    xlabel = "ρ",
    ylabel = "δnᴰ (m⁻³)",
    title  = "LFS Z=0.4 injection",
    legend = :topright,
    lw = 2,
    grid = false,
    framestyle = :box,
    tickfont = font(10, "Arial"),
    guidefont = font(12, "Arial"),
    titlefont = font(14, "Arial"),
)

for model in models
    for (i, angle) in enumerate(angles)
        println("🚀 模型=$(model), 角度=$(angle)° ...")
        local pellet = PAM.run_PAM(
            dd;
            t_start  = t_start,
            t_finish = t_finish,
            time_step = time_step,
            drift_model = model == :PAM ? :Parks : :HPI2,
            Bt_dependance = include_B_field_dependence,
            update_plasma = update_plasma_profiles
        )

        ρ_depo, depo = deposition_profile_vs_rho(dd, pellet)
        ρ_abl , abl  = ablation_profile_vs_rho(pellet)

        scale = maximum(depo) > 0 ? maximum(depo) : 1.0
        depo_n = depo ./ scale
        abl_n  = (maximum(abl) > 0 ? abl ./ maximum(abl) : abl)

        # 绘制虚线 + 三角标记（ablation）
        plot!(
            ρ_abl, abl_n,
            color = colors[model][i],
            linestyle = :dash,
            marker = markers[model],
            markersize = 4,
            markerstrokecolor = colors[model][i],
            label = "$(model), ablation, InjAngle=$(angle)°"
        )

        # 绘制实线 + 反三角标记（deposition）
        plot!(
            ρ_depo, depo_n,
            color = colors[model][i],
            linestyle = :solid,
            marker = markers[model],
            markersize = 4,
            fillalpha = 0.4,
            label = "$(model), deposition, InjAngle=$(angle)°"
        )
    end
end

savefig(plt, "pellet_density_profile_pretty.png")
println("✅ 已生成图像: pellet_density_profile_pretty.png")