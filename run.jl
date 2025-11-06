using Plots
using JSON
using IMAS
using PAM   # 确保已经在 LOAD_PATH 或 dev 安装好 PAM 模块

# 1️⃣ 读取 IMAS 数据（包含 equilibrium, core_profiles, pellets）
dd = IMAS.hdf2imas("examples/template_D3D_1layer_2species.h5")
##dd.pellets.time_slice[].pellet[1].velocity_initial = 200.0 
# 2️⃣ 设置模拟参数
params = (
    t_start = 0.0,
    t_finish = 0.0045,
    time_step = 0.0001,
    drift_model = :HPI2,        # 可选 :Parks / :none
    Bt_dependance = true,
    update_plasma = false,
)

# 3️⃣ 运行 PAM 主函数
pellet = PAM.run_PAM(dd; params...)

# 4️⃣ 绘制托卡马克边界
eq = dd.equilibrium.time_slice[1]
r_bnd = eq.boundary.outline.r
z_bnd = eq.boundary.outline.z
pam.run(ablation_model='PAM', deposition_model='Gaussian2D')
pam.plot_profiles(compare='HPI2')
# 5️⃣ 绘制轨迹
plot(r_bnd, z_bnd,
     color=:blue, lw=2, label="Plasma boundary", aspect_ratio=:equal)

plot!(pellet.r, pellet.z,
      color=:red, lw=2, label="Pellet trajectory")

plot!(pellet.r .+ pellet.R_drift, pellet.z,
      color=:orange, lw=2, ls=:dot, label="Drift corrected trajectory")

# 起点标记
scatter!([pellet.r[1]], [pellet.z[1]],
         color=:green, markersize=8, label="Injection point")

xlabel!("R (m)")
ylabel!("Z (m)")
title!("Pellet Injection Trajectory Simulation")
plot!(grid = true)

# 6️⃣ 保存图像
savefig("pellet_trajectory.png")

println("✅ 图像已保存为 pellet_trajectory.png")