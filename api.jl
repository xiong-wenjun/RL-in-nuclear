# ================= api.jl =================
module PelletAPI

using JSON3
include("PAM.jl")   # 载入你的主模块
include("run.jl")

"Run a single pellet simulation and return results as Dict"
function run_pellet_sim(v::Float64, r::Float64, angle::Float64, Z::Float64)
    # 调用现有仿真逻辑，这里用 run_simulation 或 project/task 等函数
    result = run_simulation(v, r, angle, Z)   # 假设你仓库已有此函数
    # 返回关键物理量
    return Dict(
        "rho"        => result.rho,
        "deposition" => result.deposition,
        "ablation"   => result.ablation
    )
end

end # module
