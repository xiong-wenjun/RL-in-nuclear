module NuclearEnv

# 载入同目录下的 PAM.jl
include("PAM.jl")
using .PAM
import IMAS

export EnvState, init_env, reset_env, step_env, set_action_hook!

"""
EnvState 持有：
- dd:            IMAS.dd 原始数据（会被动作 hook 修改）
- pellet:        PAM.Pellet 上一轮 run_PAM 的结果（含全时间序列）
- k:             当前时间步索引（1-based）
- target_density: 目标等离子密度（用于奖励）
- t_start, t_finish, dt: 模拟时间区间与步长
- drift_model, Bt_dependance, update_plasma: 直接透传给 run_PAM 的参数
- action_hook!:  一个可选的函数 (dd, action, state) -> nothing，用于把动作写入 dd
                 若为 nothing，表示“无动作模式”，step_env 仅推进索引，不重算
"""
mutable struct EnvState
    dd::IMAS.dd
    pellet::PAM.Pellet
    k::Int
    target_density::Float64
    t_start::Float64
    t_finish::Float64
    dt::Float64
    drift_model::Symbol
    Bt_dependance::Bool
    update_plasma::Bool
    action_hook!::Union{Nothing, Function}
end

# -- 辅助：从 dd 与超参生成新的 pellet（完整时间序列）
function _rerun(dd::IMAS.dd; t_start::Float64, t_finish::Float64, dt::Float64,
                drift_model::Symbol, Bt_dependance::Bool, update_plasma::Bool)
    return PAM.run_PAM(dd;
        t_start=t_start, t_finish=t_finish, time_step=dt,
        drift_model=drift_model, Bt_dependance=Bt_dependance, update_plasma=update_plasma
    )
end

"""
init_env(dd; kwargs...)
  - 必需：dd :: IMAS.dd
  - 可选参数：
      t_start=0.0, t_finish=1e-3, dt=1e-5,
      drift_model=:none, Bt_dependance=true, update_plasma=false,
      target_density=1.5e19
"""
function init_env(dd::IMAS.dd; t_start::Float64=0.0, t_finish::Float64=1e-3, dt::Float64=1e-5,
                  drift_model::Symbol=:none, Bt_dependance::Bool=true, update_plasma::Bool=false,
                  target_density::Float64=1.5e19)
    pellet = _rerun(dd; t_start=t_start, t_finish=t_finish, dt=dt,
                    drift_model=drift_model, Bt_dependance=Bt_dependance, update_plasma=update_plasma)
    return EnvState(dd, pellet, 1, target_density, t_start, t_finish, dt,
                    drift_model, Bt_dependance, update_plasma, nothing)
end

"""
reset_env(state)
  - 重新根据当前 state 的超参/当前 dd 重跑 PAM，并把索引复位为 1
"""
function reset_env(state::EnvState)
    state.pellet = _rerun(state.dd; t_start=state.t_start, t_finish=state.t_finish, dt=state.dt,
                          drift_model=state.drift_model, Bt_dependance=state.Bt_dependance,
                          update_plasma=state.update_plasma)
    state.k = 1
    return state
end

"""
set_action_hook!(state, f!)
  - 注册动作写入函数：f!(dd::IMAS.dd, action, state::EnvState) -> nothing
  - 在 step_env 中若存在该 hook，则会先调用它把 action 写进 dd，
    随即重跑 run_PAM，再推进到下一步。
"""
function set_action_hook!(state::EnvState, f!::Function)
    state.action_hook! = f!
    return state
end

# -- 奖励：让 ne[k] 靠近 target_density（你可根据课题改写）
@inline function _reward(state::EnvState, k::Int)
    ne = state.pellet.ne[k]      # 粒子数密度 (来自 Pellet 轨迹)
    return -abs(ne - state.target_density)
end

# -- 终止条件：丸完全烧尽 或 时间走到末尾
@inline function _done(state::EnvState, k::Int)
    rad = state.pellet.radius[k]
    if rad <= 0.0
        return true
    end
    return k >= length(state.pellet.time)
end

"""
step_env(state, action)
  - 若已设置 action_hook!：先调用 hook 把动作写入 dd，然后重跑 PAM 得到新 pellet
  - 不论是否重跑：把索引推进一格，返回 (new_state, reward, done)
  - action 可以是任意可迭代对象（来自 Python 的 list / numpy array 都可以）
"""
function step_env(state::EnvState, action)
    # 1) 处理动作（可选）
    if state.action_hook! !== nothing
        # 宽松接收 Python/JL 的各种容器，转成 Vector{Float64}
        act = Float64.(collect(action))
        state.action_hook!((state.dd)::IMAS.dd, act, state)
        # 重跑以反映动作影响
        state.pellet = _rerun(state.dd; t_start=state.t_start, t_finish=state.t_finish, dt=state.dt,
                              drift_model=state.drift_model, Bt_dependance=state.Bt_dependance,
                              update_plasma=state.update_plasma)
        # 可选：把 k 对齐到动作发生后的下一帧，这里简单设为 k=k（即使用当前索引）
    end

    # 2) 推进一步（注意 Pellet 内部是 1-based 索引）
    k_next = min(state.k + 1, length(state.pellet.time))
    state.k = k_next

    # 3) 计算奖励/终止
    r = _reward(state, state.k)
    done = _done(state, state.k)

    return state, r, done
end

end # module