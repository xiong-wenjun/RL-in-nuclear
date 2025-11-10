module NuclearEnv
using PAM
using JSON3  # 如果要通过 JSON 传递状态给 Python（推荐）

export EnvState, init_env, reset_env, step_env

mutable struct EnvState
    pellet_pos::Float64
    pellet_vel::Float64
    pellet_radius::Float64
    plasma_density::Float64
    t::Float64
end

function init_env()
    return EnvState(0.0, 0.0, 1.0, 1e19, 0.0)
end

function reset_env()
    return init_env()
end

function step_env(state::EnvState, action::Vector{Float64})
    v, angle = action
    dt = 1e-4

    # 调用 PAM 模型计算
    pellet_state = PAM.step_pellet(state.pellet_pos, v, angle, dt)
    ablation_rate = PAM.compute_ablation_rate(pellet_state)
    new_density = state.plasma_density + ablation_rate * dt

    reward = -abs(new_density - 1.5e19)
    done = (pellet_state.radius < 0.01) || (state.t > 0.1)

    new_state = EnvState(
        pellet_state.pos,
        pellet_state.vel,
        pellet_state.radius,
        new_density,
        state.t + dt
    )

    return new_state, reward, done
end

end # module
