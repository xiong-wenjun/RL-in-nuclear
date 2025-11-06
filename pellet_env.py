# ================= pellet_env.py =================
import gymnasium as gym
from gymnasium import spaces
import numpy as np
import yaml
from julia.api import Julia

# 初始化 Julia 环境（推荐方式）
jl = Julia(compiled_modules=False)
from julia import Main
Main.include("api.jl")  # 直接加载你刚写的 api.jl 模块

class PelletInjectionEnv(gym.Env):
    """Gym wrapper for RL-in-nuclear pellet injection simulation."""

    metadata = {"render_modes": ["human"]}

    def __init__(self, config_path=None):
        super().__init__()
        # -------------------------
        # 读取环境背景参数
        # -------------------------
        if config_path:
            self.config = yaml.safe_load(open(config_path))
        else:
            self.config = {"Z": 0.4, "plasma_Te": 1000, "plasma_ne": 1e19}

        # 动作空间: [velocity, radius, angle]
        self.action_space = spaces.Box(
            low=np.array([100.0, 0.05, 0.0]),
            high=np.array([300.0, 0.15, 60.0]),
            dtype=np.float64
        )

        # 观测空间: ρ分布 (100个点)
        self.observation_space = spaces.Box(
            low=0.0, high=np.inf, shape=(100,), dtype=np.float64
        )

        self.state = np.zeros(100)
        self.reward_range = (0, np.inf)

    # 重置环境
    def reset(self, *, seed=None, options=None):
        super().reset(seed=seed)
        self.state = np.zeros(100)
        return self.state, {}

    # 执行动作
    def step(self, action):
        v, r, ang = map(float, action)
        Z = float(self.config.get("Z", 0.4))

        # 调用 Julia 仿真
        result = Main.API.run_pellet_sim(v, r, ang, Z)
        deposition = np.array(result["deposition"], dtype=np.float64)

        # 计算奖励（沉积积分）
        reward = float(np.sum(deposition))

        # 每次仿真为独立episode
        done = True
        info = {"velocity": v, "radius": r, "angle": ang}

        self.state = deposition
        return self.state, reward, done, False, info

    # 可选渲染函数
    def render(self):
        import matplotlib.pyplot as plt
        plt.plot(np.linspace(0, 1, len(self.state)), self.state)
        plt.xlabel("ρ")
        plt.ylabel("Δnᴰ (m⁻³)")
        plt.show()
