from juliacall import Main as jl

# 1. 载入环境模块
jl.include("/gemini/code/PAM.jl/NuclearEnv.jl")
env = jl.NuclearEnv

# 2. 假设你已有 dd 对象（示例: 从 IMAS 读取或 mock）
# 例如：
# jl.seval("using IMAS; dd = IMAS.dd()")  # 如果要创建空的 dd
# dd = jl.seval("dd")                     # 获取 Julia 变量到 Python
# 实际中应替换为你加载的真实 dd

# 这里暂时放一个占位对象以避免报错（你可以替换成真实 dd）

jl.seval("""
using IMAS, JSON3

# 假设文件路径是 /gemini/code/PAM.jl/imas_case.json
dd = IMAS.json2imas("/gemini/code/PAM.jl/examples/template_D3D_1layer_2species.json")
# 确认一下
println("✅ 成功载入 IMAS JSON 文件")
println("equilibrium.time_slice 数量: ", length(dd.equilibrium.time_slice))
println("pellets.time_slice 数量: ", length(dd.pellets.time_slice))
""")
dd = jl.seval("dd")

# 3. 初始化环境
state = jl.NuclearEnv.init_env(dd,
                     t_start=0.0,
                     t_finish=1e-3,
                     dt=1e-5,
                     drift_model=jl.Symbol("none"),
                     Bt_dependance=True,
                     update_plasma=False,
                     target_density=1.5e19)

# 4. 迭代 step（不带动作 hook）
for i in range(10):
    # 直接传 Python list 即可，Julia 会自动 Float64.(collect(action))
    action = [200.0, 0.1]
    new_state, reward, done = env.step_env(state, action)
    print(f"Step {i}: reward={reward}, done={done}")
    state = new_state
    if done:
        break