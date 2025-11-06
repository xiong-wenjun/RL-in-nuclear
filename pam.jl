using PAM 
using Plots

import IMAS

dd = IMAS.json2imas("/gemini/code/PAM.jl/examples/ods.txt")
dd.core_profiles.profiles_1d[1].time = dd.equilibrium.time_slice[1].time;
# save dd pellet structure for future use
old_pellets = deepcopy(dd.pellets);
