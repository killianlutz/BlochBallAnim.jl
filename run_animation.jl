using Pkg
Pkg.activate("./venv_BlochAnimation")
# Pkg.instantiate() # first use: resolves appropriate package versions

include("./src/BlochAnimation.jl")
using .BlochAnimation
import GLMakie.theme_dark
import GLMakie.Theme

figure_theme = theme_dark();
# figure_theme = Theme();

############ CHOOSE specific GKSL ode parameters
### 2D control ω and control hamiltonians: Pauli X/2 and Y/2
### free Hamiltonian (must be hermitian)
H0 = [-1 2; 2 1] 
### jump operators (any number of 2x2 matrices, damping rates γ_k = 1)
# J_minus and J_plus with asymmetric damping rates
h = [
    [0 3; 0 0], [0 0; 1 0]
]
# # J_minus and J_plus
# h = [
#     [0 1; 0 0], [0 0; 1 0]
# ]
## J_z
# h = [
#     [0.5 0; 0 -0.5]
# ]

############ Run animation
BlochAnimation.bloch_animation(; figure_theme, H0, h)