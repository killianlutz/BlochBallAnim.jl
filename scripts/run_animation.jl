using Pkg
Pkg.activate("./venv_BlochAnimation")
# Pkg.instantiate() # first use: resolves appropriate package versions

include("../src/BlochAnimation.jl")
using .BlochAnimation
import GLMakie.theme_dark
import GLMakie.Theme
import GLMakie.with_theme

figure_theme = theme_dark();
# figure_theme = Theme();

############ CHOOSE specific GKSL ODE parameters
free_hamiltonian = [-1 2; 2 1] # hermitian
collapse_operators = [
    [0 sqrt(3); 0 0], [0 0; 1 0]
] # any number of 2x2 matrices, unit damping rates

############ Run animation
with_theme(figure_theme) do
    BlochAnimation.bloch_animation(; H0=free_hamiltonian, h=collapse_operators)
end