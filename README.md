# BlochBallAnim.jl

## Goal: open qubit quantum control made interactive and intuitive!


https://github.com/killianlutz/BlochBallAnim.jl/assets/152091888/b39c8a53-96e6-474f-8b27-4a522892be73


By adopting the [Bloch ball](https://en.wikipedia.org/wiki/Bloch_sphere) perspective, this project showcases the dynamics of open markovian quantum systems $\varrho(t)$ undergoing classical coherent control $\Omega(t)$.

Starting from an initial density operator $\varrho_0 \in \mathscr{M}(\mathbb{C})$, the time-evolution of such a system is modelled via the GKS-Lindblad ordinary differential equation

$$\dot{\varrho} = \mathcal{L}(\Omega(t),\varrho), \quad 0 \leq t \leq T$$

where the controlled Hamiltonian is 
$H_{\Omega} = H_0 + \frac{1}{2}(\Omega_x X + \Omega_y Y)$ 
and the controlled Lindbladian in diagonal form writes

$$\mathcal{L}(\Omega, \varrho) = -i\left(H_{\Omega}\varrho-\varrho H_{\Omega}\right) + \sum_{k} \gamma_k \left(2h_k \varrho h_k^* - \left(h_k^*h_k \varrho + \varrho h_k^*h_k\right) \right).$$

Using the Pauli matrices $X, Y, Z$, the time-volution of qubit density operators $\varrho(t)$ may be embedded as a curve $x(t)$ in the closed unit ball of $\mathbb{R}^3$, also called Bloch ball.

Once embedded into the Bloch ball, the GKS-lindblad ODE is a affine and of the form 
$$\dot{x} = Ax + B(\Omega, x) + b.$$
The `Velocity` button displays the vector-field $x \mapsto Ax + B(\Omega, x) + b$ in $\mathbb{R}^3$.

## This animation provides answers to
* Interpretation of coherent controls $\Omega(t)$?
* Effect of a particular choice of collapse operators $h_k$?
* To what extent are Lindbladians $\mathcal{L}$ `dissipative' (unitality)?
* Action of most important quantum qubit gates $\mathcal{Q}\in SU(4)$?

## For whom were those animations designed?
Physicists or mathematicians, students learning about quantum mechanics as well as researchers studying qubits and seeking to gain intuition.

## Source code
The code in the `src/` folder is written in `Julia` and builds upon the plotting library [GLMakie](https://docs.makie.org/stable/).

## Getting started 
Create a directory, open the terminal and run
```
git clone https://github.com/killianlutz/BlochBallAnim.jl.git
```
Then open and run the file `scripts/run_animation.jl` in your favorite code editor. Enjoy!

## Outlooks
Enlarging the set of compatible quantum gates by parametrizing $SO(3)$ using Euler angles?

## Last words
The source code is mainly provided for reproducibility purposes and the user's convenience. Any suggestions are more than welcome.

If you enjoy those animations, I would be grateful if you could refer to this [repository](https://github.com/killianlutz/BlochBallAnim.jl) or its author [Killian Lutz](https://github.com/killianlutz).

Thank you for your time and enjoy!
