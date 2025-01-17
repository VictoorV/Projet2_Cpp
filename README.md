# Résolution du problème de diffusion thermique dans un four

Ce projet traite de la résolution du problème de diffusion thermique dans un four.  
L'objectif est de **modéliser numériquement** le problème discrétisé par différences finies, puis de **résoudre le problème linéaire** $A \mathbf{x} = \mathbf{b}$ qui en découle.

## Description du domaine étudié

Le four est représenté par un domaine $\Omega = ]0,1[^2$, avec :  
- $\partial \Omega_D$ : les bords de Dirichlet.  
- $\partial \Omega_N$ : les bords de Neumann.  

L'objet dans le four est modélisé par le domaine $S = [0.25, 0.75] \times [0.4, 0.6]$.  
On note $\mathbf{n}$ le vecteur normal extérieur unitaire.

## Équations du problème

La solution $T : \Omega \to \mathbb{R}$ doit vérifier le système d'équations suivant :

$$
\begin{cases}
-\text{div}(\rho \nabla T(x)) = 0, & \forall x \in \Omega, \\
T = T_D, & \forall x \in \partial \Omega_D, \\
\nabla T \cdot \mathbf{n} = 0, & \forall x \in \partial \Omega_N.
\end{cases}
$$

### Variables et notations
- $T(x)$ : Température au point $x$.
- $\rho$ : Conductivité thermique (supposée constante ici).
- $\nabla T$ : Gradient de la température.
- $\text{div}$ : Divergence d'un champ vectoriel.
- $\mathbf{n}$ : Vecteur normal unitaire extérieur.
