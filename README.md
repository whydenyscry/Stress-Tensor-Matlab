# Stress Tensor in Matlab
Stress tensor, hydrostatic stress tensor, deviatoric stress tensor, directional stress tensor, invariants of stress tensor and deviatoric stress tensor, principal stresses and shear stresses, von Mises stress, stress triaxiality, Lode angle, octahedral stresses and Mohr's diagram in Matlab for the case when all components of stress tensor are given.
The script is programmed with extended formulas for the purpose of stress-free export to other languages, if necessary.

## Table of Contents

- [Theoretical Background](#theoretical-background)
    - [Cauchy Stress Tensor](#cauchy-stress-tensor)
    - [Volumetric & Mean Stress](#volumetric--mean-stress)
    - [Hydrostatic Stress Tensor](#hydrostatic-stress-tensor)
    - [Deviatoric Stress Tensor](#deviatoric-stress-tensor)
    - [Invariants of the Stress Tensor](#invariants-of-the-stress-tensor)
    - [Invariants of the Deviatoric Stress Tensor](#invariants-of-the-deviatoric-stress-tensor)
    - [Principal Stresses](#principal-stresses)
	- [Invariants in Terms of Prinicipal Stresses](#invariants-in-terms-of-prinicipal-stresses)
    - [Von Mises Equivalent Stress](#von-mises-equivalent-stress)
	- [Lode Angle Dependency](#lode-angle-dependency)
	- [Stress Triaxiality](#stress-triaxiality)
    - [Octahedral Stress](#octahedral-stress)
    - [Directional Stress Tensor](#directional-stress-tensor)
	- [Stress Norm](#stress-norm)
    - [Stress Total Measure](#stress-total-measure)
    - [Principal Shear Stresses](#principal-shear-stresses)
    - [Mohr's Diagram](#mohrs-diagram)
- [Example](#example)
- [References](#references)

 ## Theoretical background
 
 ### Cauchy Stress Tensor
 
 $$
 \boldsymbol{\sigma} = 
	\begin{bmatrix}
		\sigma_x&\tau_{xy}&\tau_{xz}\\
		\tau_{yx}&\sigma_y&\tau_{yz}\\
		\tau_{zx}&\tau_{zy}&\sigma_z
	\end{bmatrix}.
 $$
 
 ### Volumetric & Mean Stress
 
```math
\begin{gather}
\sigma_\text{v} = \text{tr}\boldsymbol{\sigma}, \\
\sigma_\text{m} = \frac{1}{3}\sigma_\text{v}.
\end{gather}
```
 
 ### Hydrostatic Stress Tensor
 
 $$
\boldsymbol{\sigma}\_\text{hyd} = \sigma\_\text{m}\mathbf{I}.
 $$

 ### Deviatoric Stress Tensor
 
 $$
\boldsymbol{\sigma}' =\boldsymbol{\sigma} - \boldsymbol{\sigma}_\text{hyd}.
 $$
 
  ### Invariants of the Stress Tensor
 
 $$
\begin{bmatrix}
			I_1 \\
			I_2 \\
			I_3 
		\end{bmatrix} =
		\begin{bmatrix}
			\text{tr} \boldsymbol{\sigma}\\
			\frac{1}{2}\left[\text{tr}^2 \boldsymbol{\sigma} - \text{tr} \left(\boldsymbol{\sigma}^2\right)\right]\\
			\det \boldsymbol{\sigma}
		\end{bmatrix}.
 $$
 
### Invariants of the Deviatoric Stress Tensor

 $$
\begin{bmatrix}
			J_1 \\
			J_2 \\
			J_3 
		\end{bmatrix}
		=\begin{bmatrix}
			\text{tr}\boldsymbol{\sigma}'\\
			\frac{1}{2}\text{tr}\left(\boldsymbol{\sigma}'^2\right)\\
			\det \boldsymbol{\sigma}'
		\end{bmatrix}=\begin{bmatrix}
		0\\
		\frac{1}{3}I_1^2 - I_2\\
		\frac{2}{27}I_1^3 - \frac{1}{3}I_1 I_2 + I_3
		\end{bmatrix}.
 $$

### Principal Stresses

The principal stresses are eigenvalues of the stress tensor, therefore

```math
\begin{gather}
		\det\left(\boldsymbol{\sigma}-\sigma\mathbf{I}\right)=0\Rightarrow\sigma^3-I_1\sigma^2+I_2\sigma-I_3=0,\\
		\det\left(\boldsymbol{\sigma}'-\sigma'\mathbf{I}\right)=0\Rightarrow \sigma'^3-J_2\sigma'-J_3=0.
	\end{gather}
```

### Invariants in Terms of Prinicipal Stresses

$$
\begin{bmatrix}
			I_1 \\
			I_2 \\
			I_3 
		\end{bmatrix} = 
		\begin{bmatrix}
			\sigma_1+\sigma_2+\sigma_3\\
			\sigma_1\sigma_2+\sigma_2\sigma_3+\sigma_3\sigma_1\\
			\sigma_1\sigma_2\sigma_3
		\end{bmatrix},
$$

$$
\begin{bmatrix}
			J_1 \\
			J_2 \\
			J_3 
		\end{bmatrix} = 
		\begin{bmatrix}
			0\\
			-\left(\sigma_1'\sigma_2'+\sigma_2'\sigma_3'+\sigma_3'\sigma_1'\right)\\
			\sigma_1'\sigma_2'\sigma_3'
		\end{bmatrix}.
$$

### Von Mises Equivalent Stress

$$
\sigma_\text{vM} = \sqrt{\frac{3}{2}\boldsymbol{\sigma}':\boldsymbol{\sigma}'}= \sqrt{3J_2}= \frac{1}{\sqrt{2}}\sqrt{\left(\sigma_1-\sigma_2\right)^2+\left(\sigma_2-\sigma_3\right)^2+\left(\sigma_3-\sigma_1\right)^2}.
$$

### Lode Angle Dependency

Closed-form solution for ordered eigenvalues (principal stresses)

```math
\begin{gather}
		\sigma_1\geq\sigma_2\geq\sigma_3,\\
		\sigma_1'\geq\sigma_2'\geq\sigma_3'
	\end{gather}
```

can be obtained via Lode angle $\Theta$ (or azimuthal Lode angle $\bar\Theta$). Define Lode parameter $\mathcal{L}$ as

$$
	\mathcal{L}=\frac{27}{2} \frac{J_3}{\sigma_\text{vM}^3}=\cos 3\Theta=-\sin 3\bar\Theta,
$$

then

```math
	\begin{gather}
		\Theta = \frac{1}{3}\arccos \mathcal{L} \in \left[0,\dfrac{\pi}{3}\right],\\
		\bar\Theta = \Theta - \frac{\pi}{6} = -\frac{1}{3}\arcsin\mathcal{L} \in \left[-\frac{\pi}{6},\frac{\pi}{6}\right],
	\end{gather}
```

so principal stresses $\sigma_1,\sigma_2,\sigma_3$

$$
\begin{bmatrix}
			\sigma_1\\
			\sigma_2\\
			\sigma_3
		\end{bmatrix} = \begin{bmatrix}
			\sigma_\text{m}\\
			\sigma_\text{m}\\
			\sigma_\text{m}
		\end{bmatrix} + \frac{2}{3}\sigma_\text{vM} \cos\begin{bmatrix}
			\Theta\\
			\Theta-\dfrac{2\pi}{3}\\
			\Theta+\dfrac{2\pi}{3}
		\end{bmatrix}= \begin{bmatrix}
			\sigma_\text{m}\\
			\sigma_\text{m}\\
			\sigma_\text{m}
		\end{bmatrix} + \frac{2}{3}\sigma_\text{vM} \sin\begin{bmatrix}
			\bar\Theta+\dfrac{2\pi}{3}\\
			\bar\Theta\\
			\bar\Theta-\dfrac{2\pi}{3}
		\end{bmatrix},
$$	

and deviatoric principal stresses $\sigma_1',\sigma_2',\sigma_3'$

$$
\begin{bmatrix}
			\sigma_1'\\
			\sigma_2'\\
			\sigma_3'
		\end{bmatrix} = \begin{bmatrix}
			\sigma_1\\
			\sigma_2\\
			\sigma_3
		\end{bmatrix} - \begin{bmatrix}
			\sigma_\text{m}\\
			\sigma_\text{m}\\
			\sigma_\text{m}
		\end{bmatrix} = \frac{2}{3}\sigma_\text{vM} \cos\begin{bmatrix}
			\Theta\\
			\Theta-\dfrac{2\pi}{3}\\
			\Theta+\dfrac{2\pi}{3}
		\end{bmatrix} = \frac{2}{3}\sigma_\text{vM} \sin\begin{bmatrix}
			\bar\Theta+\dfrac{2\pi}{3}\\
			\bar\Theta\\
			\bar\Theta-\dfrac{2\pi}{3}
		\end{bmatrix}.
$$	

### Stress Triaxiality

$$
T_X = \frac{\sigma_\text{m}}{\sigma_\text{vM}}.
$$

### Octahedral Stress

```math
\begin{gather}
		\sigma_\text{oct.} = \sigma_\text{m},\\
		\tau_\text{oct.} = \sqrt{\frac{2}{3}J_2} = \frac{1}{3}\sqrt{\left(\sigma_1-\sigma_2\right)^2+\left(\sigma_2-\sigma_3\right)^2+\left(\sigma_3-\sigma_1\right)^2}.
	\end{gather}
```

### Directional Stress Tensor
The directional stress tensor determines only the principal stresses and the relationship between the components of the stress tensor directions of stresses and the relationship between the components of the stress tensor, but doesn't determine their values, since the components of the directional stress tensor are dimensionless.

$$
\overline{\boldsymbol{\sigma}}' = \frac{1}{\tau_\text{oct.}}\boldsymbol{\sigma}'.
$$

### Stress Norm

$$
\sigma_\text{norm} = \left\lVert\boldsymbol{\sigma}\right\rVert = \sqrt{\boldsymbol{\sigma}:\boldsymbol{\sigma}}=\sqrt{\text{tr}\left(\boldsymbol{\sigma}^2\right)}.
$$

### Stress Total Measure

$$
\sigma_\text{tm} = \sqrt{\sigma_1^2+\sigma_2^2+\sigma_3^2}=\sqrt{\frac{1}{3}I_1^2+2J_2}.
$$

### Principal Shear Stresses

```math
\begin{gather}
		\begin{bmatrix}
			\tau_{13}\\
			\tau_{12}\\
			\tau_{23}
		\end{bmatrix}
		=\frac{1}{2}
		\begin{bmatrix}
			\sigma_1-\sigma_3\\
			\sigma_1-\sigma_2\\
			\sigma_2-\sigma_3\\
		\end{bmatrix},
		\\
		\tau_\text{max} = \tau_{13}.
	\end{gather}
```

The corresponding normal stresses $\sigma_{13}, \sigma_{12}, \sigma_{23}$ acting on the sections where
	$\tau_{13}, \tau_{12}, \tau_{23}$ are acting, respectively, are
	
```math
\begin{gather}
		\begin{bmatrix}
			\sigma_{13}\\
			\sigma_{12}\\
			\sigma_{23}
		\end{bmatrix}
		=\frac{1}{2}
		\begin{bmatrix}
			\sigma_1+\sigma_3\\
			\sigma_1+\sigma_2\\
			\sigma_2+\sigma_3
		\end{bmatrix}
	\end{gather}
```

### Mohr's Diagram

Radii of circles:

$$
\begin{bmatrix}
		R_1\\
		R_2 \\
		R_3
	\end{bmatrix}=
	\begin{bmatrix}
		\tau_{23}\\
		\tau_{13}\\
		\tau_{12}
\end{bmatrix}.
$$

Centers of circles:

$$
O_1=\left(\sigma_{23},  0\right),
$$

$$
O_2=  \left(\sigma_{13}, 0\right),
$$

$$
O_3= \left(\sigma_{12}, 0\right).
$$

## Example 

In the script _StressTensor.m_ the plotting of the Mohr's diagram for the tensor is considered

$$
\boldsymbol{\sigma} =
\begin{bmatrix}
			-22.2 & 9.1 & 7.3 \\
			9.1 & -16.9 & -4.6 \\
			7.3 & -4.6 & 31.8 \\
		\end{bmatrix}
$$

![The Mohr's Diagram](https://github.com/whydenyscry/Stress-Tensor-Matlab/blob/main/images/The_Mohrs_Diagram.png)

## References
1. Chakrabarty, J.(2016). Theory of Plasticity: Third edition
2. Brannon, R.M. (2009).  KAYENTA: Theory and User's Guide
3. Yu, M.-H. (2004). Unified Strength Theory and Its Applications. https://doi.org/10.1007/978-3-642-18943-2 