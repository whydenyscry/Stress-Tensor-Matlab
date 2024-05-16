# Stress Tensor in Matlab
Stress tensor, hydrostatic stress tensor, deviatoric stress tensor, invariants of stress tensor and deviatoric stress tensor, principal stresses, von Mises stress and Mohr's diagram in Matlab for the case when all components of stress tensor are given.
The script is programmed with extended formulas for the purpose of stress-free export to other languages, if necessary.

 ## Theoretical background
 
 ### Stress Tensor
 
 $$
 \boldsymbol{\sigma} = 
	\begin{bmatrix}
		\sigma_x&\tau_{xy}&\tau_{xz}\\
		\tau_{yx}&\sigma_y&\tau_{yz}\\
		\tau_{zx}&\tau_{zy}&\sigma_z
	\end{bmatrix}.
 $$
 
 ### Thermodynamic Pressure
 
 $$
P = -\frac{1}{3}\text{tr}\boldsymbol{\sigma}.
 $$
 
 ### Hydrostatic Stress Tensor
 
 $$
\boldsymbol{\sigma}_\text{hyd} = -P\mathbf{I}.
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

$$
\det\left(\boldsymbol{\sigma}-\sigma\mathbf{I}\right)=0\Rightarrow\sigma^3-I_1\sigma^2+I_2\sigma-I_3=0,
$$

$$
\det\left(\boldsymbol{\sigma}'-\sigma'\mathbf{I}\right)=0\Rightarrow \sigma'^3-J_2\sigma'-J_3=0.
$$

The closed-form solution for the ordered eigenvalues can be obtained via the Lode angle $\theta$:

$$
\sin\left(3\theta\right)=-\frac{J_3}{2}\left(\frac{3}{J_2}\right)^{3/2},
$$

$$
\begin{bmatrix}
			\sigma_1\\
			\sigma_2\\
			\sigma_3
		\end{bmatrix} = \frac{1}{3} \begin{bmatrix}
		I_1\\
		I_1\\
		I_1\\
		\end{bmatrix}+\sqrt{J_2} \begin{bmatrix}
		\cos\theta-\frac{1}{\sqrt{3}}\sin\theta\\
		\frac{2}{\sqrt{3}}\sin\theta\\
		\cos\theta+\frac{1}{\sqrt{3}}\sin\theta
		\end{bmatrix},
$$

$$
\begin{bmatrix}
			\sigma_1'\\
			\sigma_2'\\
			\sigma_3'
		\end{bmatrix} = \begin{bmatrix}
		\sigma_1\\
		\sigma_2\\
		\sigma_3
		\end{bmatrix} - \frac{1}{3} \begin{bmatrix}
		I_1\\
		I_1\\
		I_1\\
		\end{bmatrix}=\sqrt{J_2} \begin{bmatrix}
		\cos\theta-\frac{1}{\sqrt{3}}\sin\theta\\
		\frac{2}{\sqrt{3}}\sin\theta\\
		\cos\theta+\frac{1}{\sqrt{3}}\sin\theta
		\end{bmatrix}.
$$

The principal stresses are sorted as follows

$$
\sigma_1\geq\sigma_2\geq\sigma_3,
$$
	
$$
\sigma_1'\geq\sigma_2'\geq\sigma_3'.
$$
	
Then invariants can be expressed in terms of principal stresses

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

### Von Mises Stress

$$
\sigma_\text{vM} = \sqrt{3J_2}= \sqrt{\frac{3}{2}\text{tr}\left(\boldsymbol{\sigma}'^2\right)} =\frac{1}{2}\sqrt{\left(\sigma_1-\sigma_2\right)^2+\left(\sigma_2-\sigma_3\right)^2+\left(\sigma_3-\sigma_1\right)^2}.
$$

### Octahedral Stress

$$
\sigma_\text{oct.} = \frac{1}{3}I_1,
$$

$$
\tau_\text{oct.} = \sqrt{\frac{2}{3}J_2}=\frac{1}{3}\sqrt{\left(\sigma_1-\sigma_2\right)^2+\left(\sigma_2-\sigma_3\right)^2+\left(\sigma_3-\sigma_1\right)^2}.
$$

### Directional Stress Tensor
The directional stress tensor determines only the principal stresses and the relationship between the components of the stress tensor directions of stresses and the relationship between the components of the stress tensor, but doesn't determine their values, since the components of the directional stress tensor are dimensionless.

$$
\overline{\boldsymbol{\sigma}}' = \frac{1}{\tau_\text{oct.}}\boldsymbol{\sigma}'.
$$

### Principal Shear Stresses

$$
\tau_{13} = \frac{1}{2}\left(\sigma_1-\sigma_3\right),
$$

$$
\tau_{12} = \frac{1}{2}\left(\sigma_1-\sigma_2\right),
$$

$$
\tau_{23} = \frac{1}{2}\left(\sigma_2-\sigma_3\right).
$$

$$
\tau_\text{max} = \tau_{13}.
$$

The corresponding normal stresses $\sigma_{13}, \sigma_{12}, \sigma_{23}$ acting on the sections where
	$\tau_{13}, \tau_{12}, \tau_{23}$ are acting, respectively, are
	
$$
\sigma_{13} = \frac{1}{2}\left(\sigma_1+\sigma_3\right),
$$

$$
\sigma_{12} = \frac{1}{2}\left(\sigma_1+\sigma_2\right),
$$

$$
\sigma_{23} = \frac{1}{2}\left(\sigma_2+\sigma_3\right).
$$

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
\begin{pmatrix}
			-22.2 & 9.1 & 7.3 \\
			9.1 & -16.9 & -4.6 \\
			7.3 & -4.6 & 31.8 \\
		\end{pmatrix}
$$

![The Mohr's Diagram](https://github.com/whydenyscry/Stress-Tensor-Matlab/blob/main/images/The_Mohrs_Diagram.png)

## References
1. Chakrabarty, J.(2016). Theory of Plasticity: Third edition
2. Brannon, R.M. (2009).  KAYENTA: Theory and User's Guide
3. Yu, M.-H. (2004). Unified Strength Theory and Its Applications. https://doi.org/10.1007/978-3-642-18943-2 