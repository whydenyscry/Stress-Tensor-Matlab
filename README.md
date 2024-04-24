# Stress Tensor in Matlab
Stress tensor, hydrostatic stress tensor, deviatoric stress tensor, invariants of stress tensor and deviatoric stress tensor, principal stresses, von Mises stress and Mohr's diagram in Matlab for the case when all components of stress tensor are given.
The script is programmed with extended formulas for the purpose of stress-free export to other languages, if necessary.

 ## Theoretical background
 
 ### Stress Tensor
 
 $$
 \boldsymbol{\sigma} = 
	\begin{pmatrix}
		\sigma_x&\tau_{xy}&\tau_{xz}\\
		\tau_{yx}&\sigma_y&\tau_{yz}\\
		\tau_{zx}&\tau_{zy}&\sigma_z
	\end{pmatrix}.
 $$
 
 ### Thermodynamic Pressure
 
 $$
P = -\frac{1}{3}\text{tr}\boldsymbol{\sigma}.
 $$
 
 ### Hydrostatic Stress Tensor
 
 $$
\boldsymbol{\sigma}_\text{hyd} = -P\mathbf{I}.
 $$

 ### Hydrostatic Stress Tensor
 
 $$
\boldsymbol{\sigma}' =\boldsymbol{\sigma} - \boldsymbol{\sigma}_\text{hyd}.
 $$
 
  ### Invariants of the Stress Tensor
 
 $$
	I_1 = \text{tr} \boldsymbol{\sigma},
 $$
 
  $$
	I_2 = \frac{1}{2}\left[\text{tr}^2 \boldsymbol{\sigma} - \text{tr} \left(\boldsymbol{\sigma}^2\right)\right],
 $$
 
  $$
	I_3 = \det \boldsymbol{\sigma}.
 $$
 
### Invariants of the Deviatoric Stress Tensor

 $$
	J_1 = \text{tr}\boldsymbol{\sigma}'=0,
 $$
 
 $$
	J_2 = \frac{1}{2}\text{tr}\left(\boldsymbol{\sigma}'^2\right) = \frac{1}{3}I_1^2 - I_2,
 $$
 
 $$
	J_3 = \det \boldsymbol{\sigma}' =\frac{2}{27}I_1^3 - \frac{1}{3}I_1 I_2 + I_3.
 $$

### Principal Stresses

The principal stresses are eigenvalues of the stress tensor, therefore

$$
\det\left(\boldsymbol{\sigma}-\sigma\mathbf{I}\right)=0\Rightarrow\sigma^3-I_1\sigma^2+I_2\sigma-I_3=0.
$$

The solution of this equation can be given in trigonometric form:

$$
p=-\frac{1}{3}I_1^2+I_2,\quad q=-2\left(\frac{I_1}{3}\right)^3+\frac{1}{3}I_1I_2-I_3,\quad \alpha=\arccos{\left(-\frac{q}{2\sqrt{-\left(\frac{p}{3}\right)^3}}\right)},
$$

$$
\sigma_1=\frac{1}{3}I_1 + 2\sqrt{-\frac{p}{3}}\cos{\frac{\alpha}{3}},
$$

$$
\sigma_2= \frac{1}{3}I_1 -2\sqrt{-\frac{p}{3}}\cos{\left(\frac{\alpha}{3}+\frac{\pi}{3}\right)},
$$

$$
\sigma_3= \frac{1}{3}I_1 -2\sqrt{-\frac{p}{3}}\cos{\left(\frac{\alpha}{3}-\frac{\pi}{3}\right)}.
$$

Further, it is assumed that the principal stresses are sorted as follows

$$
\sigma_1\geq\sigma_2\geq\sigma_3.
$$

Then invariants can be expressed in terms of principal stresses

$$
I_1=\sigma_1+\sigma_2+\sigma_3,
$$

$$
I_2= \sigma_1\sigma_2+\sigma_2\sigma_3+\sigma_3\sigma_1,
$$

$$
I_3= \sigma_1\sigma_2\sigma_3.
$$

### Von Mises Stress

$$
\sigma_\text{vM} = \sqrt{3J_2}= \sqrt{\frac{3}{2}\text{tr}\left(\boldsymbol{\sigma}'^2\right)} =\sqrt{\frac{\left(\sigma_1-\sigma_2\right)^2+\left(\sigma_2-\sigma_3\right)^2+\left(\sigma_3-\sigma_1\right)^2}{2}}.
$$

### Mohr's Diagram

Radii of circles:

$$
R_1=\frac{1}{2}\left(\sigma_2 - \sigma_3\right),
$$

$$
R_2= \frac{1}{2}\left(\sigma_1 - \sigma_3\right),
$$

$$
R_3= \frac{1}{2}\left(\sigma_1 - \sigma_2\right).
$$

Centers of circles:

$$
O_1=\left(\frac{1}{2}\left(\sigma_2 + \sigma_3\right),  0\right),
$$

$$
O_2=  \left(\frac{1}{2}\left(\sigma_1 + \sigma_3\right), 0\right),
$$

$$
O_3= \left(\frac{1}{2}\left(\sigma_1 + \sigma_2\right), 0\right).
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