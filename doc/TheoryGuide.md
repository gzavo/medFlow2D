************************
# medFlow2D Theory Guide
************************


> Written by Dr. Gábor Závodszky.
>
> All rights reserved, 2015.


## General information on the lattice Boltzmann method

### Books
- Sukop, M. C., DT Thorne, Jr. (2006). Lattice Boltzmann Modeling: An introduction for geoscientists and engineers, Springer.
- Mohamad, A. A. (2011). Lattice Boltzmann method: fundamentals and engineering applications with computer codes. Springer Science & Business Media.


### Web
- http://www.scholarpedia.org/article/Lattice_Boltzmann_Method
- http://wiki.palabos.org/Lattice%20Boltzmann%20Method


## Grid numbering and the coordinate system
The coordinate system was chosen to match the coordinate system of the images which describe the geometry.

  Coordinate system :

    (0,0)
      +--------> x
      |
      |
      |
      |  
      ˇ
      y

  Lattice velocity directions for D2Q9 grid ($\vec{c_i}$):

      4  3  2
       \ | /  
        \|/
    5 -- 0 -- 1
        /|\
       / | \
      6  7  8


## Basics, BGK dynamics

> Reference: Chen, S., & Doolen, G. D. (1998). Lattice Boltzmann method for fluid flows. Annual review of fluid mechanics, 30(1), 329-364.

Directions from the grid numbering in the given coordinate system:

$$ \vec{c} = [(0,0); (1,0); (1,-1); (-1,0); (-1,-1); (-1,0); (-1,1); (0,1); (1,1)] $$

Weighting factor for the discrete velocities:

$$ \vec{w} = [4/9; 1/9; 1/36; 1/9; 1/36; 1/9; 1/36; 1/9; 1/36;] $$

LBM equation:
$$
f_i(x+\vec{c_i} \delta t, t+\delta t) - f_i(x,t) =  \Omega_i
$$

BGK approximation of the collision operator $\Omega$:

$$ f_i(x+\vec{c_i} \delta t, t+\delta t) = f_i(x,t) - \omega(f_i-f_i^e) $$

with a single relaxation time $\omega = \frac{1}{\tau}$.

The equilibrium term (based on the power series expansion of the Maxwell-Boltzmann distribution):

$$ f_i^e = w_i \rho [1 + \frac{\vec{c_i} \vec{u}}{c_s^2} + \frac{(\vec{c_i} \vec{u})^2}{2 c_s^4}-\frac{u^2}{2 c_s^2}] = w_i \rho [1 + 3 \vec{c_i} \vec{u} + 9(\vec{c_i} \vec{u})^2 - 3/2 \vec{u}^2] $$

For the D2Q9 grid the lattice speed of sound $c_s^2 = \frac{1}{3}$.

Regaining the macroscopic variables:

$$
\rho = \sum_{i=0}^8 f_i \\
\vec{u} = \frac{1}{\rho} \sum_{i=0}^8\vec{c_i} f_i \\
p=\rho c_s^2 \\
\nu = 1/3(\tau -1/2)
$$

## Multiple Relaxation Time (MRT) dynamics

> Reference: d'Humières, D. (2002). Multiple–relaxation–time lattice Boltzmann models in three dimensions. Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and Engineering Sciences, 360(1792), 437-451.

> Reference 2: Guo, X., Zhong, C., Zhuo, C., & Cao, J. (2014). Multiple-relaxation-time lattice Boltzmann method for study of two-lid-driven cavity flow solution multiplicity. Theoretical and Computational Fluid Dynamics, 28(2), 215-231.

> Reference 3: Lallemand, P., & Luo, L. (2000). Theory of the lattice boltzmann method: dispersion, dissipation, isotropy, galilean invariance, and stability. Physical Review. E, Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics, 61(6 Pt A), 6546–62.

With this dynamics the collision takes place in momentum space to allow the different moments to relaxate differently.

The collision operator: $\Omega = -M^{-1} \cdot S \cdot [m(x,t) -m^e(x,t) ]$.

The distribution componenets cast into moment space:
$m = (\rho, e, \epsilon, j_x = \rho u_x, q_x, j_y = \rho u_y, q_y, p_{xx}, p_{xy})$

The momentum equilibrium:
$m^e = \rho \cdot (1, -2+3u^2, 1-3u^2, u_x, -u_x, u_y, -u_y, u_x^2-u_y^2, u_x u_y)$

The matrix to transform to momentum space:

		   {{ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.},
			{-4., -1.,  2., -1.,  2., -1.,  2., -1.,  2.},
			{ 4., -2.,  1., -2.,  1., -2.,  1., -2.,  1.},
			{ 0.,  1.,  1.,  0., -1., -1., -1.,  0.,  1.},
	M =		{ 0., -2.,  1.,  0., -1.,  2., -1.,  0.,  1.},
			{ 0.,  0., -1., -1., -1.,  0.,  1.,  1.,  1.},
			{ 0.,  0., -1.,  2., -1.,  0.,  1., -2.,  1.},
			{ 0.,  1.,  0., -1.,  0.,  1.,  0., -1.,  0.},
			{ 0.,  0., -1.,  0.,  1.,  0., -1.,  0.,  1.}}


Since the rows of $M$ are mutually orthogonal $M M^T$ is a diagonal matrix, thus, $M^{-1} = M^T (M M^T)^{-1}$.

		   {{ 1./9, -1./ 9,  1./ 9,  0   ,  0    ,  0   ,  0    ,  0   ,  0   },
			{ 1./9, -1./36, -1./18,  1./6, -1./ 6,  0   ,  0    ,  1./4,  0   },
			{ 1./9,  1./18,  1./36,  1./6,  1./12, -1./6, -1./12,  0   , -1./4},
			{ 1./9, -1./36, -1./18,  0   ,  0    , -1./6,  1./ 6, -1./4,  0   },
	Minv =	{ 1./9,  1./18,  1./36, -1./6, -1./12, -1./6, -1./12,  0   ,  1./4},
			{ 1./9, -1./36, -1./18, -1./6,  1./6 ,  0   ,  0    ,  1./4,  0   },
			{ 1./9,  1./18,  1./36, -1./6, -1./12,  1./6,  1./12,  0   , -1./4},
			{ 1./9, -1./36, -1./18,  0   , 0     ,  1./6, -1./ 6, -1./4,  0   },
			{ 1./9,  1./18,  1./36,  1./6,  1./12,  1./6,  1./12,  0   ,  1./4}}

The relaxation matrix $ S $ is a diagonal matrix:

$ S = diag(0, s_e, s_\epsilon, 0, s_q, 0, s_q, s_\nu, s_\nu) $

The relaxation parameters $ s_0, s_3, s_5 $ were set to zero (to obey the conservation of mass and momentum).
According to (Lallemand, 2000), though there is no definite optimum that matches all setup, a generally acceptable setting of the free relaxations to set them slightly higher than one: $ s_e = 1.001 $, $ s_\epsilon = 1.001 $, $ s_q = 1.001 $.

The kinematic $ \nu $ and bulk viscosity $ \zeta $:

$$
\nu = c_s^2(\frac{1}{s_\nu}-1/2) \\
\zeta = c_s^2(\frac{1}{s_e}-1/2)
$$

## Evaluation of the stress tensor

In the Chapman-Enskog multiscale analysis the populations are expanded with repect to a small number $ \epsilon $ (identified with the Knudsen number $ \epsilon = \frac{\lambda}{L}$, the mean free path of particles over the reference lenbgth):

$ f_i = f_i^{(0)} + \epsilon f_i^{(1)} + \epsilon^2 f_i^{(2)} + ... $

To reach the governing equations the populations must be expanded around the equilibrium functions, so the zeroth order term of each population is identified with the matching equilibrium distribution function $ f_i^{(0)} = f^e $. The rest of the expanded series is the non-equilibrium part.

Thus, $ f_i - f_i^e = f_i^n $ and $ \sum_{i=0}^{8}f_i^n = 0 $

In practice $ f^{(1)} $ can be approximated with $ f^n $.

The stress tensor $ \tau $ is defined in terms of the rate of strain tensor $ S $:

$ \tau = 2 \nu \rho S $

The third order moment of the distributions is:
$ \prod = \sum_{i=0}^{8} \vec{c_i}\vec{c_i}f_i $

Its first order expanded part corresponds to the rate of strain tensor:

$ \prod^{(1)}= \sum_{i=0}^{8} \vec{c_i}\vec{c_i}f_i^{(1)} = -\frac{2c_s^2}{\omega}\rho S$

## Boundary conditions (Zou-He)
> Reference: Zou, Q., & He, X. (1997). On pressure and velocity boundary conditions for the lattice Boltzmann BGK model. Physics of Fluids (1994-present), 9(6), 1591-1598.

Note: medFlow2D only allows straight openings located along one of the four sides of the domain.

The additional assumption for the Zou-He method is that the non-equilibrium parts are equal between the directions perpendiculair to the boundary surface.

### von Neumann (flux/velocity) and Dirichlet (density/pressure) boundaries

#### Left side
Unknown distributions: 2; 1; 8

Equations to start with:
1. The equality of the non-equilibrium part of the distributions perpendicular to the boundary.
2. The macroscopic velocity in x direction.
3. The macroscopic velocity in y direction.
4. The macroscopic density.

Incase of velocity normal to the boundary $ \vec{u} = (u_x; 0) $:

$ f_5 - f_5^e = f_1 - f_1^e 					\tag 1 \\ $
$ \rho u_x = f_1 + f_2+f_8-f_4-f_5-f_6 			\tag 2 \\ $
$ 0 = f_6+f_7+f_8-f_2-f_3-f_4  					\tag 3 \\$
$ \rho = f_0+f_1+f_2+f_3+f_4+f_5+f_6+f_7+f_8  	\tag 4 $

Steps to solve:
1. Express the three common components from (2) and (4), then equate the right sides.
2. From this equation solve for $\rho$ or $u_x/u_y$ depending on the boundary type.
3. Using the expression of the equilibrium term (1) can be expressed with $\rho$ and $u_x/u_y$.
4. Using this new expression and (3) in (2) one can solve for either of the two remaining unknown components.


Resulting equations for the unknown distribution components in general case ($ \vec{u} = (u_x; u_y) $):

$ \rho = \frac{f_0 + f_3 + f_7 + 2(f_4 + f_5 + f_6)}{1-u_x} $

$ u_x = 1 - \frac{f_0 + f_3 + f_7 + 2(f_4 + f_5 + f_6)}{\rho} $

$ f_1 = f_5 + 2/3\rho u_x $
$ f_2 = 1/6\rho u_x - 1/2\rho u_y + f_6 + 1/2(f_7-f_3) $
$ f_8 = 1/6\rho u_x + 1/2\rho u_y + f_4 + 1/2(f_3-f_7) $

#### Rigth side

$ \rho = \frac{f_0 + f_3 + f_7 + 2(f_1 + f_2 + f_8)}{1+u_x} $

$ u_x = -1 + \frac{f_0 + f_3 + f_7 + 2(f_1 + f_2 + f_8)}{\rho} $

$ f_5 = f_1 - 2/3\rho u_x $
$ f_6 = -1/6\rho u_x + 1/2\rho u_y + f_2 + 1/2(f_3-f_7) $
$ f_4 = -1/6\rho u_x - 1/2\rho u_y + f_8 + 1/2(f_7-f_3) $


#### Top boundary

$ \rho = \frac{f_0 + f_1 + f_5 + 2(f_2 + f_3 + f_4)}{1-u_y} $

$ u_y = 1 - \frac{f_0 + f_1 + f_5 + 2(f_2 + f_3 + f_4)}{\rho} $

$ f_7 = f_3 + 2/3\rho u_y $
$ f_6 = 1/6\rho u_y - 1/2\rho u_x + f_2 + 1/2(f_1-f_5) $
$ f_8 = 1/6\rho u_y + 1/2\rho u_x + f_4 + 1/2(f_5-f_1) $

#### Bottom boundary

$ \rho = \frac{f_0 + f_1 + f_5 + 2(f_6 + f_7 + f_8)}{1+u_y} $

$ u_y = -1 + \frac{f_0 + f_1 + f_5 + 2(f_6 + f_7 + f_8)}{\rho}$

$ f_3 = f_7 - 2/3\rho u_y $
$ f_2 = -1/6\rho u_y + 1/2\rho u_x + f_6 + 1/2(f_5-f_1) $
$ f_4 = -1/6\rho u_y - 1/2\rho u_x + f_8 + 1/2(f_1-f_5) $

## Passive solute component

> Reference: Sukop, M. C. (2006). DT Thorne, Jr. Lattice Boltzmann Modeling: An introduction for geoscientists and engineers, Springer.

### Describing equations

The equilibrium term: $ f_{s,i}^e = w_i \rho_s(1 + 3\vec{c_i} \vec{u}) $

The concentration: $ \rho_s = \sum_{i=0}^8 f_{s,i} $

The diffusion coefficient: $ D_s = 1/3(\tau_s - 1/2) $

### Boundary conditions - Constant concentration bd.

#### Left side

Equations to start with:
1. The required concentration: $ \rho_s = \sum_{i=0}^8 f_{s,i} $


Steps to solve:

1. Since the velocity is provided by the underlying flow field, the unknown densities can be calculated using the weighting factors, thus the unknown densities can be in the form of $ f_{s,i} = w_i \rho_{residual} $.
2. Two since this is one equation with one unknown it can be solved for $ \rho_{residual} $.

$ \rho_{residual}= \frac{\rho_s - (f_{s,0} + f_{s,3} + f_{s,4} + f_{s,5} + f_{s,6} + f_{s,7})}{w_1+w_2+w_8} $

From this, the missing distributions:

$$
f_1 = w_1*\rho_{residual} \\
f_2 = w_2*\rho_{residual} \\
f_8 = w_8*\rho_{residual}
$$

#### Right side

$ \rho_{residual}= \frac{\rho_s - (f_{s,0} + f_{s,1} + f_{s,2} + f_{s,3} + f_{s,7} + f_{s,8})}{w_4+w_5+w_6} $

$$
f_4 = w_4*\rho_{residual} \\
f_5 = w_5*\rho_{residual} \\
f_6 = w_6*\rho_{residual}
$$

#### Top

$ \rho_{residual}= \frac{\rho_s - (f_{s,0} + f_{s,1} + f_{s,2} + f_{s,3} + f_{s,4} + f_{s,5})}{w_6+w_7+w_8} $

$$
f_6 = w_6*\rho_{residual} \\
f_7 = w_7*\rho_{residual} \\
f_8 = w_8*\rho_{residual}
$$

#### Bottom

$ \rho_{residual}= \frac{\rho_s - (f_{s,0} + f_{s,1} + f_{s,5} + f_{s,6} + f_{s,7} + f_{s,8})}{w_2+w_3+w_4} $

$$
f_2 = w_2*\rho_{residual} \\
f_3 = w_3*\rho_{residual} \\
f_4 = w_4*\rho_{residual}
$$

### Boundary conditions - Zero diffusive flux bd.

#### Left side

Equations to start with:
1. Require the diffusive solute mass flux to be zero normal to the boundary: $ \sum_{i=0}^8 f_{s,i}\vec{c_i}\vec{n} = 0 $

Steps to solve:

1. Expand the equation using the previous  $ f_{s,i} = w_i \rho_{residual} $ assumption for the unknown distributions.
2. Solve for $ \rho_{residual} $.
3. Express the unknown terms with $ \rho_{residual} $.

$ \rho_{residual}= \frac{f_{s,4} + f_{s,5} + f_{s,6} }{w_1+w_2+w_8} $

From this, the missing distributions:

$$
f_1 = w_1*\rho_{residual} \\
f_2 = w_2*\rho_{residual} \\
f_8 = w_8*\rho_{residual}
$$

#### Right side

$ \rho_{residual}= \frac{f_{s,1} + f_{s,2} + f_{s,8} }{w_4+w_5+w_6} $

$$
f_4 = w_4*\rho_{residual} \\
f_5 = w_5*\rho_{residual} \\
f_6 = w_6*\rho_{residual}
$$

#### Top

$ \rho_{residual}= \frac{f_{s,2} + f_{s,3} + f_{s,4} }{w_6+w_7+w_8} $

$$
f_6 = w_6*\rho_{residual} \\
f_7 = w_7*\rho_{residual} \\
f_8 = w_8*\rho_{residual}
$$

#### Bottom

$ \rho_{residual}= \frac{f_{s,6} + f_{s,7} + f_{s,8} }{w_2+w_3+w_4} $

$$
f_2 = w_2*\rho_{residual} \\
f_3 = w_3*\rho_{residual} \\
f_4 = w_4*\rho_{residual}
$$

## Boundary conditions (Regularized)

> Reference: Latt, J., Chopard, B., Malaspinas, O., Deville, M., & Michler, A. (2008). Straight velocity boundaries in the lattice Boltzmann method. Physical Review E, 77(5), 056703.

Note: medFlow2D only allows straight openings located along one of the four sides of the domain.

This type of boundary replaces all of the population (in contrast of Zou-He), typically providing higher stability in return for lower accuracy.

1. Compute the missing macroscopic value (component of $ \vec{u}$ or $ \rho $) the same way as with the Zou-He method.
2. Apply the bounce back of the off-equilibrium part to all unknown populations ($ f_{i,unknown} = f_i^e + f_{opp(i)}-f_{opp(i)}^e$).
3. Use the resulting temporary distributions to evaluate $ \prod^{(1)} $
4. Replace all the distribution components with the one reconstructable from $ \prod^{(1)} $.

## Calculation of wall normals

> Reference: Matyka, M., Koza, Z., & Mirosław, Ł. (2013). Wall orientation and shear stress in the lattice Boltzmann model. Computers & Fluids, 73, 115–123. Computational Physics; Fluid Dynamics.

## Non-Newtonian dynamics

> Reference: Ashrafizaadeh, M., & Bakhshaei, H. (2009). A comparison of non-Newtonian models for lattice Boltzmann blood flow simulations. Computers & Mathematics with Applications, 58(5), 1045–1054.

The non-Newtonian dynamics is implemented by varying the local viscosity according to the model fitted to blood viscosity measurements.
medFlow2D uses the Carreau-Yasuda description.

### Carreau-Yasuda

$$
\frac{\mu - \mu_\infty}{\mu_0-\mu_\infty} = [1+ (\lambda \dot{\gamma})^a]^{\frac{n-1}{a}}
$$

The usual parameters for simulating blood flows are:
$a=0.644, n=0.392, \lambda=0.11, \mu_0=22\times 10^{-3}, \mu_\infty=2.2\times 10^{-3} $.

## Macroscopic porous material

The effects of the porous media are implemented as momentum losses imposed upon the fluid field.

### Darcy-Forcheimer equation

> Guo, Z., & Zhao, T. S. (2002). Lattice Boltzmann model for incompressible flows through porous media. Physical Review E, 66(3), 036304. doi:10.1103/PhysRevE.66.036304

For porous media, it is known that the pressure drop is nearly linearly proportional
to the flow velocity in case of low velocities.
For higher velocities, the pressure drop exceeds the one predicted by this linear approximation.
To account for this effect, one can extend the linear model
with a quadratic term called Forcheimer's term:

$$
\bigtriangledown p = -\frac{\mu}{\kappa} \vec{v} - \frac{\rho}{\kappa_1} \lvert \vec{v}\rvert \vec{v}.
$$

## Fast numerical calculation of local shear stress properties

The shear rate tensor has the following components:

$$
\vec{S}_{\alpha \beta} = \frac{1}{2} \left( \partial_\alpha u_\beta + \partial_\beta u_\alpha \right).
$$

This tensor can be related to the second-order moment of the $f_i^{(1)}$ part of the distribution function:
$$
	\vec{\prod}^{(1)}=-\frac{2 c_s^2 \rho}{\omega} \vec{S},
$$
where the tensor $\vec{\prod}^{(1)}$ denotes this second-order moment.

$ f_i^e $ can be regarded as a first-order accurate approximation of $f_i^{(1)}$.
It is easy to calculate it with the use of the equilibrium part $f_i^{eq}$ as $f_i^{neq} = f_i - f_i^{eq}$.
This leads to the following expression, which is computable locally at every lattice site:
$$
	\vec{S}=-\frac{\omega}{2 C_s^2 \rho}\vec{\prod}^{(1)}=-\frac{\omega}{2 c_s^2 \rho} \sum\limits_{i} (f_i -f_i^{(eq)}) \vec{c_i} \vec{c_i}.
$$

Using the technique of Mohr's Circle, one can extract the direction of the principal stresses and the maximum shear stresses acting on inclined planes.
For Newtonian fluids the viscous stress tensor can be related to the shear rate tensor as:
$$
	 \vec{ \sigma } = 2 \mu \vec{S},
$$

where $\mu$ denotes the dynamic viscosity. In two dimension the components of this tensor can be written as:
$$
	\vec{\sigma} =  \left[
	  \begin{array}{ c c }
	     \sigma_x & \tau_{xy} \\
	     \tau_{xy} & \sigma_y
	  \end{array} \right].
$$

Using these notations the magnitude of the maximum shear stress equals the radius of Mohr's Circle:
$$
	\tau_{max} = \sqrt{\left[ \frac{1}{2}(\sigma_x - \sigma_y)\right]^2 + \tau_{xy}^2},
$$

while the orientation of the principal planes (which are $90^{\circ}$ apart) can be recovered as:
$$
	\tan 2 \theta_P = \frac{2 \tau_{xy}}{\sigma_x - \sigma_y}.
$$

Therefore, using these pieces of information, the whole computation of the virtual force $\vec{F_M}$ can be reduced to a few algebraic operations and a trigonometric function per lattice site.


## Margination process

- This $ \vec{F_M} $ force is perpendicular to the inclined plane of the maximum shear stress at a given location and points towards that side of the plane along where the velocity gradient is negative.
- Its magnitude is proportional to the magnitude of the emerging shear stress acting in the aforementioned plane. I note here that this proportionality ratio definitely depends on the material properties and the relative sizes of the two particle types (i.e., the platelets and the RBC cells),  and that it will require a proper parameter study later. In the current work, it takes the value of unity.


## Simple model for hemostasis

$$
	P_{coag} = \frac{\rho_{platelet} \cdot \rho_{ADP}}{\tau_{max}},
$$

where $ P_{coag} $ is a probability in the sense that its value decides whether a fluid cell should come to stasis or not at any given time, based on the current local platelet concentration $ \rho_{platelet} $, the local ADP concentration $ \rho_{ADP} $, and the local maximum shear component of the stress tensor denoted by $ \tau_{max} $.
The threshold level of $ P_{coag} $ that a numerical lattice has to reach for coagulation is an empirical parameter of the model for now. A later work should explore the deeper relationship among these parameters.
For a numerical fluid cell to come to stasis, these three parameters have to remain in the coagulation zone for a $ t_{window} $ time. This time-window was chosen to be $ t_{window}=20~ms $ because this is the usual time-frame for ADP to activate a platelet. This also means that a newly registered near-wall lattice cannot turn into a solid cell sooner than 20 ms.

## Smagorinsky subgrid model

> Reference: Yu, H., Girimaji, S. S., & Luo, L. S. (2005). Lattice Boltzmann simulations of decaying homogeneous isotropic turbulence. Physical Review E, 71(1), 016708.

> Reference 2: Hou, S., Sterling, J., Chen, S., & Doolen, G. D. (1996). A Lattice Boltzmann Subgrid Model for High Reynolds Number Flows. Pattern Formation and Lattice Gas Automata, 6, 1–18.

The role of this procedure is to parameterize the turbulent energy dissipation in turbulent flows, where the larger eddies extract energy from the mean flow and ultimately transfer some of it to the smaller eddies which, in turn, pass the energy to even smaller eddies, and so on up to the smallest scales, where the eddies convert the kinetic energy into internal energy of the fluid. At this scales (also known as Kolmogorov scale), the viscous friction dominates the flow. In classical LB applications, a convenient method of modelling turbulent dissipation is through a locally-enhanced collision, which effectively stabilizes the simulation.

After selecting a proper Smagorinsky constant $ C $ one can compute the strain rate tensor as
$$
\|S\| = \frac{1}{6C^2}(\sqrt{\nu_0^2+18C^2 \sqrt{\prod^{(1)} : \prod^{(1)}} } - \nu_0).
$$

Using this tensor the new relaxation frequency can be written as:

$$
\tau_S=3(\nu_0 + C^2 \|S\|)+1/2.
$$

A suitable value for $ C $ can be $ 0.038 $ (naturally, its appropriate value depends on the problem at hand).

The formalism for the implementation is the following:

$$
\tau_0 = \frac{1}{\omega} \\
Q = \sqrt{\prod^{(1)} : \prod^{(1)} } \\
\tau_{turb} = 1/2 (\sqrt{\tau_0^2 + 18*C^2*Q}) - \tau_0 \\
\omega_{turb} = \frac{1}{\tau_0+\tau_{turb}}
$$

### Driest damping of Smagorinsky coefficient

> Reference: Ferziger, Joel H., and Milovan Peric. Computational methods for fluid dynamics. Springer Science & Business Media, 2012.

A usual technique to reduce the near-wall eddy viscosity in RANS models is to use the following damping term:

$$
C=C_0(1-e^{-n^+/A^+})^2,
$$

where $ n^+ $ is the distance from the wall in viscous wall units: $ n^+ = n u_\tau/\nu$, where $ u_\tau = \sqrt{\tau_w/\rho} $ is the shear velocity, $ \tau $ is the wall shear stress and $ A^+ $ is a constant usually taken to be approximately 25.

> Note: Altough this model produces the desired results, it is hard to justify it in the contexts of LES, since it depends on the distance from the wall that is a non-local property.



## Forcing term

> Reference 1: Guo, Z., Zheng, C., & Shi, B. (2002). Discrete lattice effects on the forcing term in the lattice Boltzmann method. Physical Review E, 65(4), 046308.

> Reference 2: Mohamad, a. a., & Kuzmin, A. (2010). A critical evaluation of force term in lattice Boltzmann method, natural convection problem. International Journal of Heat and Mass Transfer, 53(5-6), 990–996.

Given the force $ \vec{F} $, the term accounting for forcing effects is added to the collision process:

$$ f_i^e = w_i (1- \frac{1}{2 \tau}) [\frac{\vec{c_i} - \vec{u}}{c_s^2} + \frac{\vec{c_i}(\vec{c_i}\cdot \vec{u})}{c_s^4}] \cdot \vec{F} $$

The macroscopic velocity also needs to be shifted in accordance:

$$
\vec{u} =  \frac{1}{\rho} \sum_{i=0}^8\vec{c_i} f_i + \frac{\vec{F}\triangle t}{2 \rho}
$$

## Fluid ageing

Fluid age is approximated through the use of a passive scalar field.
The age distribution function is represented as:

$$
fa_i(x+\vec{c_i} \delta t, t+\delta t) = (A+\delta t)\cdot f_i(x,t),
$$

where $ A $ is the macroscopic age of a fluid cell.

The macroscopic age can be calculated as:
$$
A = \frac{1}{\rho} \sum_{i=0}^8 fa_i.
$$

The inflow boundaries are considered as $ A = 0 $.
The outflow boundaries implemented using the hypothesis that the unknown distributions can not change the age of a cell. Thus, the macroscopic age is calculated from the known distributions.

The solid boundaries are simple bounce-back surfaces for fluid ages.
