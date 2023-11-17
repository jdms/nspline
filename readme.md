# Building curves with the NSpline class

Let $\gamma : I\subset \mathbb{R} \rightarrow \mathbb{R}^3$ be a smooth curve
for which we know a set of N samples: 

$$\\{(t_1, \gamma(t_1)), (t_2,\gamma(t_2)), \ldots, (t_N, \gamma(t_N))\\} \subset I \times \mathbb{R}^3.$$

We are interested in interpolating $\gamma$, *i.e.* in finding a new curve
$\alpha$ such that 

$$\alpha(t_i) = \gamma(t_i), ~ \text{for all} ~ 1 \leq i \leq N.$$
{#eq:points_and_centers}

Alternatively, we may only know a sampling of the curve's trace, *i.e.*:

$$\\{\gamma_1, \gamma_2, \ldots, \gamma_N\\} \subset \mathbb{R}^3,$$

and, in this case, we could enforce the condition:

$$\alpha(i) = \gamma_i, ~ \text{for all} ~ 1 \leq i \leq N.$$
{#eq:only_points} 

For concreteness, lets tackle the second case.  We need to
interpolate the three coordinates of $\alpha$, namely we need to find the functions: 

$$\alpha(t) = ( x(t), y(t), z(t) ) : [1, N] \rightarrow \mathbb{R}^3,$$
{#eq:components}

subject to the conditions given in @eq:only_points.  This is accomplished with the following code: 

```{.cpp}
#include <numeric> // for std::iota

// (...)

int N; // number of samples (given)

// First component of points sampled from gamma (given)
std::vector<double> gamma_x = { gamma1_x, gamma2_x, ..., gammaN_x};

// Second component of points sampled from gamma (given)
std::vector<double> gamma_y = { gamma1_y, gamma2_y, ..., gammaN_y};

// Third component of points sampled from gamma (given)
std::vector<double> gamma_z = { gamma1_z, gamma2_z, ..., gammaN_z};

std::vector<double> t(N);
std::iota(t.begin(), t.end(), 1); // makes t = { 1, 2, ..., N }

// Create the actual interpolants

jdms::NSpline x(t, gamma_x), y(t, gamma_y), z(t, gamma_z);
```

The interpolants are evaluated with `NSpline::operator()`, thus `x(1) ==
gamma1_x` and so on.  Derivatives and second derivatives can be computed with
methods `NSpline::D()` and `NSpline::D2()`. Now, a few useful formulas.

### Tangent vector to $\alpha$ at a point $t$, $T(t)$: 

$$T(t) = \displaystyle \frac{\alpha^\prime(t)}{|\alpha^\prime(t)|},$$
{#eq:tangents}

provided that $|\alpha^\prime(t)| \neq 0$.


### Curvature of $\alpha$ at a point $t$, $k(t)$:

$$k(t) = \displaystyle \frac{|\alpha^\prime(t) \wedge \alpha^{\prime\prime}(t)|}{|\alpha^\prime(t)|^3},$$
{#eq:curvatures}

where $\wedge$ is the vector product.


### Normal vector of $\alpha$ at a point $t$, $N(t)$:

$$N(t) = \frac{T^\prime(t)}{k(t)}; ~ \text{where} ~
    T^\prime(t) = \displaystyle \frac{\alpha^{\prime\prime}(t)}{|\alpha^\prime(t)|}
    - \frac{T(t)}{|\alpha^\prime(t)|} \langle\alpha^{\prime\prime}(t), T(t)\rangle,$$
{#eq:normals}

provided that $k(t) \neq 0$ and $|\alpha^\prime(t)| \neq 0$, see @eq:curvatures and @eq:tangents.  Here $\langle, \rangle$ stands for the inner product
$\langle x, y \rangle = \sum_{k} x_k y_k.$

