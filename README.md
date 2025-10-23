# Conjugate Heat Transfer (CHT) — Plate with Hollow Border and Channel Flow

![Geometry overview](docs/cht-geometry.png)

> **TL;DR**  
> 2D **conjugate heat transfer** (CHT) setup: an outer channel with two inlet/outlet “tabs” on each vertical wall and a centered **solid plate** (inner square).  
> This repo includes:
> - A clean **mathematical model** (strong & weak forms).
> - **FreeFEM** scripts to generate **matching meshes** for fluid and solid (shared interface nodes).
> - A LaTeX project you can extend into a report.

---

## 1) Context & Motivation

Conjugate heat transfer couples:
- **Fluid mechanics** (incompressible Navier–Stokes) in the **fluid domain** $\Omega_f$, and
- **Heat conduction** in the **solid** $\Omega_s$,

with **thermal exchange** across the common interface $\Gamma_{fs}$. This configuration models a plate cooled by channel flow: the plate acts as an internal boundary (a hole from the fluid’s point of view, and the solid itself for conduction).

---

## 2) Geometry & Domains

- Outer rectangle (channel): $[x_L, x_R] \times [y_B, y_T] = [0,1.5]\times[-1,0]$, with **two rectangular tabs** on each vertical side (slots you can later label as inlets/outlets).
- Inner **solid plate**: centered square $[x_s, x_e]\times[y_s, y_e]=[0.10,1.40]\times[-0.90,-0.10]$.

Notation:
- $\Omega_f$: fluid domain (outer region **minus** inner square).
- $\Omega_s$: solid domain (the inner square).
- $\Gamma_{fs} = \overline{\Omega_f}\cap\overline{\Omega_s}$: fluid–solid interface (square boundary).
- $\Gamma_w$: remaining channel walls (outer boundary).
- Optionally: $\Gamma_{\text{in}}^{(1,2)}$, $\Gamma_{\text{out}}^{(1,2)}$ on tab segments.

**Boundary labels in meshes**
- `30` → $\Gamma_w$ (outer walls)
- `40` → $\Gamma_{fs}$ (fluid–solid interface)

**Orientation (important for BAMG)**
- Outer boundary: **CCW**.
- Inner boundary for the fluid (hole): **CW**.
- Inner boundary for the solid (filled): **CCW**.
- Same discretization count on $\Gamma_{fs}$ for fluid and solid ⇒ **shared interface nodes**.

---

## 3) Governing Equations (Strong Form)

Unknowns (in $d=2$):
- Fluid: velocity $\mathbf{u}(x,t)$, pressure $p(x,t)$, and temperature $T_f(x,t)$ in $\Omega_f$.
- Solid: temperature $T_s(x,t)$ in $\Omega_s$.

Parameters:
- $\rho$: density; $\mu$: dynamic viscosity.  
- $\kappa$: thermal diffusivity (solid); $\hat\kappa$: thermal diffusivity (fluid).  
- $\alpha$: interface heat transfer coefficient (Robin coupling).

**Solid (heat conduction) in $\Omega_s$**
$$
\partial_t T_s - \kappa \Delta T_s = 0 \quad \text{in } \Omega_s\times(0,T).
$$

**Fluid (incompressible Navier–Stokes) in $\Omega_f$**

$$
\begin{aligned}
\rho(\partial_t \mathbf{u} + \mathbf{u}\cdot\nabla \mathbf{u})
&= -\nabla p + \mu \Delta \mathbf{u} + T_f\,\mathbf{f}_{ext}
\quad \text{in } \Omega_f\times(0,T),\\
\nabla\cdot \mathbf{u} &= 0 \quad \text{in } \Omega_f\times(0,T).
\end{aligned}
$$

with no-slip $\mathbf{u}=\mathbf{0}$ on $\Gamma_w\cup\Gamma_{fs}$.  
On inlet/outlet windows you can use **do-nothing** tractions or prescribe profiles.

**Fluid temperature (advection–diffusion) in $\Omega_f$**
$$
\partial_t T_f - \hat\kappa \Delta T_f + \mathbf{u}\cdot\nabla T_f = 0
\quad \text{in } \Omega_f\times(0,T).
$$

**Coupling and BCs**
- Interface (Robin–Robin) on $\Gamma_{fs}$:
$$
\partial_{\mathbf{n}}T_s = \alpha (T_f - T_s), 
\qquad
\partial_{\mathbf{n}}T_f = \alpha (T_s - T_f).
$$
- Walls $\Gamma_w$: $\mathbf{u}=\mathbf{0}$, $\partial_{\mathbf{n}}T_f=0$.
- Inlet/Outlet tabs: optional temperature Dirichlet $T_f=T_{\text{in/out}}$ or natural.
- Initial data: $ \mathbf{u}(x,0)=\mathbf{u}_0(x)$, $T_f(x,0)=T_{f,0}(x)$, $T_s(x,0)=T_{s,0}(x)$.

---

## 4) Weak Formulation (Variational)

Spaces:
$$
\begin{aligned}
V &:= \{\mathbf{v}\in H^1(\Omega_f)^d:\ \mathbf{v}=\mathbf{0}\ \text{on }\Gamma_w\cup\Gamma_{fs}\},\\
Q &:= L^2_0(\Omega_f) = \{q\in L^2(\Omega_f): \int_{\Omega_f} q = 0\},\\
W_s &:= H^1(\Omega_s),\qquad
W_f := \{\phi\in H^1(\Omega_f): \phi=0 \text{ on Dirichlet parts}\}.
\end{aligned}
$$

Con $D\mathbf{u}=\tfrac12(\nabla\mathbf{u}+\nabla\mathbf{u}^\top)$ y tests $\mathbf{v}\in V$, $q\in Q$, $\phi_s\in W_s$, $\phi_f\in W_f$:

**Solid heat**
$$
\int_{\Omega_s}\partial_t T_s\,\phi_s
+\kappa\int_{\Omega_s}\nabla T_s\cdot\nabla\phi_s
+\kappa\alpha\int_{\Gamma_{fs}} T_s\,\phi_s
= \kappa\alpha\int_{\Gamma_{fs}} T_f\,\phi_s.
$$

**Fluid temperature**
$$
\int_{\Omega_f}\partial_t T_f\,\phi_f
+\int_{\Omega_f}(\mathbf{u}\cdot\nabla T_f)\,\phi_f
+\hat\kappa\int_{\Omega_f}\nabla T_f\cdot\nabla\phi_f
+\hat\kappa\alpha\int_{\Gamma_{fs}} T_f\,\phi_f
= \hat\kappa\alpha\int_{\Gamma_{fs}} T_s\,\phi_f.
$$

**Navier–Stokes**
$$
\begin{aligned}
&\rho\int_{\Omega_f}\partial_t\mathbf{u}\cdot\mathbf{v}
+\rho\int_{\Omega_f}(\mathbf{u}\cdot\nabla)\mathbf{u}\cdot\mathbf{v}
+2\mu\int_{\Omega_f} D\mathbf{u}:D\mathbf{v}
-\int_{\Omega_f} p\,\nabla\cdot\mathbf{v}
= \int_{\Omega_f} T_f\,\mathbf{f}_{ext}\cdot\mathbf{v},\\
&\int_{\Omega_f} q\,\nabla\cdot\mathbf{u}=0.
\end{aligned}
$$

This matches the LaTeX derivation included in the project.

