"""
03_modal_analysis.py - Modal Analysis (Eigenvalue Problem)

This example demonstrates:
- Creating a simple cantilever beam with lumped mass
- Running modal (eigenvalue) analysis
- Extracting natural frequencies and mode shapes

Structure:
                      m (lumped mass)
    ▓==================|
    ▓                  ↓
    ▓   L = 3000 mm
  Fixed

Cantilever beam with tip mass

"""
import sys
sys.path.insert(0, '..')

from femnet import *
import math

print("=" * 60)
print("FEMNet Example 03: Modal Analysis (Eigenvalue Problem)")
print("=" * 60)

# ============================================================
# 1. Create model and define geometry
# ============================================================
print("\n1. Creating cantilever beam model...")

model = FEModel()

# Geometry parameters
L = 3000.0  # Beam length [mm]
n_elements = 5  # Number of beam elements

# Create nodes along the beam
print("   Adding nodes...")
for i in range(n_elements + 1):
    x = i * L / n_elements
    model.AddNode(i, x, 0, 0)

print(f"   Number of nodes: {model.NodeNum()}")

# ============================================================
# 2. Define supports (boundary conditions)
# ============================================================
print("\n2. Setting boundary conditions...")

# Fixed support at node 0 (all DOFs constrained)
model.GetNode(0).Fix.FixAll()

# For 2D beam in XY plane, constrain out-of-plane for all other nodes
# Only keep Uy free (translational DOF for bending)
# Note: Rz is constrained because mass matrix only includes translational mass
for i in range(1, model.NodeNum()):
    node = model.GetNode(i)
    # Free: Uy only (bending displacement)
    # Fixed: All others (Ux, Uz, Rx, Ry, Rz)
    node.Fix = Support(True, False, True, True, True, True)

print(f"   Total DOF: {model.DOFNum()}")
print(f"   Free DOF: {model.FreeDOFNum()}")

# ============================================================
# 3. Define material and section
# ============================================================
print("\n3. Adding material and section...")

# Steel material with density
E = 205e3       # Young's modulus [N/mm^2]
nu = 0.3        # Poisson's ratio
rho = 7.85e-9   # Density [kg/mm^3] (7850 kg/m^3)
model.AddMaterialWithDensity(E, nu, rho)

# Rectangular section: 100mm x 200mm
b = 100.0  # width [mm]
h = 200.0  # height [mm]
A = b * h  # Cross-sectional area
Iy = b * h**3 / 12  # Moment of inertia (strong axis)
Iz = h * b**3 / 12  # Moment of inertia (weak axis - for bending in XY plane)
K = b * h**3 / 3 * (1 - 0.63 * b / h)
model.AddSection(A, Iy, Iz, K)

print(f"   Material: Steel (E={E/1e3:.0f} GPa, rho={rho*1e9:.0f} kg/m^3)")
print(f"   Section: Rectangular {b:.0f}x{h:.0f}mm")
print(f"     A = {A:.0f} mm^2")
print(f"     Iz = {Iz:.2e} mm^4 (for bending in XY plane)")

# ============================================================
# 4. Create beam elements
# ============================================================
print("\n4. Creating beam elements...")

for i in range(n_elements):
    model.add_beam_element(i, i, i + 1, 0, 0, 0.0)

print(f"   Number of elements: {len(model.Elements)}")

# ============================================================
# 5. Add lumped mass to all free nodes
# ============================================================
print("\n5. Adding lumped masses...")

# Distribute mass to all free nodes
total_mass = 1000.0  # 1000 kg total
mass_per_node = total_mass / (n_elements)  # Distribute among free nodes

for i in range(1, model.NodeNum()):
    node = model.GetNode(i)
    node.MassData.Mass = mass_per_node

print(f"   Mass per node: {mass_per_node:.0f} kg")
print(f"   Total model mass: {model.SumNodeMass():.0f} kg")

# ============================================================
# 6. Run modal analysis
# ============================================================
print("\n6. Running modal analysis...")

# Number of modes to extract
n_modes = 3

# Create output containers
eigen_values = VectorDouble()
mode_vectors = VectorMode()

# Solve eigenvalue problem
result = model.SolveVibration(n_modes, eigen_values, mode_vectors)

if result > 0:
    print(f"   Modal analysis completed! {result} modes extracted.")
elif result == 0:
    print("   Modal analysis returned 0 modes.")
    sys.exit(1)
else:
    print(f"   Modal analysis failed with error code: {result}")
    sys.exit(1)

# ============================================================
# 7. Display results
# ============================================================
print("\n7. Results:")

# Natural frequencies and periods
print("\n   === Natural Frequencies and Periods ===")
print(f"   {'Mode':<6} {'omega [rad/s]':<15} {'f [Hz]':<12} {'T [s]':<12}")
print("   " + "-" * 50)

for i in range(len(eigen_values)):
    omega = eigen_values[i]
    f = omega / (2 * math.pi)  # Frequency in Hz
    T = 1 / f if f > 0 else float('inf')  # Period in seconds
    print(f"   {i+1:<6} {omega:<15.4f} {f:<12.4f} {T:<12.4f}")

# Mode shapes
print("\n   === Mode Shapes (Y-displacement along beam) ===")
if len(mode_vectors) > 0:
    print(f"   {'Node':<6} {'x [mm]':<10} {'Mode 1':<12} {'Mode 2':<12} {'Mode 3':<12}")
    print("   " + "-" * 56)

    for node_id in range(model.NodeNum()):
        node = model.GetNode(node_id)
        mode1_dy = mode_vectors[0][node_id].Dy() if len(mode_vectors) > 0 else 0
        mode2_dy = mode_vectors[1][node_id].Dy() if len(mode_vectors) > 1 else 0
        mode3_dy = mode_vectors[2][node_id].Dy() if len(mode_vectors) > 2 else 0
        print(f"   {node_id:<6} {node.Location.x:<10.0f} {mode1_dy:<12.6f} {mode2_dy:<12.6f} {mode3_dy:<12.6f}")

# Theoretical estimate for first mode of cantilever with distributed mass
# For uniform beam: omega_n = (beta_n)^2 * sqrt(E*I / (rho*A*L^4))
# beta_1 = 1.875 for first mode
print("\n   === Theoretical Comparison ===")
I = Iz  # For bending in XY plane

# Approximate formula for cantilever with lumped masses (SDOF approximation)
# Using equivalent mass at tip: m_eff = 0.23 * m_beam + m_tip
# For our case: total_mass distributed uniformly, treat as m_eff ~ 0.5 * total_mass
m_eff = 0.5 * total_mass
omega_theory = math.sqrt(3 * E * I / (m_eff * L**3))
f_theory = omega_theory / (2 * math.pi)
T_theory = 1 / f_theory

print(f"   1st mode (approximate SDOF model):")
print(f"     Equivalent tip mass: {m_eff:.0f} kg")
print(f"     Theoretical omega: {omega_theory:.4f} rad/s")
print(f"     Theoretical f:     {f_theory:.4f} Hz")
print(f"     Theoretical T:     {T_theory:.4f} s")
print(f"")
print(f"     FEM omega: {eigen_values[0]:.4f} rad/s")
print(f"     FEM f:     {eigen_values[0]/(2*math.pi):.4f} Hz")
print(f"     FEM T:     {2*math.pi/eigen_values[0]:.4f} s")
print("   (Note: Difference expected due to distributed mass approximation)")

print("\n" + "=" * 60)
print("Example 03 completed successfully!")
print("=" * 60)
