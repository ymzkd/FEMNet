"""
03_modal_analysis.py - Modal Analysis (Eigenvalue Problem)

This example demonstrates:
- Creating a cantilever beam whose distributed mass is lumped to the nodes
- Running modal (eigenvalue) analysis
- Extracting natural frequencies and mode shapes
- Comparing with the exact solution of a uniform cantilever beam

Structure:
      m/2n    m/n   m/n   m/n   m/2n   (lumped mass)
    ▓==o======o=====o=====o=====o
    ▓   L = 3000 mm, n elements
  Fixed

Units: N, mm, s. Mass is t (= N·s^2/mm).
FEMNet takes node mass as weight [N] (MassData.Mass) and divides it by
the gravitational acceleration (FEModel.GraityAccel [mm/s^2]).

"""
import sys
sys.path.insert(0, '..')

from femnet import *

# 支持条件の指定(ConstraintType: Free / Fix。ばね支持は SupportSpringElement を要素として追加)
FREE = ConstraintType_Free
FIX = ConstraintType_Fix
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

# For 2D beam in XY plane, keep only the in-plane bending DOFs free
# Note: Rz has no mass, but the vibration analysis condenses massless DOFs
#       statically, so it does not need to be constrained
for i in range(1, model.NodeNum()):
    node = model.GetNode(i)
    # Free: Uy, Rz (in-plane bending)
    # Fixed: Ux, Uz, Rx, Ry
    node.Fix = Support(FIX, FREE, FIX, FIX, FIX, FREE)

print(f"   Total DOF: {model.DOFNum()}")
print(f"   Free DOF: {model.FreeDOFNum()}")

# ============================================================
# 3. Define material and section
# ============================================================
print("\n3. Adding material and section...")

# Steel material with unit weight
# (element self-weight is not used as mass here: ComputeElementNodeMass() is not called)
E = 205e3        # Young's modulus [N/mm^2]
nu = 0.3         # Poisson's ratio
gamma = 7.85e-5  # Unit weight [N/mm^3] (78.5 kN/m^3)
model.AddMaterialWithDensity(E, nu, gamma)

# Rectangular section: 100mm x 200mm
b = 100.0  # width [mm]
h = 200.0  # height [mm]
A = b * h  # Cross-sectional area
Iy = b * h**3 / 12  # Moment of inertia (strong axis)
Iz = h * b**3 / 12  # Moment of inertia (weak axis - for bending in XY plane)
K = b * h**3 / 3 * (1 - 0.63 * b / h)
model.AddSection(A, Iy, Iz, K)

print(f"   Material: Steel (E={E/1e3:.0f} GPa, gamma={gamma*1e6:.1f} kN/m^3)")
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
# 5. Lump the distributed mass to the nodes
# ============================================================
print("\n5. Adding lumped masses...")

# Uniformly distributed mass of 1000 kg (= 1.0 t) over the beam.
# Each element puts half of its mass on each end node, so the end nodes
# get m/2n and the interior nodes get m/n.
total_mass = 1.0                     # [t] (1000 kg)
m_elem = total_mass / n_elements     # [t] per element
g = model.GraityAccel                # [mm/s^2]

for i in range(model.NodeNum()):
    node_mass = m_elem if 0 < i < n_elements else 0.5 * m_elem
    # FEMNet takes the node mass as weight [N]
    model.GetNode(i).MassData.Mass = node_mass * g

print(f"   Mass per interior node: {m_elem * 1e3:.0f} kg (end nodes: {0.5 * m_elem * 1e3:.0f} kg)")
print(f"   Total model mass: {model.SumNodeMass() / g * 1e3:.0f} kg "
      f"(weight {model.SumNodeMass() / 1e3:.2f} kN)")

# ============================================================
# 6. Run modal analysis
# ============================================================
print("\n6. Running modal analysis...")

# Number of modes to extract
n_modes = 3

# Solve eigenvalue problem (operator owns the computation)
vib = FEVibrationAnalysis(model)
result = vib.Compute(n_modes)
eigen_values = vib.EigenValues()
mode_vectors = vib.ModeVectors()

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

# Exact solution of a uniform cantilever beam with distributed mass:
#   omega_n = (beta_n*L)^2 * sqrt(E*I / (m_bar * L^4)),  m_bar = total_mass / L [t/mm]
#   beta_n*L = 1.8751, 4.6941, 7.8548 for modes 1-3
print("\n   === Theoretical Comparison (uniform cantilever beam) ===")
I = Iz  # For bending in XY plane
m_bar = total_mass / L
beta_L = [1.8751, 4.6941, 7.8548]

print(f"   {'Mode':<6} {'Theory [rad/s]':<16} {'FEM [rad/s]':<14} {'Error':<8}")
print("   " + "-" * 44)
for i in range(min(len(beta_L), len(eigen_values))):
    omega_theory = beta_L[i] ** 2 * math.sqrt(E * I / (m_bar * L ** 4))
    error = (eigen_values[i] - omega_theory) / omega_theory * 100
    print(f"   {i+1:<6} {omega_theory:<16.4f} {eigen_values[i]:<14.4f} {error:+.2f}%")
print(f"   (Lumping the mass to {n_elements} elements lowers the accuracy of higher modes)")

print("\n" + "=" * 60)
print("Example 03 completed successfully!")
print("=" * 60)
