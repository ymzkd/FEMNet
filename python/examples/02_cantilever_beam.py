"""
02_cantilever_beam.py - Cantilever Beam Analysis

This example demonstrates:
- Creating a cantilever beam with BeamElement
- Applying point loads and distributed loads
- Running linear static analysis
- Extracting displacements, reactions, and beam stresses

Structure:
                w (distributed load)
    ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓ ↓
    ============================= → P (axial)
    ▓                           |
    ▓   L = 3000 mm             ↓ F (point load)
    ▓
  Fixed

"""
import sys
sys.path.insert(0, '..')

from femnet import *
import math

print("=" * 60)
print("FEMNet Example 02: Cantilever Beam Analysis")
print("=" * 60)

# ============================================================
# 1. Create model and define geometry
# ============================================================
print("\n1. Creating cantilever beam model...")

model = FEModel()

# Geometry parameters
L = 3000.0  # Beam length [mm]
n_elements = 6  # Number of beam elements

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
for i in range(1, model.NodeNum()):
    node = model.GetNode(i)
    # Free: Ux, Uy, Rz (in-plane DOFs)
    # Fixed: Uz, Rx, Ry (out-of-plane DOFs)
    node.Fix = Support(False, False, True, True, True, False)

print(f"   Total DOF: {model.DOFNum()}")
print(f"   Free DOF: {model.FreeDOFNum()}")

# ============================================================
# 3. Define material and section
# ============================================================
print("\n3. Adding material and section...")

# Steel material: E = 205 GPa, nu = 0.3
E = 205e3  # Young's modulus [N/mm^2]
nu = 0.3
model.AddMaterial(E, nu)

# Rectangular section: 100mm x 200mm (width x height)
b = 100.0  # width [mm]
h = 200.0  # height [mm]
A = b * h  # Cross-sectional area
Iy = b * h**3 / 12  # Moment of inertia about y-axis (strong axis)
Iz = h * b**3 / 12  # Moment of inertia about z-axis (weak axis)
K = b * h**3 / 3 * (1 - 0.63 * b / h)  # Torsion constant (approximate)
model.AddSection(A, Iy, Iz, K)

print(f"   Material: Steel (E={E/1e3:.0f} GPa)")
print(f"   Section: Rectangular {b:.0f}x{h:.0f}mm")
print(f"     A = {A:.0f} mm^2")
print(f"     Iy = {Iy:.2e} mm^4")
print(f"     Iz = {Iz:.2e} mm^4")

# ============================================================
# 4. Create beam elements
# ============================================================
print("\n4. Creating beam elements...")

# add_beam_element(id, n1_id, n2_id, sec_id, mat_id, beta)
# beta: rotation angle about the beam axis (0 for standard orientation)
for i in range(n_elements):
    model.add_beam_element(i, i, i + 1, 0, 0, 0.0)

print(f"   Number of elements: {len(model.Elements)}")

# ============================================================
# 5. Apply loads
# ============================================================
print("\n5. Applying loads...")

# Create load vector
loads = VectorLoad()

# Point load at free end (node n_elements)
F = -10000.0  # 10 kN downward [N]
tip_node_id = n_elements
point_load = NodeLoad(tip_node_id, 0, F, 0)
loads.append(point_load)
print(f"   Point load: {-F/1000:.1f} kN downward at tip (node {tip_node_id})")

# Distributed load on all beam elements
w = -5.0  # 5 N/mm downward (uniform)
print(f"   Distributed load: {-w:.1f} N/mm on all elements")

for i in range(n_elements):
    beam_elem = model.GetBeamElement(i)
    # Create uniform distributed load using BeamPolyLoad
    # Parameters: [w_start, w_end], [param_start, param_end], element, axis
    dist_load = BeamPolyLoad(VectorDouble([w, w]), VectorDouble([0.0, 1.0]), beam_elem, YAxis)
    loads.append(dist_load)

# ============================================================
# 6. Run linear static analysis
# ============================================================
print("\n6. Running linear static analysis...")

solver = FELinearStaticOp(model, loads)
solver.Compute()

if solver.Computed():
    print("   Analysis completed successfully!")
else:
    print("   Analysis failed!")
    sys.exit(1)

# ============================================================
# 7. Display results
# ============================================================
print("\n7. Results:")

# Nodal displacements
print("\n   === Nodal Displacements ===")
displacements = solver.GetDisplacements()
print(f"   {'Node':<6} {'x [mm]':<10} {'Dy [mm]':<14} {'Rz [rad]':<14}")
print("   " + "-" * 50)

for i, disp in enumerate(displacements):
    node = model.GetNode(i)
    print(f"   {i:<6} {node.Location.x:<10.0f} {disp.Dy():<14.4f} {disp.Rz():<14.6f}")

# Maximum displacement at tip
tip_disp = displacements[n_elements]
print(f"\n   Maximum tip deflection: {abs(tip_disp.Dy()):.4f} mm")

# Reaction forces at fixed support
print("\n   === Reaction Forces at Fixed Support (Node 0) ===")
reactions = solver.GetReactForces()
for react in reactions:
    if react.id == 0:
        print(f"   Px = {react.Px():>12.2f} N")
        print(f"   Py = {react.Py():>12.2f} N")
        print(f"   Mz = {react.Mz():>12.2f} N-mm")

# Theoretical comparison for point load only:
# delta = P*L^3 / (3*E*I)
# theta = P*L^2 / (2*E*I)
# Note: For Y-direction bending, we use Iz (bending about Z-axis)
I = Iz  # Use Iz for bending in XY plane (Y-direction loading)
delta_point_theory = abs(F) * L**3 / (3 * E * I)
theta_point_theory = abs(F) * L**2 / (2 * E * I)

# For distributed load:
# delta = w*L^4 / (8*E*I)
# theta = w*L^3 / (6*E*I)
delta_dist_theory = abs(w) * L**4 / (8 * E * I)
theta_dist_theory = abs(w) * L**3 / (6 * E * I)

total_delta_theory = delta_point_theory + delta_dist_theory
total_theta_theory = theta_point_theory + theta_dist_theory

print("\n   === Theoretical Comparison ===")
print(f"   Tip deflection (theory): {total_delta_theory:.4f} mm")
print(f"   Tip deflection (FEM):    {abs(tip_disp.Dy()):.4f} mm")
print(f"   Error: {abs(abs(tip_disp.Dy()) - total_delta_theory) / total_delta_theory * 100:.2f}%")

print(f"\n   Tip rotation (theory): {total_theta_theory:.6f} rad")
print(f"   Tip rotation (FEM):    {abs(tip_disp.Rz()):.6f} rad")
print(f"   Error: {abs(abs(tip_disp.Rz()) - total_theta_theory) / total_theta_theory * 100:.2f}%")

# Equilibrium check
print("\n   === Equilibrium Check ===")
total_dist_load = abs(w) * L  # Total distributed load
total_applied_Py = F + w * L
print(f"   Total applied Py: {total_applied_Py:.2f} N")

reaction_Py = sum(r.Py() for r in reactions if r.id == 0)
print(f"   Sum Py (reactions + applied): {reaction_Py + total_applied_Py:.6f} N (should be 0)")

print("\n" + "=" * 60)
print("Example 02 completed successfully!")
print("=" * 60)
