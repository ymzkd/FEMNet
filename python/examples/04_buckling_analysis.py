"""
04_buckling_analysis.py - Linear Buckling Analysis

This example demonstrates:
- Creating an Euler column (simply supported at both ends)
- Applying axial compression load
- Running linear static analysis first
- Running linear buckling analysis
- Extracting critical load factor and buckling mode

Structure:
        P (compression)
        ↓
        o  <- Pinned (Uy, Uz, Rx, Ry, Rz fixed)
        |
        |  L = 3000 mm
        |
        o  <- Pinned (Ux, Uy, Uz, Rx, Ry, Rz fixed)

Euler column with pinned-pinned ends

"""
import sys
sys.path.insert(0, '..')

from femnet import *
import math

print("=" * 60)
print("FEMNet Example 04: Linear Buckling Analysis")
print("=" * 60)

# ============================================================
# 1. Create model and define geometry
# ============================================================
print("\n1. Creating Euler column model...")

model = FEModel()

# Geometry parameters
L = 3000.0  # Column length [mm]
n_elements = 10  # Number of beam elements

# Create nodes along the column (vertical)
print("   Adding nodes...")
for i in range(n_elements + 1):
    y = i * L / n_elements
    model.AddNode(i, 0, y, 0)

print(f"   Number of nodes: {model.NodeNum()}")

# ============================================================
# 2. Define supports (boundary conditions)
# ============================================================
print("\n2. Setting boundary conditions...")

# For pinned-pinned column buckling in XY plane:
# - Bottom node: Pin support (Ux, Uy, Uz fixed, rotations free for bending)
# - Top node: Roller (Ux, Uz fixed, Uy free for axial shortening, rotations free)
# - Intermediate: All translational DOFs free for buckling, out-of-plane fixed

# Bottom node: Pinned (translations fixed, Rz free for bending)
model.GetNode(0).Fix = Support(True, True, True, True, True, False)

# Top node: Roller in Y direction (Ux, Uz fixed; Uy, Rz free)
top_node = model.GetNode(n_elements)
top_node.Fix = Support(True, False, True, True, True, False)

# Intermediate nodes: Free in Ux, Uy, Rz for buckling mode
# Fixed: Uz, Rx, Ry (out-of-plane constraints)
for i in range(1, n_elements):
    node = model.GetNode(i)
    node.Fix = Support(False, False, True, True, True, False)
    node.Fix.UnlockAllRot()

print(f"   Total DOF: {model.DOFNum()}")
print(f"   Free DOF: {model.FreeDOFNum()}")

# ============================================================
# 3. Define material and section
# ============================================================
print("\n3. Adding material and section...")

# Steel material
E = 205e3  # Young's modulus [N/mm^2]
nu = 0.3
model.AddMaterial(E, nu)

# Rectangular section: 50mm x 100mm
b = 50.0   # width [mm] (weak direction)
h = 100.0  # height [mm] (strong direction)
A = b * h  # Cross-sectional area
Iy = b * h**3 / 12  # Strong axis (about y)
Iz = h * b**3 / 12  # Weak axis (about z) - buckling occurs about this axis
K = b * h**3 / 3 * (1 - 0.63 * b / h)
model.AddSection(A, Iy, Iz, K)

print(f"   Material: Steel (E={E/1e3:.0f} GPa)")
print(f"   Section: Rectangular {b:.0f}x{h:.0f}mm")
print(f"     A = {A:.0f} mm^2")
print(f"     Iz = {Iz:.2e} mm^4 (weak axis)")

# ============================================================
# 4. Create beam elements
# ============================================================
print("\n4. Creating beam elements...")

for i in range(n_elements):
    model.add_beam_element(i, i, i + 1, 0, 0, 0.0)

print(f"   Number of elements: {len(model.Elements)}")

# ============================================================
# 5. Apply reference axial load
# ============================================================
print("\n5. Applying reference compression load...")

# Reference load (unit load for buckling factor calculation)
P_ref = -1000.0  # 1 kN compression [N]

loads = VectorLoad()
node_load = NodeLoad(n_elements, 0, P_ref, 0)  # Apply at top node
loads.append(node_load)

print(f"   Reference load: {abs(P_ref)/1000:.1f} kN compression at top")

# ============================================================
# 6. Run linear static analysis first
# ============================================================
print("\n6. Running linear static analysis...")

static_solver = FELinearStaticOp(model, loads)
static_solver.Compute()

if static_solver.Computed():
    print("   Static analysis completed successfully!")
else:
    print("   Static analysis failed!")
    sys.exit(1)

# ============================================================
# 7. Run buckling analysis
# ============================================================
print("\n7. Running linear buckling analysis...")

buckling_solver = FEBucklingAnalysis(static_solver)
buckling_solver.mode_num = 3  # Extract first 3 buckling modes
result = buckling_solver.SolveBuckling()

if result > 0:
    print(f"   Buckling analysis completed! {result} modes extracted.")
else:
    print(f"   Buckling analysis failed with error code: {result}")
    sys.exit(1)

# ============================================================
# 8. Display results
# ============================================================
print("\n8. Results:")

# Critical load factors
print("\n   === Critical Load Factors ===")
print(f"   {'Mode':<6} {'Load Factor':<15} {'Critical Load [kN]':<20}")
print("   " + "-" * 45)

eigen_values = buckling_solver.EigenValues()
for i in range(len(eigen_values)):
    load_factor = eigen_values[i]
    P_cr = abs(P_ref) * load_factor / 1000  # Convert to kN
    print(f"   {i+1:<6} {load_factor:<15.4f} {P_cr:<20.4f}")

# First buckling mode shape
print("\n   === 1st Buckling Mode Shape ===")
mode_vectors = buckling_solver.ModeVectors()
if len(mode_vectors) > 0:
    print(f"   {'Node':<6} {'y [mm]':<10} {'Mode 1 Ux':<12}")
    print("   " + "-" * 30)

    for node_id in range(model.NodeNum()):
        node = model.GetNode(node_id)
        mode1_dx = mode_vectors[0][node_id].Dx()
        print(f"   {node_id:<6} {node.Location.y:<10.0f} {mode1_dx:<12.6f}")

# Theoretical comparison (Euler's formula for pinned-pinned column)
# P_cr = pi^2 * E * I / L^2
print("\n   === Theoretical Comparison (Euler's Formula) ===")
I = Iz  # Buckling about weak axis
P_cr_euler = math.pi**2 * E * I / L**2

print(f"   Euler critical load (1st mode):")
print(f"     P_cr = pi^2*EI/L^2 = {P_cr_euler/1000:.4f} kN")
print(f"")
print(f"   FEM critical load (1st mode):")
P_cr_fem = abs(P_ref) * eigen_values[0] / 1000
print(f"     P_cr = {P_cr_fem:.4f} kN")
print(f"")
print(f"   Error: {abs(P_cr_fem - P_cr_euler/1000) / (P_cr_euler/1000) * 100:.2f}%")

# Higher modes comparison
print("\n   === Higher Modes (Theoretical) ===")
for n in range(1, min(4, len(eigen_values) + 1)):
    P_cr_n = n**2 * P_cr_euler / 1000
    P_cr_fem_n = abs(P_ref) * eigen_values[n-1] / 1000 if n <= len(eigen_values) else 0
    print(f"   Mode {n}: P_cr = {n}^2 x {P_cr_euler/1000:.2f} = {P_cr_n:.2f} kN (FEM: {P_cr_fem_n:.2f} kN)")

print("\n" + "=" * 60)
print("Example 04 completed successfully!")
print("=" * 60)
