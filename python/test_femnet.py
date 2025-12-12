"""
FEMNet Python bindings test script - Full solver example
"""
import sys
sys.path.insert(0, '.')

from femnet import *

print("=" * 60)
print("FEMNet Python Bindings - Linear Static Analysis Test")
print("=" * 60)

# ============================================================
# 1. Create a cantilever beam model
# ============================================================
print("\n1. Creating cantilever beam model...")

model = FEModel()

# Beam parameters
L = 3000.0       # Total length [mm]
n_elem = 3       # Number of elements
elem_L = L / n_elem

# Add nodes
print("   Adding nodes...")
for i in range(n_elem + 1):
    model.AddNode(i, i * elem_L, 0, 0)
print(f"   Number of nodes: {model.NodeNum()}")

# Fix first node (cantilever support)
print("   Setting boundary conditions (fixed at node 0)...")
model.GetNode(0).Fix.FixAll()

# Add material: Steel
print("   Adding material (Steel: E=205GPa, nu=0.3)...")
model.AddMaterial(205e3, 0.3)

# Add section: Rectangular 100x10 mm
print("   Adding section (A=1000mm2, Iy=8333mm4, Iz=833333mm4, K=100mm4)...")
A = 100 * 10            # Area
Iy = 100 * 10**3 / 12   # I about y-axis (weak axis)
Iz = 10 * 100**3 / 12   # I about z-axis (strong axis)
K = 100                 # Torsion constant
model.AddSection(A, Iy, Iz, K)

# Add beam elements
print("   Adding beam elements...")
for i in range(n_elem):
    model.add_beam_element(i, i, i + 1, 0, 0, 0.0)
print(f"   Number of elements: {len(model.Elements)}")

print(f"\n   Model summary:")
print(f"     Total DOF: {model.DOFNum()}")
print(f"     Free DOF: {model.FreeDOFNum()}")

# ============================================================
# 2. Create loads
# ============================================================
print("\n2. Creating loads...")

# Apply point load at tip (node 3): P = 1000 N downward (-Y direction)
P = 1000.0
tip_node_id = n_elem  # Last node
print(f"   Applying {P} N downward load at node {tip_node_id}...")

# Create load vector
loads = VectorLoad()

# Create NodeLoad and add to vector
node_load = NodeLoad(tip_node_id, 0, -P, 0)  # Py = -1000 N
loads.append(node_load)

print(f"   Number of load cases: {len(loads)}")

# ============================================================
# 3. Run linear static analysis
# ============================================================
print("\n3. Running linear static analysis...")

# Create solver
solver = FELinearStaticOp(model, loads)

# Compute
solver.Compute()

if solver.Computed():
    print("   Analysis completed successfully!")
else:
    print("   Analysis failed!")
    sys.exit(1)

# ============================================================
# 4. Get and display results
# ============================================================
print("\n4. Results:")

# Get nodal displacements
print("\n   === Nodal Displacements ===")
displacements = solver.GetDisplacements()
print(f"   {'Node':<6} {'Dx [mm]':<12} {'Dy [mm]':<12} {'Dz [mm]':<12} {'Rx [rad]':<12} {'Ry [rad]':<12} {'Rz [rad]':<12}")
print("   " + "-" * 78)

for i, disp in enumerate(displacements):
    print(f"   {i:<6} {disp.Dx():<12.6f} {disp.Dy():<12.6f} {disp.Dz():<12.6f} {disp.Rx():<12.6e} {disp.Ry():<12.6e} {disp.Rz():<12.6e}")

# Theoretical tip deflection for cantilever beam: delta = P*L^3 / (3*E*Iz)
E = 205e3
delta_theory = P * L**3 / (3 * E * Iz)
delta_fem = displacements[tip_node_id].Dy()
print(f"\n   Theoretical tip deflection (Dy): {delta_theory:.6f} mm")
print(f"   FEM tip deflection (Dy):         {delta_fem:.6f} mm")
print(f"   Error: {abs(delta_fem - (-delta_theory)) / delta_theory * 100:.4f} %")

# Get reaction forces using GetReactionData() helper method
print("\n   === Reaction Forces ===")
try:
    reactions = solver.GetReactionData()
    if len(reactions) > 0:
        for i in range(len(reactions)):
            react = reactions[i]
            # NodeLoadData has id, Px(), Py(), Pz(), Mx(), My(), Mz() accessors
            print(f"   Reaction {i}: node_id={react.id}")
            print(f"     Px={react.Px():.2f} N, Py={react.Py():.2f} N, Pz={react.Pz():.2f} N")
            print(f"     Mx={react.Mx():.2f} N-mm, My={react.My():.2f} N-mm, Mz={react.Mz():.2f} N-mm")
    else:
        print("   No reaction forces returned (check fixed DOFs)")
except Exception as e:
    print(f"   Could not get reaction forces: {e}")
    import traceback
    traceback.print_exc()

# Get beam stress at specific locations
print("\n   === Beam Stress (Element 0, mid-span p=0.5) ===")
try:
    stress_data = solver.GetBeamStress(0, 0.5)
    print(f"   Axial force Nx:      {stress_data.Nx:.2f} N")
    print(f"   Shear force Qy:      {stress_data.Qy:.2f} N")
    print(f"   Shear force Qz:      {stress_data.Qz:.2f} N")
    print(f"   Torsion Mx:          {stress_data.Mx:.2f} N-mm")
    print(f"   Bending moment My:   {stress_data.My:.2f} N-mm")
    print(f"   Bending moment Mz:   {stress_data.Mz:.2f} N-mm")
except Exception as e:
    print(f"   Error accessing stress data: {e}")

# Get beam displacement along element
print("\n   === Beam Displacement along Element 2 ===")
print(f"   {'p':<8} {'Dx [mm]':<12} {'Dy [mm]':<12} {'Rz [rad]':<12}")
print("   " + "-" * 44)
for p in [0.0, 0.25, 0.5, 0.75, 1.0]:
    beam_disp = solver.GetBeamDisplace(2, p)
    print(f"   {p:<8.2f} {beam_disp.Dx():<12.6f} {beam_disp.Dy():<12.6f} {beam_disp.Rz():<12.6e}")

print("\n" + "=" * 60)
print("Linear static analysis completed successfully!")
print("=" * 60)
