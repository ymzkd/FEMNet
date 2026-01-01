"""
01_simple_truss.py - Simple Truss Structure Analysis

This example demonstrates:
- Creating a simple 2D truss structure
- Adding truss elements (axial force only)
- Applying point loads
- Running linear static analysis
- Extracting displacements and reactions

Structure:

        P
        |
        v
       (2)
       /\
      /  \
     /    \
   (0)-----(1)  <- Fixed supports

"""
import sys
sys.path.insert(0, '..')

from femnet import *

print("=" * 60)
print("FEMNet Example 01: Simple Truss Structure")
print("=" * 60)

# ============================================================
# 1. Create model and define geometry
# ============================================================
print("\n1. Creating truss model...")

model = FEModel()

# Geometry parameters
span = 2000.0    # Horizontal span [mm]
height = 1000.0  # Height [mm]

# Add nodes
print("   Adding nodes...")
model.AddNode(0, 0, 0, 0)           # Node 0: Left support
model.AddNode(1, span, 0, 0)        # Node 1: Right support
model.AddNode(2, span/2, height, 0) # Node 2: Top (load point)

print(f"   Number of nodes: {model.NodeNum()}")

# ============================================================
# 2. Define supports (boundary conditions)
# ============================================================
print("\n2. Setting boundary conditions...")

# For 2D truss in XY plane, we need to:
# 1. Fix support nodes for Ux, Uy (pin supports)
# 2. Fix Uz for all nodes (out-of-plane constraint)
# Note: Rx, Ry, Rz are already locked by default in Support constructor

# Support nodes: Fix all translations (Ux, Uy, Uz)
# The Support constructor takes: (Ux, Uy, Uz, Rx, Ry, Rz)
# Rx, Ry, Rz are locked by default, but we explicitly set them for clarity
model.GetNode(0).Fix = Support(True, True, True, True, True, True)  # Node 0: Fixed support
model.GetNode(1).Fix = Support(True, True, True, True, True, True)  # Node 1: Fixed support

# Free node (node 2): Only fix out-of-plane (Uz) and rotations
# Ux, Uy are free for displacement
model.GetNode(2).Fix = Support(False, False, True, True, True, True)

print(f"   Total DOF: {model.DOFNum()}")
print(f"   Free DOF: {model.FreeDOFNum()}")

# ============================================================
# 3. Define material and section
# ============================================================
print("\n3. Adding material and section...")

# Steel material: E = 205 GPa, nu = 0.3
model.AddMaterial(205e3, 0.3)  # id = 0

# Circular section: diameter = 20mm
d = 20.0
A = 3.14159 * (d/2)**2  # Area = pi*r^2
Iy = 3.14159 * (d/4)**4 / 4  # Not used for truss but required
Iz = Iy
K = 0  # Torsion constant (not used for truss)
model.AddSection(A, Iy, Iz, K)  # id = 0

print(f"   Material: Steel (E=205GPa)")
print(f"   Section: Circular d={d}mm, A={A:.2f}mm2")

# ============================================================
# 4. Create truss elements
# ============================================================
print("\n4. Creating truss elements...")

# add_truss_element(id, n1_id, n2_id, sec_id, mat_id)
model.add_truss_element(0, 0, 2, 0, 0)  # Element 0: Node 0 to 2 (left diagonal)
model.add_truss_element(1, 1, 2, 0, 0)  # Element 1: Node 1 to 2 (right diagonal)
model.add_truss_element(2, 0, 1, 0, 0)  # Element 2: Node 0 to 1 (bottom chord)

print(f"   Number of elements: {len(model.Elements)}")

# ============================================================
# 5. Apply loads
# ============================================================
print("\n5. Applying loads...")

# Apply vertical load at top node
P = -5000.0  # 5000 N downward (negative Y)

loads = VectorLoad()
node_load = NodeLoad(2, 0, P, 0)  # Node 2, Px=0, Py=-5000, Pz=0
loads.append(node_load)

print(f"   Applied load: {-P} N downward at node 2")

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
print(f"   {'Node':<6} {'Dx [mm]':<14} {'Dy [mm]':<14} {'Dz [mm]':<14}")
print("   " + "-" * 50)

for i, disp in enumerate(displacements):
    print(f"   {i:<6} {disp.Dx():<14.6f} {disp.Dy():<14.6f} {disp.Dz():<14.6f}")

# Reaction forces
print("\n   === Reaction Forces ===")
reactions = solver.GetReactionData()
for react in reactions:
    print(f"   Node {react.id}:")
    print(f"     Px = {react.Px():>10.2f} N")
    print(f"     Py = {react.Py():>10.2f} N")
    print(f"     Pz = {react.Pz():>10.2f} N")

# Verify equilibrium
print("\n   === Equilibrium Check ===")
sum_Px = sum(r.Px() for r in reactions)
sum_Py = sum(r.Py() for r in reactions) + P  # Add applied load (P is negative)
sum_Pz = sum(r.Pz() for r in reactions)
print(f"   Sum Px = {sum_Px:.6f} N (should be 0)")
print(f"   Sum Py = {sum_Py:.6f} N (should be 0)")
print(f"   Sum Pz = {sum_Pz:.6f} N (should be 0)")

print("\n" + "=" * 60)
print("Example 01 completed successfully!")
print("=" * 60)
