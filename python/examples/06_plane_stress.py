"""
06_plane_stress.py - 2D Plane Stress Analysis

This example demonstrates:
- Creating 2D quadrilateral plane stress elements
- Modeling a plate with a hole (simplified as rectangular with mesh)
- Applying in-plane tensile loading
- Running linear static analysis
- Extracting membrane stresses

Structure:
    +--------+--------+--------+
    |        |        |        |   ^
    |   e0   |   e1   |   e2   |   | Fy (tension)
    +--------+--------+--------+   |
    |        |        |        |
    |   e3   |   e4   |   e5   |
    +--------+--------+--------+
    ^
    Fixed (Ux, Uy = 0 at bottom left corner, Uy = 0 along bottom edge)

Rectangular plate under uniaxial tension

"""
import sys
sys.path.insert(0, '..')

from femnet import *
import math

print("=" * 60)
print("FEMNet Example 06: 2D Plane Stress Analysis")
print("=" * 60)

# ============================================================
# 1. Create model and define geometry
# ============================================================
print("\n1. Creating plane stress model...")

model = FEModel()

# Geometry parameters
width = 300.0   # Plate width [mm]
height = 200.0  # Plate height [mm]
thickness = 10.0  # Plate thickness [mm]

nx = 3  # Number of elements in x direction
ny = 2  # Number of elements in y direction

dx = width / nx
dy = height / ny

# Create nodes in a grid
print("   Adding nodes...")
node_id = 0
for j in range(ny + 1):
    for i in range(nx + 1):
        x = i * dx
        y = j * dy
        model.AddNode(node_id, x, y, 0)
        node_id += 1

print(f"   Number of nodes: {model.NodeNum()}")

# ============================================================
# 2. Define supports (boundary conditions)
# ============================================================
print("\n2. Setting boundary conditions...")

# For plane stress in XY plane:
# - Uz, Rx, Ry, Rz must be fixed for all nodes (only Ux, Uy free)
# - Bottom left corner: Fix Ux, Uy (prevent rigid body motion)
# - Bottom edge: Fix Uy only (roller supports)

for i in range(model.NodeNum()):
    node = model.GetNode(i)
    x = node.Location.x
    y = node.Location.y

    # All nodes: Fix out-of-plane DOFs (Uz, Rx, Ry, Rz)
    if abs(y) < 1e-6:  # Bottom edge (y = 0)
        if abs(x) < 1e-6:  # Bottom left corner
            # Fix Ux, Uy, Uz, Rx, Ry, Rz (pinned support - prevents rigid body motion)
            node.Fix = Support(True, True, True, True, True, True)
        else:
            # All other bottom nodes: Fix Uy (roller supports)
            # This constrains the bottom edge to remain on the x-axis
            node.Fix = Support(False, True, True, True, True, True)
    else:
        # Interior and top nodes: Free in Ux, Uy
        node.Fix = Support(False, False, True, True, True, True)

print(f"   Total DOF: {model.DOFNum()}")
print(f"   Free DOF: {model.FreeDOFNum()}")

# ============================================================
# 3. Define material
# ============================================================
print("\n3. Adding material...")

# Steel material
E = 205e3  # Young's modulus [N/mm^2]
nu = 0.3   # Poisson's ratio
model.AddMaterial(E, nu)

print(f"   Material: Steel (E={E/1e3:.0f} GPa, nu={nu})")

# ============================================================
# 4. Create quadrilateral plane stress elements
# ============================================================
print("\n4. Creating plane stress elements...")

# Get material
mat = model.Materials[0]

elem_id = 0
for j in range(ny):
    for i in range(nx):
        # Node indices for this element (counterclockwise)
        n0_idx = j * (nx + 1) + i
        n1_idx = j * (nx + 1) + (i + 1)
        n2_idx = (j + 1) * (nx + 1) + (i + 1)
        n3_idx = (j + 1) * (nx + 1) + i

        # Get node references
        n0 = model.GetNode(n0_idx)
        n1 = model.GetNode(n1_idx)
        n2 = model.GetNode(n2_idx)
        n3 = model.GetNode(n3_idx)

        # Create QuadPlaneElement and add to model
        elem = QuadPlaneElement(elem_id, n0, n1, n2, n3, thickness, mat)
        model.add_element(elem)
        elem_id += 1

print(f"   Number of elements: {len(model.Elements)}")

# ============================================================
# 5. Apply loads
# ============================================================
print("\n5. Applying loads...")

# Apply uniform tension on top edge
# Total force distributed among top edge nodes
total_force = 50000.0  # 50 kN total [N]
top_nodes = nx + 1  # Number of nodes on top edge

# Force per node (end nodes get half for uniform distribution)
force_per_node = total_force / (top_nodes - 1)

loads = VectorLoad()

# Apply forces to top edge nodes
for i in range(nx + 1):
    node_idx = ny * (nx + 1) + i

    if i == 0 or i == nx:
        # End nodes get half force
        f = force_per_node / 2
    else:
        # Interior nodes get full force
        f = force_per_node

    node_load = NodeLoad(node_idx, 0, f, 0)  # Fy in +Y direction (tension)
    loads.append(node_load)

print(f"   Total applied force: {total_force/1000:.1f} kN tension in Y direction")
print(f"   Applied stress (nominal): {total_force / (width * thickness):.2f} N/mm^2")

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
print(f"   {'Node':<6} {'x [mm]':<10} {'y [mm]':<10} {'Dx [mm]':<12} {'Dy [mm]':<12}")
print("   " + "-" * 55)

for i in range(model.NodeNum()):
    node = model.GetNode(i)
    disp = displacements[i]
    print(f"   {i:<6} {node.Location.x:<10.1f} {node.Location.y:<10.1f} {disp.Dx():<12.6f} {disp.Dy():<12.6f}")

# Maximum displacement
max_dx = max(abs(d.Dx()) for d in displacements)
max_dy = max(abs(d.Dy()) for d in displacements)
print(f"\n   Max Dx: {max_dx:.6f} mm")
print(f"   Max Dy: {max_dy:.6f} mm")

# Element stresses (theoretical for uniaxial tension)
print("\n   === Theoretical Element Stresses ===")
sigma_y_theory = total_force / (width * thickness)
sigma_x_theory = 0.0  # No stress in X direction for uniaxial tension
tau_xy_theory = 0.0  # No shear stress

print(f"   For uniaxial tension:")
print(f"   Sxx = {sigma_x_theory:.2f} N/mm^2 (MPa)")
print(f"   Syy = {sigma_y_theory:.2f} N/mm^2 (MPa)")
print(f"   Sxy = {tau_xy_theory:.2f} N/mm^2 (MPa)")
print(f"   (Note: Stress extraction from elements requires casting which is not directly available)")

# Theoretical comparison
print("\n   === Theoretical Comparison ===")
sigma_y_theory = total_force / (width * thickness)
eps_y_theory = sigma_y_theory / E
eps_x_theory = -nu * eps_y_theory

# Theoretical displacements
dy_top_theory = eps_y_theory * height
dx_edge_theory = eps_x_theory * width / 2  # From center

print(f"   Applied stress (sigma_y): {sigma_y_theory:.2f} N/mm^2")
print(f"   Axial strain (eps_y): {eps_y_theory:.6f}")
print(f"   Lateral strain (eps_x): {eps_x_theory:.6f}")
print(f"")
print(f"   Theoretical Dy at top: {dy_top_theory:.6f} mm")
print(f"   FEM Dy at top center: {displacements[ny * (nx + 1) + nx // 2].Dy():.6f} mm")
print(f"")
# Lateral contraction at right edge (x = width)
dx_edge_theory_full = eps_x_theory * width
print(f"   Theoretical Dx at right edge: {dx_edge_theory_full:.6f} mm")

# Get Dx at right edge
max_dx_edge = displacements[nx].Dx()  # Node at (300, 0)
print(f"   FEM Dx at right edge: {max_dx_edge:.6f} mm")
print(f"   Error: {abs(max_dx_edge - dx_edge_theory_full) / abs(dx_edge_theory_full) * 100:.2f}%")

# Reaction forces
print("\n   === Reaction Forces ===")
reactions = solver.GetReactionData()
for react in reactions:
    if abs(react.Px()) > 1e-6 or abs(react.Py()) > 1e-6:
        print(f"   Node {react.id}: Px = {react.Px():>10.2f} N, Py = {react.Py():>10.2f} N")

# Equilibrium check
sum_ry = sum(r.Py() for r in reactions)
print(f"\n   Sum of Ry reactions: {sum_ry:.2f} N")
print(f"   Applied force: {total_force:.2f} N")
print(f"   Equilibrium check: {abs(sum_ry + total_force):.6f} N (should be ~0)")

print("\n" + "=" * 60)
print("Example 06 completed successfully!")
print("=" * 60)
