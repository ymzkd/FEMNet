"""
07_plate_element.py - Plate Bending Analysis

This example demonstrates:
- Creating plate elements (combined membrane + bending)
- Modeling a simply supported square plate
- Applying uniform pressure loading
- Running linear static analysis
- Extracting bending displacements

Structure:
    Simply supported edges (Uz = 0)
    +--------+--------+--------+--------+
    |        |        |        |        |
    |   e0   |   e1   |   e2   |   e3   |
    +--------+--------+--------+--------+
    |        |        |        |        |  Uniform pressure q
    |   e4   |   e5   |   e6   |   e7   |     (downward)
    +--------+--------+--------+--------+
    |        |        |        |        |
    |   e8   |   e9   |  e10   |  e11   |
    +--------+--------+--------+--------+
    |        |        |        |        |
    |  e12   |  e13   |  e14   |  e15   |
    +--------+--------+--------+--------+

Square plate under uniform pressure with simply supported edges

"""
import sys
sys.path.insert(0, '..')

from femnet import *
import math

print("=" * 60)
print("FEMNet Example 07: Plate Bending Analysis")
print("=" * 60)

# ============================================================
# 1. Create model and define geometry
# ============================================================
print("\n1. Creating plate model...")

model = FEModel()

# Geometry parameters
a = 1000.0  # Plate side length [mm] (square plate)
thickness = 10.0  # Plate thickness [mm]

n_elem = 4  # Number of elements per side

dx = a / n_elem
dy = a / n_elem

# Create nodes in a grid
print("   Adding nodes...")
node_id = 0
for j in range(n_elem + 1):
    for i in range(n_elem + 1):
        x = i * dx
        y = j * dy
        model.AddNode(node_id, x, y, 0)
        node_id += 1

print(f"   Number of nodes: {model.NodeNum()}")

# ============================================================
# 2. Define supports (boundary conditions)
# ============================================================
print("\n2. Setting boundary conditions...")

# Simply supported plate:
# - All edge nodes: Fix Uz (transverse displacement) = 0
# - One corner: Fix Ux, Uy (prevent rigid body motion)
# - Interior nodes: All in-plane DOFs fixed, Uz and rotations Rx, Ry free

for i in range(model.NodeNum()):
    node = model.GetNode(i)
    x = node.Location.x
    y = node.Location.y

    # Check if node is on the boundary
    on_left = abs(x) < 1e-6
    on_right = abs(x - a) < 1e-6
    on_bottom = abs(y) < 1e-6
    on_top = abs(y - a) < 1e-6
    on_boundary = on_left or on_right or on_bottom or on_top

    if on_boundary:
        if abs(x) < 1e-6 and abs(y) < 1e-6:
            # Bottom left corner: Pin (Ux, Uy, Uz fixed)
            node.Fix = Support(True, True, True, False, False, True)
        elif abs(x - a) < 1e-6 and abs(y) < 1e-6:
            # Bottom right corner: Roller in X (Uy, Uz fixed)
            node.Fix = Support(False, True, True, False, False, True)
        else:
            # Other boundary nodes: Simply supported (Uz fixed)
            # Ux, Uy free for membrane behavior
            # Rx, Ry free for rotation
            # Rz fixed (drilling rotation)
            node.Fix = Support(False, False, True, False, False, True)
    else:
        # Interior nodes: Only fix drilling rotation (Rz)
        # Ux, Uy, Uz, Rx, Ry are free
        node.Fix = Support(False, False, False, False, False, True)

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
# 4. Create quadrilateral plate elements
# ============================================================
print("\n4. Creating plate elements...")

elem_id = 0
for j in range(n_elem):
    for i in range(n_elem):
        # Node indices for this element (counterclockwise)
        n0_idx = j * (n_elem + 1) + i
        n1_idx = j * (n_elem + 1) + (i + 1)
        n2_idx = (j + 1) * (n_elem + 1) + (i + 1)
        n3_idx = (j + 1) * (n_elem + 1) + i

        # Create QuadPlateElement using index-based method
        model.add_quad_plate_element(elem_id, n0_idx, n1_idx, n2_idx, n3_idx, thickness, 0)
        elem_id += 1

print(f"   Number of elements: {len(model.Elements)}")

# ============================================================
# 5. Apply uniform pressure load
# ============================================================
print("\n5. Applying uniform pressure load...")

# Uniform pressure on the plate
q = 0.01  # Pressure [N/mm^2] = 10 kPa

# For plate elements, we can convert area load to nodal loads
# Total force = q * a * a = q * a^2
total_force = q * a * a
print(f"   Uniform pressure: {q} N/mm^2 ({q*1000:.0f} kPa)")
print(f"   Total force: {total_force/1000:.1f} kN")

# Apply equivalent nodal forces (simplified: distribute to all nodes)
# For a uniform mesh, each node gets force proportional to tributary area
loads = VectorLoad()

n_nodes = (n_elem + 1) ** 2
for i in range(model.NodeNum()):
    node = model.GetNode(i)
    x = node.Location.x
    y = node.Location.y

    # Check if node is on corner, edge, or interior
    on_left = abs(x) < 1e-6
    on_right = abs(x - a) < 1e-6
    on_bottom = abs(y) < 1e-6
    on_top = abs(y - a) < 1e-6

    # Corner nodes: tributary area = dx*dy/4
    # Edge nodes: tributary area = dx*dy/2
    # Interior nodes: tributary area = dx*dy

    corners = (on_left or on_right) and (on_bottom or on_top)
    edges = (on_left or on_right or on_bottom or on_top) and not corners

    if corners:
        trib_area = dx * dy / 4
    elif edges:
        trib_area = dx * dy / 2
    else:
        trib_area = dx * dy

    Fz = -q * trib_area  # Negative for downward load

    node_load = NodeLoad(i, 0, 0, Fz)
    loads.append(node_load)

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

# Nodal displacements - show center and key points
print("\n   === Key Nodal Displacements ===")
displacements = solver.GetDisplacements()

# Find center node
center_idx = (n_elem // 2) * (n_elem + 1) + (n_elem // 2)

print(f"   {'Location':<20} {'Dz [mm]':<12} {'Rx [rad]':<12} {'Ry [rad]':<12}")
print("   " + "-" * 60)

# Center node
center_node = model.GetNode(center_idx)
center_disp = displacements[center_idx]
print(f"   {'Center':<20} {center_disp.Dz():<12.6f} {center_disp.Rx():<12.6f} {center_disp.Ry():<12.6f}")

# Mid-edge nodes
mid_bottom_idx = n_elem // 2
mid_top_idx = n_elem * (n_elem + 1) + n_elem // 2
mid_left_idx = (n_elem // 2) * (n_elem + 1)
mid_right_idx = (n_elem // 2) * (n_elem + 1) + n_elem

for name, idx in [("Mid-bottom", mid_bottom_idx), ("Mid-top", mid_top_idx),
                   ("Mid-left", mid_left_idx), ("Mid-right", mid_right_idx)]:
    disp = displacements[idx]
    print(f"   {name:<20} {disp.Dz():<12.6f} {disp.Rx():<12.6f} {disp.Ry():<12.6f}")

# Maximum displacement
max_dz = max(abs(d.Dz()) for d in displacements)
print(f"\n   Maximum deflection |Dz|: {max_dz:.6f} mm")

# Theoretical comparison for simply supported square plate
# Maximum deflection at center: w_max = alpha * q * a^4 / D
# where D = E * t^3 / (12 * (1 - nu^2))
# For simply supported square plate, alpha = 0.00406 (Roark's formula)

print("\n   === Theoretical Comparison ===")
D = E * thickness**3 / (12 * (1 - nu**2))  # Flexural rigidity
alpha = 0.00406  # Coefficient for simply supported square plate

w_max_theory = alpha * q * a**4 / D

print(f"   Plate flexural rigidity D: {D:.2e} N-mm")
print(f"   Theoretical max deflection (Roark): {w_max_theory:.6f} mm")
print(f"   FEM max deflection: {abs(center_disp.Dz()):.6f} mm")
print(f"   Error: {abs(abs(center_disp.Dz()) - w_max_theory) / w_max_theory * 100:.2f}%")
print(f"   (Note: Coarse mesh may result in some error)")

# Equilibrium check
print("\n   === Equilibrium Check ===")
reactions = solver.GetReactionData()
sum_rz = sum(r.Pz() for r in reactions)
print(f"   Sum of Rz reactions: {sum_rz:.2f} N (upward)")
print(f"   Applied force: {-total_force:.2f} N (downward)")
print(f"   Equilibrium: {abs(sum_rz - total_force):.6f} N (should be ~0)")

print("\n" + "=" * 60)
print("Example 07 completed successfully!")
print("=" * 60)
