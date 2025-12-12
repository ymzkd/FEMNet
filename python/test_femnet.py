"""
FEMNet Python bindings test script
"""
import sys
sys.path.insert(0, '.')

from femnet import *

print("=" * 60)
print("FEMNet Python Bindings Test")
print("=" * 60)

# Create a simple beam model
print("\n1. Creating FEM model...")
model = FEModel()

# Add nodes using the helper method
print("   Adding nodes...")
model.AddNode(0, 0, 0, 0)
model.AddNode(1, 1000, 0, 0)
model.AddNode(2, 2000, 0, 0)
model.AddNode(3, 3000, 0, 0)

print(f"   Number of nodes: {model.NodeNum()}")

# Set boundary conditions (fix first node)
print("   Setting boundary conditions...")
model.GetNode(0).Fix.FixAll()

# Add material (Material has default constructor so append works)
print("   Adding material...")
mat = model.AddMaterial(205e3, 0.3)  # Steel: E=205GPa, nu=0.3

# Add section using helper method
print("   Adding section...")
model.AddSection(100*10, 1000, 500, 100)  # A, Iy, Iz, K

# Add beam elements
print("   Adding beam elements...")
for i in range(3):
    model.add_beam_element(
        i,          # element id
        i,          # node i
        i + 1,      # node j
        0,          # section id
        0,          # material id
        0.0         # beta angle
    )

print(f"   Number of elements: {len(model.Elements)}")

# Create loads
print("   Creating loads...")
loads = VectorLoad()
node_load = NodeLoad(3, 0, -1000, 0)  # 1000N downward at node 3
# NodeLoad needs to be wrapped in shared_ptr for VectorLoad
# For now, let's try a different approach using make_shared

# Create linear static operator using shared model
print("\n2. Creating shared model pointer...")
shared_model = FEModel()

# Copy nodes
for i in range(4):
    shared_model.AddNode(i, i * 1000, 0, 0)
shared_model.GetNode(0).Fix.FixAll()

# Add material and section
shared_model.AddMaterial(205e3, 0.3)
shared_model.AddSection(100*10, 1000, 500, 100)

# Add elements
for i in range(3):
    shared_model.add_beam_element(i, i, i + 1, 0, 0, 0.0)

print("   Model setup complete")
print(f"   Nodes: {shared_model.NodeNum()}")
print(f"   DOF: {shared_model.DOFNum()}")
print(f"   Free DOF: {shared_model.FreeDOFNum()}")

# Test basic functionality
print("\n3. Testing basic functionality...")
print("   Testing Material:")
test_mat = Material(205e3, 0.3)
print(f"   Young's modulus: {test_mat.Young}")
print(f"   Poisson's ratio: {test_mat.Poisson}")
print(f"   Shear modulus G: {test_mat.G()}")

print("\n   Testing Displacement:")
disp = Displacement(1.0, 2.0, 3.0, 0.1, 0.2, 0.3)
print(f"   Dx: {disp.Dx()}, Dy: {disp.Dy()}, Dz: {disp.Dz()}")
print(f"   Rx: {disp.Rx()}, Ry: {disp.Ry()}, Rz: {disp.Rz()}")

print("\n   Testing Vector:")
v1 = Vector(1, 2, 3)
print(f"   Vector: ({v1.x}, {v1.y}, {v1.z})")
print(f"   Norm: {v1.norm()}")

print("\n   Testing Point:")
p1 = Point(0, 0, 0)
p2 = Point(3, 4, 0)
print(f"   Point 1: ({p1.x}, {p1.y}, {p1.z})")
print(f"   Point 2: ({p2.x}, {p2.y}, {p2.z})")
print(f"   Distance: {p1.distance_to(p2)}")

print("\n   Testing NodeLoadData:")
nld = NodeLoadData(1, 100, 200, 300)
print(f"   NodeLoadData ID: {nld.id}")

print("\n" + "=" * 60)
print("Basic tests completed successfully!")
print("=" * 60)
