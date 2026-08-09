"""
08_response_spectrum.py - Response Spectrum Analysis

This example demonstrates:
- Creating a multi-story frame structure
- Running modal analysis to get natural frequencies and mode shapes
- Defining a design response spectrum
- Running response spectrum analysis using SRSS combination
- Extracting maximum displacements and forces

Structure:
    Level 3  o===o===o===o
             |   |   |   |
    Level 2  o===o===o===o
             |   |   |   |
    Level 1  o===o===o===o
             |   |   |   |
    Base     X===X===X===X  (Fixed supports)

3-story, 3-bay frame under seismic response spectrum loading

"""
import sys
sys.path.insert(0, '..')

from femnet import *
import math

print("=" * 60)
print("FEMNet Example 08: Response Spectrum Analysis")
print("=" * 60)

# ============================================================
# 1. Define a custom Response Spectrum class
# ============================================================

class DesignSpectrum(IResponseSpectrum):
    """
    Simple design response spectrum based on typical code values.

    Spectral acceleration Sa = f(T) where T is natural period
    - Constant acceleration plateau: T < T1
    - Descending branch: T1 < T < T2 (1/T relationship)
    - Long period descending: T > T2 (1/T^2 relationship)
    """

    def __init__(self, Sa_max=2.5, T1=0.3, T2=2.0, Sa_base=0.4):
        """
        Initialize design spectrum parameters.

        Args:
            Sa_max: Maximum spectral acceleration ratio (g)
            T1: End of constant acceleration plateau (s)
            T2: End of 1/T descending branch (s)
            Sa_base: Base acceleration at T=0 (g)
        """
        super().__init__()
        self.Sa_max = Sa_max
        self.T1 = T1
        self.T2 = T2
        self.Sa_base = Sa_base
        self.g = 9806.65  # mm/s^2

    def Acceleration(self, T):
        """Return spectral acceleration in mm/s^2 for period T (seconds)."""
        if T <= 0:
            return self.Sa_base * self.g
        elif T < self.T1:
            # Linear increase from Sa_base to Sa_max
            return (self.Sa_base + (self.Sa_max - self.Sa_base) * T / self.T1) * self.g
        elif T < self.T2:
            # Constant acceleration plateau then 1/T descending
            return self.Sa_max * (self.T1 / T) * self.g
        else:
            # 1/T^2 descending branch
            return self.Sa_max * self.T1 * self.T2 / (T * T) * self.g

    def Velocity(self, T):
        """Return spectral velocity in mm/s for period T (seconds)."""
        # Sv = Sa * T / (2*pi)
        Sa = self.Acceleration(T)
        return Sa * T / (2 * math.pi)

    def Displacement(self, T):
        """Return spectral displacement in mm for period T (seconds)."""
        # Sd = Sa * (T / (2*pi))^2
        Sa = self.Acceleration(T)
        return Sa * (T / (2 * math.pi)) ** 2


# ============================================================
# 2. Create frame model
# ============================================================
print("\n1. Creating 3-story frame model...")

model = FEModel()

# Frame geometry
n_stories = 3
n_bays = 3
story_height = 3500.0  # mm
bay_width = 6000.0     # mm

# Create nodes
print("   Adding nodes...")
node_id = 0
for level in range(n_stories + 1):  # 0 to n_stories
    for col in range(n_bays + 1):   # 0 to n_bays
        x = col * bay_width
        y = level * story_height
        model.AddNode(node_id, x, y, 0)
        node_id += 1

n_nodes_per_level = n_bays + 1
print(f"   Number of nodes: {model.NodeNum()}")

# ============================================================
# 3. Define supports (boundary conditions)
# ============================================================
print("\n2. Setting boundary conditions...")

# Base nodes are fully fixed (level 0)
for col in range(n_bays + 1):
    node_idx = col  # Level 0 nodes
    model.GetNode(node_idx).Fix.FixAll()

# All other nodes: constrain out-of-plane DOFs only
# Free: Ux, Uy, Rz (in-plane DOFs for 2D frame)
for i in range(n_nodes_per_level, model.NodeNum()):
    node = model.GetNode(i)
    # Fix: Uz, Rx, Ry (out-of-plane)
    # Free: Ux, Uy, Rz (in-plane for frame)
    node.Fix = Support(False, False, True, True, True, False)

print(f"   Total DOF: {model.DOFNum()}")
print(f"   Free DOF: {model.FreeDOFNum()}")

# ============================================================
# 4. Define material and section
# ============================================================
print("\n3. Adding material and section...")

# Steel material with density
E = 205e3       # Young's modulus [N/mm^2]
nu = 0.3        # Poisson's ratio
rho = 7.85e-9   # Density [kg/mm^3]
model.AddMaterialWithDensity(E, nu, rho)

# Column section: H-300x300
A_col = 11730.0
Iy_col = 2.04e8
Iz_col = 6.76e7
K_col = 1.03e6
model.AddSection(A_col, Iy_col, Iz_col, K_col)  # Section 0: Columns

# Beam section: H-400x200
A_beam = 8410.0
Iy_beam = 2.37e8
Iz_beam = 1.74e7
K_beam = 4.45e5
model.AddSection(A_beam, Iy_beam, Iz_beam, K_beam)  # Section 1: Beams

print(f"   Column: H-300x300 (A={A_col:.0f} mm^2)")
print(f"   Beam: H-400x200 (A={A_beam:.0f} mm^2)")

# ============================================================
# 5. Create beam/column elements
# ============================================================
print("\n4. Creating frame elements...")

elem_id = 0

# Columns (vertical elements)
for level in range(n_stories):
    for col in range(n_bays + 1):
        n1_idx = level * n_nodes_per_level + col
        n2_idx = (level + 1) * n_nodes_per_level + col
        model.add_beam_element(elem_id, n1_idx, n2_idx, 0, 0, 0.0)  # Section 0
        elem_id += 1

# Beams (horizontal elements)
for level in range(1, n_stories + 1):
    for bay in range(n_bays):
        n1_idx = level * n_nodes_per_level + bay
        n2_idx = level * n_nodes_per_level + bay + 1
        model.add_beam_element(elem_id, n1_idx, n2_idx, 1, 0, 0.0)  # Section 1
        elem_id += 1

print(f"   Number of elements: {len(model.Elements)}")

# ============================================================
# 6. Add floor masses
# ============================================================
print("\n5. Adding floor masses...")

# Typical floor mass (distributed to nodes)
floor_mass = 50000.0  # kg per floor (dead + live load)
mass_per_node = floor_mass / n_nodes_per_level

for level in range(1, n_stories + 1):
    for col in range(n_nodes_per_level):
        node_idx = level * n_nodes_per_level + col
        model.GetNode(node_idx).MassData.Mass = mass_per_node

print(f"   Mass per floor: {floor_mass:.0f} kg")
print(f"   Total model mass: {model.SumNodeMass():.0f} kg")

# ============================================================
# 7. Run modal analysis
# ============================================================
print("\n6. Running modal analysis...")

n_modes = 5  # Number of modes to extract

vib = FEVibrationAnalysis(model)
result = vib.Compute(n_modes)
eigen_values = vib.EigenValues()
mode_vectors = vib.ModeVectors()

if result > 0:
    print(f"   Extracted {result} modes successfully!")
else:
    print(f"   Modal analysis failed with error code: {result}")
    sys.exit(1)

# Display natural periods
print("\n   === Natural Frequencies ===")
print(f"   {'Mode':<6} {'omega [rad/s]':<15} {'f [Hz]':<12} {'T [s]':<12}")
print("   " + "-" * 50)

for i in range(len(eigen_values)):
    omega = eigen_values[i]
    f = omega / (2 * math.pi)
    T = 1 / f if f > 0 else float('inf')
    print(f"   {i+1:<6} {omega:<15.4f} {f:<12.4f} {T:<12.4f}")

# ============================================================
# 8. Create FEVibrationAnalysis and response spectrum
# ============================================================
print("\n7. Setting up response spectrum analysis...")

# Create FEVibrationAnalysis object
vibrate_result = FEVibrationAnalysis(model, mode_vectors, eigen_values)

# Get natural periods
periods = vibrate_result.NaturalPeriods()
print(f"   Natural periods: {[f'{T:.4f}' for T in periods]}")

# Create design spectrum
spectrum = DesignSpectrum(Sa_max=0.4, T1=0.3, T2=2.0, Sa_base=0.15)
print(f"   Design spectrum: Sa_max={spectrum.Sa_max}g, T1={spectrum.T1}s, T2={spectrum.T2}s")

# Show spectral acceleration at each mode period
print("\n   === Spectral Accelerations ===")
for i, T in enumerate(periods):
    Sa = spectrum.Acceleration(T)
    print(f"   Mode {i+1}: T = {T:.4f} s, Sa = {Sa:.2f} mm/s^2 ({Sa/spectrum.g:.3f}g)")

# ============================================================
# 9. Run response spectrum analysis
# ============================================================
print("\n8. Running response spectrum analysis...")

# Create response spectrum method (SRSS combination)
direction = Vector(1.0, 0.0, 0.0)  # X-direction earthquake

# Use SRSS method for combining modal responses
rs_method = ResponseSpectrumMethod(model, vibrate_result, direction, spectrum, SRSS)

# Compute
rs_method.Compute()

if rs_method.Computed():
    print("   Response spectrum analysis completed!")
else:
    print("   Response spectrum analysis failed!")
    sys.exit(1)

# ============================================================
# 10. Display results
# ============================================================
print("\n9. Results:")

# Get displacements
displacements = rs_method.GetDisplacements()

# Show story displacements (at first column)
print("\n   === Story Displacements (at left column) ===")
print(f"   {'Level':<8} {'Height [m]':<12} {'Dx [mm]':<12} {'Dy [mm]':<12}")
print("   " + "-" * 45)

for level in range(n_stories + 1):
    node_idx = level * n_nodes_per_level  # First column
    node = model.GetNode(node_idx)
    disp = displacements[node_idx]
    print(f"   {level:<8} {node.Location.y/1000:<12.1f} {disp.Dx():<12.4f} {disp.Dy():<12.4f}")

# Maximum displacements
print("\n   === Maximum Response ===")
max_dx = max(d.Dx() for d in displacements)
max_dy = max(d.Dy() for d in displacements)
print(f"   Max Dx (horizontal): {max_dx:.4f} mm")
print(f"   Max Dy (vertical): {max_dy:.4f} mm")

# Story drift (inter-story displacement)
print("\n   === Inter-Story Drifts ===")
print(f"   {'Story':<8} {'Drift [mm]':<12} {'Drift ratio':<12}")
print("   " + "-" * 35)

for level in range(1, n_stories + 1):
    node_idx_upper = level * n_nodes_per_level
    node_idx_lower = (level - 1) * n_nodes_per_level

    dx_upper = displacements[node_idx_upper].Dx()
    dx_lower = displacements[node_idx_lower].Dx()

    drift = dx_upper - dx_lower
    drift_ratio = drift / story_height

    print(f"   {level:<8} {drift:<12.4f} {drift_ratio:<12.6f}")

# Effective mass participation
print("\n   === Modal Participation (X-direction) ===")
participation = vibrate_result.ParticipationFactors(Vector(1.0, 0.0, 0.0))
mass_rates = vibrate_result.EffectiveMassRates()

print(f"   {'Mode':<6} {'Part. Factor':<15} {'Eff. Mass %':<15}")
print("   " + "-" * 40)
total_mass_rate = 0.0
for i in range(len(participation)):
    pf = participation[i]
    mr = mass_rates[i] * 100  # Convert to percentage
    total_mass_rate += mr
    print(f"   {i+1:<6} {pf:<15.4f} {mr:<15.2f}")

print(f"\n   Total effective mass ratio: {total_mass_rate:.1f}%")

print("\n" + "=" * 60)
print("Example 08 completed successfully!")
print("=" * 60)
