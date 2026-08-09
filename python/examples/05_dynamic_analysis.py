"""
05_dynamic_analysis.py - Dynamic Time History Analysis

This example demonstrates:
- Creating a cantilever beam structure
- Defining acceleration time history input
- Running dynamic analysis using Newmark-beta method
- Extracting displacement, velocity, and acceleration response

Structure:
                      m (distributed mass)
    Fixed =================== <- Cantilever beam
          ^
          Ground acceleration input (Y-direction)

Cantilever beam subjected to harmonic base motion

"""
import sys
sys.path.insert(0, '..')

from femnet import *
import math

print("=" * 60)
print("FEMNet Example 05: Dynamic Time History Analysis")
print("=" * 60)

# ============================================================
# 1. Create model (Cantilever beam)
# ============================================================
print("\n1. Creating cantilever beam model...")

model = FEModel()

# Geometry parameters
L = 6000.0  # Beam length [mm] - long flexible beam
n_elements = 10  # Number of beam elements

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

# Rectangular section: 50mm x 100mm (more flexible)
b = 50.0  # width [mm]
h = 100.0  # height [mm]
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
total_mass = 5000.0  # 5000 kg total (large lumped mass for low frequency)
mass_per_node = total_mass / n_elements  # Distribute among free nodes

for i in range(1, model.NodeNum()):
    node = model.GetNode(i)
    node.MassData.Mass = mass_per_node

print(f"   Mass per node: {mass_per_node:.0f} kg")
print(f"   Total model mass: {model.SumNodeMass():.0f} kg")

# ============================================================
# 6. First, perform modal analysis to get natural frequency
# ============================================================
print("\n6. Running modal analysis first...")

# Number of modes to extract
n_modes = 1

# Solve eigenvalue problem (operator owns the computation)
vib = FEVibrationAnalysis(model)
result = vib.Compute(n_modes)
eigen_values = vib.EigenValues()

if result > 0:
    omega_1 = eigen_values[0]
    f_1 = omega_1 / (2 * math.pi)
    T_1 = 1 / f_1
    print(f"   1st mode: omega = {omega_1:.4f} rad/s, f = {f_1:.4f} Hz, T = {T_1:.4f} s")
else:
    print(f"   Modal analysis failed with error code: {result}")
    sys.exit(1)

# ============================================================
# 7. Define dynamic input (harmonic acceleration)
# ============================================================
print("\n7. Defining dynamic input...")

# Harmonic ground motion parameters
duration = 2.0  # seconds
dt = 0.01  # time step [s]
n_steps = int(duration / dt)
excitation_freq = f_1 * 0.9  # Near resonance frequency (90% of natural freq)
amplitude = 500.0  # mm/s^2

# Generate harmonic acceleration time history
accel_data = VectorDouble()
for i in range(n_steps):
    t = i * dt
    # Ramped harmonic: A * sin(omega*t) * (1 - exp(-t/0.5))
    omega_excitation = 2 * math.pi * excitation_freq
    ramp = 1 - math.exp(-t / 0.2)  # Ramp to avoid initial shock
    a = amplitude * math.sin(omega_excitation * t) * ramp
    accel_data.append(a)

# Create DynamicAccelLoad (direction: Y-direction)
accel_load = DynamicAccelLoad(dt * 1000, 0, 1, 0)  # dt in ms, direction vector
accel_load.Accels = accel_data

print(f"   Duration: {duration:.1f} s")
print(f"   Time step: {dt*1000:.0f} ms")
print(f"   Excitation frequency: {excitation_freq:.4f} Hz (90%% of natural)")
print(f"   Peak acceleration: {amplitude:.1f} mm/s^2")

# ============================================================
# 8. Run dynamic analysis
# ============================================================
print("\n8. Running dynamic analysis...")

# Create damping initializer (5% damping)
damp_init = FEDynamicStiffDampInitializer(0.05)

# Create dynamic analysis object
dynamic_solver = DynamicAnalysis(model, accel_load, damp_init)

# Initialize
if not dynamic_solver.Initialize():
    print("   Failed to initialize dynamic analysis!")
    sys.exit(1)

print(f"   Initialized successfully")

# Run simulation
print(f"   Computing {n_steps} time steps...")

# Store response history at tip node
tip_node_id = n_elements  # Last node (tip)
time_history = []
disp_history = []
vel_history = []
accel_history = []

for step in range(n_steps):
    dynamic_solver.ComputeStep()

    disps = dynamic_solver.GetDisplacements()
    vels = dynamic_solver.GetVelocities()
    accels = dynamic_solver.GetAccelerations()

    time_history.append(step * dt)
    disp_history.append(disps[tip_node_id].Dy())  # Tip node Y displacement
    vel_history.append(vels[tip_node_id].Dy())     # Tip node Y velocity
    accel_history.append(accels[tip_node_id].Dy()) # Tip node Y acceleration

print(f"   Analysis completed!")

# ============================================================
# 9. Display results
# ============================================================
print("\n9. Results:")

# Summary statistics
max_disp = max(abs(d) for d in disp_history)
max_vel = max(abs(v) for v in vel_history)
max_accel = max(abs(a) for a in accel_history)

print("\n   === Response Summary (at tip node) ===")
print(f"   Max displacement: {max_disp:.4f} mm")
print(f"   Max velocity: {max_vel:.4f} mm/s")
print(f"   Max acceleration: {max_accel:.4f} mm/s^2")

# Sample of time history
print("\n   === Time History Sample (every 20 steps) ===")
print(f"   {'Time [s]':<12} {'Disp [mm]':<14} {'Vel [mm/s]':<14} {'Accel [mm/s^2]':<14}")
print("   " + "-" * 55)

for i in range(0, len(time_history), 20):
    t = time_history[i]
    d = disp_history[i]
    v = vel_history[i]
    a = accel_history[i]
    print(f"   {t:<12.3f} {d:<14.4f} {v:<14.4f} {a:<14.4f}")

# Dynamic magnification factor
print("\n   === Dynamic Analysis Summary ===")
print(f"   Excitation frequency ratio (r): {excitation_freq / f_1:.2f}")
print(f"   Damping ratio: 5%")

# For SDOF approximation: DMF = 1 / sqrt((1-r^2)^2 + (2*zeta*r)^2)
r = excitation_freq / f_1
zeta = 0.05
DMF = 1 / math.sqrt((1 - r**2)**2 + (2 * zeta * r)**2)
print(f"   Dynamic Magnification Factor (approx.): {DMF:.2f}")

# For SDOF approximation of cantilever with tip mass:
# omega^2 = 3*E*I / (m_eff * L^3) where m_eff ~ 0.23*beam_mass + tip_mass
# Pseudo-static displacement = F / k = (m_eff * a) / (omega^2 * m_eff) = a / omega^2
I = Iz
static_disp_estimate = amplitude / (omega_1**2)  # a/omega^2
print(f"   Pseudo-static tip displacement: {static_disp_estimate:.4f} mm")
print(f"   Max dynamic displacement: {max_disp:.4f} mm")
print(f"   Dynamic amplification ratio: {max_disp / static_disp_estimate:.2f}")

print("\n" + "=" * 60)
print("Example 05 completed successfully!")
print("=" * 60)
