"""Detailed test of LTI stability analysis for documentation."""

import numpy as np
from src.stability_analysis.lti_stability import LTIStability
from src.rotor_build import RotorBuild

print("=" * 70)
print("LTI STABILITY ANALYSIS TEST")
print("=" * 70)

# Build rotor system
print("\n1. Building rotor system...")
rotor = RotorBuild.build_all()
rotor.problem.required_solver = "LTI"
print("   [OK] Rotor system built successfully")
print(f"   - Number of analysis points: {rotor.problem.number_points}")

# Run stability analysis
print("\n2. Running LTI stability analysis...")
lti = LTIStability(rotor)
lti = lti.eigen_full_range()
print("   [OK] Analysis completed successfully")

# Display results
print("\n3. Analysis Results:")
print(f"   - Total number of solutions: {len(lti.modal_solution)}")
print(f"   - Number of eigenvalues per solution: {len(lti.modal_solution[0].damping)}")

# First point details
print(f"\n4. First Analysis Point (i=0):")
print(f"   - OMEGA: {lti.modal_solution[0].OMEGA_RPM:.2f} RPM ({lti.modal_solution[0].OMEGA:.4f} rad/s)")
print(f"   - Number of modes: {len(lti.modal_solution[0].damping)}")

# Last point details
print(f"\n5. Last Analysis Point (i={len(lti.modal_solution)-1}):")
last_idx = len(lti.modal_solution) - 1
print(f"   - OMEGA: {lti.modal_solution[last_idx].OMEGA_RPM:.2f} RPM ({lti.modal_solution[last_idx].OMEGA:.4f} rad/s)")

# Eigenvalue statistics at first point
print(f"\n6. Eigenvalue Statistics at First Point:")
damping = lti.modal_solution[0].damping
frequency = lti.modal_solution[0].frequency
print(f"   Damping:")
print(f"   - Min: {np.min(damping):.6e}")
print(f"   - Max: {np.max(damping):.6e}")
print(f"   - Mean: {np.mean(damping):.6e}")
print(f"   Frequency:")
print(f"   - Min: {np.min(np.abs(frequency)):.6e}")
print(f"   - Max: {np.max(np.abs(frequency)):.6e}")
print(f"   - Mean: {np.mean(np.abs(frequency)):.6e}")

# Sample eigenvalues
print(f"\n7. Sample Eigenvalues (First 3 modes at first point):")
damping_flat = damping.flatten()
frequency_flat = frequency.flatten()
for i in range(min(3, len(damping_flat))):
    eigenvalue = complex(damping_flat[i], frequency_flat[i])
    print(f"   Mode {i}: {eigenvalue}")

# Check stability
print(f"\n8. Stability Check at First Point:")
stable_modes = np.sum(damping < 0)
unstable_modes = np.sum(damping > 0)
print(f"   - Stable modes (damping < 0): {stable_modes}")
print(f"   - Unstable modes (damping > 0): {unstable_modes}")

# Check for critical speeds
print(f"\n9. Checking for instabilities across speed range...")
instability_points = []
for i, sol in enumerate(lti.modal_solution):
    if np.any(sol.damping > 0):
        instability_points.append(i)

if instability_points:
    print(f"   [WARNING] Found instabilities at {len(instability_points)} points:")
    for idx in instability_points[:5]:  # Show first 5
        print(f"   - Point {idx}: {lti.modal_solution[idx].OMEGA_RPM:.2f} RPM")
    if len(instability_points) > 5:
        print(f"   ... and {len(instability_points) - 5} more")
else:
    print(f"   [OK] No instabilities detected in the speed range")

print("\n" + "=" * 70)
print("TEST COMPLETED SUCCESSFULLY")
print("=" * 70)
