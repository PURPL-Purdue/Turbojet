# Main script?
import pint
from utils.ThermoState import ThermoState
from component_models.compressor import Compressor
ureg = pint.UnitRegistry()
# Steps (at a GIVEN RPM):
# 1. Pick a thermo state & RPM
# - Total pressure
# - Total temperature
# - Massflow (guess)
# 2. Go to compressor map: give it the thermo state, have it give you a power and outlet flow state.
# 3. Go to combustor, same
# 4. Go to turbine map, same.
# 5. Check the power balance, if not balanced, adjust massflow guess and repeat.

# INPUTS:
min_rpm = 100 * ureg.rpm;
max_rpm = 100 * 1000 * ureg.rpm;
start_massflow = 0.01 * ureg.kg / ureg.s;
P_amb = 101325 * ureg.pascal;
T_amb = 288.15 * ureg.kelvin;

#1: Step 1: Pick a thermo-state and RPM.
inlet_state = ThermoState(P=P_amb, T=T_amb, mdot=start_massflow);

#2: Go to the compressor map, and with the thermo state get an outlet state and power.
# Create compressor with SCALED map to custom design point
compressor = Compressor(
    "reference/G47_G50_109_CHICS.csv",
    "reference/G47_G50_109_ISOEFF.csv",
    design_pr=3.0,              # Target design PR
    design_efficiency=0.70,     # Target max efficiency
    design_rpm=80000,           # Target design RPM
    design_mass_flow=0.485,       # Target design mass flow [kg/s]
)
print("\n=== Scaled Map ===")
compressor.print_scaling_info()
compressor.plot_map()

# Example: Power vs RPM with LINEAR operating line (mass_flow ~ RPM)
import numpy as np
import matplotlib.pyplot as plt

# Design point for operating line
design_rpm = 80000  # rev/min
design_mdot = 0.485  # kg/s

# Linear operating line: mdot = k * rpm, where k = design_mdot / design_rpm
k_operating_line = design_mdot / design_rpm

# Sweep RPM from 40k to 100k
rpm_values = np.linspace(40000, 100000, 30)
power_values = []
pr_values = []
eff_values = []
mdot_values = []

for rpm in rpm_values:
    # Calculate mass flow from linear operating line
    mdot = k_operating_line * rpm  # kg/s
    mdot_values.append(mdot)

    # Create inlet state with this mass flow
    inlet_state = ThermoState(P=P_amb, T=T_amb, mdot=mdot * ureg.kg / ureg.s)

    # Query compressor
    rpm_qty = rpm * ureg.rpm
    outlet_state, power = compressor.calculate(inlet_state, rpm_qty)

    power_values.append(power.to(ureg.kilowatt).magnitude)
    pr_values.append(outlet_state.P_total.magnitude / inlet_state.P_total.magnitude)

    # Also get efficiency
    pr, eff = compressor.query(rpm, mdot)
    eff_values.append(eff)

# Convert to arrays and filter out NaN values
rpm_arr = np.array(rpm_values)
power_arr = np.array(power_values)
pr_arr = np.array(pr_values)
eff_arr = np.array(eff_values)
mdot_arr = np.array(mdot_values)

# Mask valid points
valid = ~np.isnan(power_arr) & ~np.isnan(pr_arr) & ~np.isnan(eff_arr)

# Plot
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

axes[0, 0].plot(rpm_arr[valid] / 1000, power_arr[valid], 'b-o')
axes[0, 0].set_xlabel('RPM [thousands]')
axes[0, 0].set_ylabel('Power [kW]')
axes[0, 0].set_title('Compressor Power vs RPM\n(linear operating line)')
axes[0, 0].grid(True, alpha=0.3)

axes[0, 1].plot(rpm_arr[valid] / 1000, pr_arr[valid], 'r-o')
axes[0, 1].set_xlabel('RPM [thousands]')
axes[0, 1].set_ylabel('Pressure Ratio')
axes[0, 1].set_title('Pressure Ratio vs RPM')
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].plot(rpm_arr[valid] / 1000, eff_arr[valid], 'g-o')
axes[1, 0].set_xlabel('RPM [thousands]')
axes[1, 0].set_ylabel('Efficiency')
axes[1, 0].set_title('Efficiency vs RPM')
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].plot(rpm_arr[valid] / 1000, mdot_arr[valid], 'm-o')
axes[1, 1].set_xlabel('RPM [thousands]')
axes[1, 1].set_ylabel('Mass Flow [kg/s]')
axes[1, 1].set_title('Mass Flow vs RPM (linear operating line)')
axes[1, 1].grid(True, alpha=0.3)

plt.suptitle(f'Compressor Performance Along Linear Operating Line\n(mdot = {k_operating_line:.2e} × RPM)', fontsize=12)
plt.tight_layout()
plt.show()

print(f"\nValid RPM range: {rpm_arr[valid].min()/1000:.1f}k - {rpm_arr[valid].max()/1000:.1f}k")
print(f"Power range: {power_arr[valid].min():.1f} - {power_arr[valid].max():.1f} kW")
print(f"PR range: {pr_arr[valid].min():.2f} - {pr_arr[valid].max():.2f}")
