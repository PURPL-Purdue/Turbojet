# Main script?
import time
import matplotlib.pyplot as plt
from numpy import linspace
from utils.units import ureg
from utils.ThermoState import ThermoState
from component_models.compressor import Compressor
from component_models.combustor import Combustor
from component_models.turbine import Turbine

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
min_rpm = 10 * 1000 * ureg.rpm;
max_rpm = 90 * 1000 * ureg.rpm;
start_massflow = 0.485 * ureg.kg / ureg.s;
P_amb = 101325 * ureg.pascal;
T_amb = 288.15 * ureg.kelvin;

# Check these RPM's:
RPM = linspace(min_rpm.magnitude, max_rpm.magnitude, 50) * ureg.rpm;

# Create compressor with SCALED map to custom design point
compressor = Compressor(
    "reference/G47_G50_109_CHICS.csv",
    "reference/G47_G50_109_ISOEFF.csv",
    design_pr=3.0,              # Target design PR
    design_efficiency=0.60,     # Target max efficiency
    design_rpm=80000,           # Target design RPM
    design_mass_flow=0.485,       # Target design mass flow [kg/s]
)
compressor.print_scaling_info();
# compressor.plot_map();

#Create combustor:
combustor = Combustor();

# Create turbine:
turbine = Turbine(
    design_efficiency=0.70,
    back_pressure=P_amb*2.0,
    design_rpm=80000,
);

# Storage for results
results = {
    'rpm': [],
    'mass_flow': [],
    'pressure_ratio': [],
    'compressor_power': [],
    'turbine_power': [],
    'net_power': [],
}

for rpm in RPM:
    print(f"RPM: {rpm:.1f}")

    power_tol = 500  # 500W tolerance
    net_power = 1000 * ureg.watt  # Initialize net power to enter the loop
    massflow_guess = start_massflow  # Start with initial guess for each RPM
    iterations = 0;

    while (abs(net_power.magnitude) > power_tol and iterations < 1000):
        #1: Step 1: Pick a thermo-state and RPM.
        compressor_inlet_state = ThermoState(P=P_amb, T=T_amb, mdot=massflow_guess)

        #2: Go to the compressor map, and with the thermo state get an outlet state and power.
        combustor_inlet_state, compressor_power = compressor.calculate(compressor_inlet_state, rpm=rpm)

        #3: Go to combustor, and with thermo state get an outlet state and power:
        turbine_inlet_state = combustor.calculate(combustor_inlet_state)

        #4: Go to turbine:
        turbine_exit_state, turbine_power = turbine.calculate(turbine_inlet_state, rpm=rpm)

        # Check the power balance, and print the compressor pressure ratio:
        pr_compressor = combustor_inlet_state.P_total / compressor_inlet_state.P_total
        # print(f"Compressor Pressure Ratio: {pr_compressor:.3f}")
        # print(f"Compressor Power Demand: {compressor_power:.2f}")
        # print(f"Turbine Power Generated: {turbine_power:.2f}")
        net_power = turbine_power - compressor_power

        # print(f"Net Power (Turbine - Compressor): {net_power:.2f}")

        if (abs(net_power.magnitude) < power_tol):
            print("Power Balanced!")
        else:
            if (net_power.magnitude > 0):
                # print("Excess Power from Turbine. Increase Massflow Guess.")
                massflow_guess *= 1.01;
            else:
                # print("Insufficient Power from Turbine. Decrease Massflow Guess.")
                massflow_guess *=  0.99;
        # print(f"New Massflow Guess: {massflow_guess:.4f}")

        # time.sleep(0.1)
        iterations += 1;
        # print(iterations);
    
    print(f"Converged in {iterations} iterations.")
    print(f"Final Massflow: {massflow_guess:.4f}")
    # Save converged results
    results['rpm'].append(rpm.magnitude)
    results['mass_flow'].append(massflow_guess.to(ureg.kg / ureg.s).magnitude)
    results['pressure_ratio'].append(pr_compressor.magnitude)
    results['compressor_power'].append(compressor_power.to(ureg.watt).magnitude)
    results['turbine_power'].append(turbine_power.to(ureg.watt).magnitude)
    results['net_power'].append(net_power.to(ureg.watt).magnitude)

# Plot results
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Mass flow vs RPM
axes[0, 0].plot(results['rpm'], results['mass_flow'], 'b-o', markersize=3)
axes[0, 0].set_xlabel('RPM')
axes[0, 0].set_ylabel('Mass Flow [kg/s]')
axes[0, 0].set_title('Mass Flow vs RPM')
axes[0, 0].grid(True, alpha=0.3)

# Pressure ratio vs RPM
axes[0, 1].plot(results['rpm'], results['pressure_ratio'], 'r-o', markersize=3)
axes[0, 1].set_xlabel('RPM')
axes[0, 1].set_ylabel('Pressure Ratio')
axes[0, 1].set_title('Compressor Pressure Ratio vs RPM')
axes[0, 1].grid(True, alpha=0.3)

# Powers vs RPM
axes[1, 0].plot(results['rpm'], [p/1000 for p in results['compressor_power']], 'b-o', markersize=3, label='Compressor')
axes[1, 0].plot(results['rpm'], [p/1000 for p in results['turbine_power']], 'r-o', markersize=3, label='Turbine')
axes[1, 0].set_xlabel('RPM')
axes[1, 0].set_ylabel('Power [kW]')
axes[1, 0].set_title('Component Power vs RPM')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# Net power vs RPM
axes[1, 1].plot(results['rpm'], [p/1000 for p in results['net_power']], 'g-o', markersize=3)
axes[1, 1].axhline(y=0, color='k', linestyle='--', alpha=0.5)
axes[1, 1].set_xlabel('RPM')
axes[1, 1].set_ylabel('Net Power [kW]')
axes[1, 1].set_title('Net Power (Turbine - Compressor) vs RPM')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()





