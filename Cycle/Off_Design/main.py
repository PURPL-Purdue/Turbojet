# Main script?
import pint
from utils.ThermoState import ThermoState
from component_models.compressor import CompressorMap
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

# Right now let's just plot the existing map from references:
CompressorModel.plot_from_file("reference/etas_import_scaled_full.txt");
