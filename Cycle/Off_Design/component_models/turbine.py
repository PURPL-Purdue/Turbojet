import pint
from typing import Tuple
from utils.ThermoState import ThermoState
from utils.units import ureg


class Turbine:

    def __init__(
        self,
        design_efficiency: float = 0.85,
        back_pressure: pint.Quantity = 2 * 101325 * ureg.pascal,
        design_rpm: float = 80000,
        efficiency_falloff: float = 0.5,
    ):
        self.design_efficiency = design_efficiency
        self.back_pressure = back_pressure
        self.design_rpm = design_rpm
        self.efficiency_falloff = efficiency_falloff  # quadratic penalty coefficient

    def get_efficiency(self, rpm: float) -> float:

        # Random slop code to scale the efficiency based on RPM deviation from design
        delta_n = (rpm - self.design_rpm) / self.design_rpm
        efficiency = self.design_efficiency * (1 - self.efficiency_falloff * delta_n**2)
        # return 1;
        return max(efficiency, 0.1)  # floor to prevent negative/zero efficiency

    def calculate(self, inlet_state: ThermoState, rpm: pint.Quantity) -> Tuple[ThermoState, pint.Quantity]:
        rpm_val = rpm.to(ureg.rpm).magnitude
        efficiency = self.get_efficiency(rpm_val)

        gamma = inlet_state.gamma
        T_in = inlet_state.T_total.magnitude
        P_in = inlet_state.P_total.magnitude
        mdot_val = inlet_state.mass_flow.to(ureg.kg / ureg.s).magnitude
        cp = inlet_state.cp.to(ureg.joule / (ureg.kg * ureg.kelvin)).magnitude

        # Outlet pressure (fixed back-pressure)
        P_out = self.back_pressure.to(ureg.pascal).magnitude
        expansion_ratio = P_in / P_out

        # Isentropic temperature drop
        T_out_isentropic = T_in * (1 / expansion_ratio) ** ((gamma - 1) / gamma)

        # Actual temperature drop (accounting for efficiency)
        # For turbine: T_out = T_in - eta * (T_in - T_out_isentropic)
        T_out = T_in - efficiency * (T_in - T_out_isentropic)

        outlet_state = ThermoState(
            P=P_out * ureg.pascal,
            T=T_out * ureg.kelvin,
            mdot=inlet_state.mass_flow,
            composition=inlet_state.composition,
            gamma=inlet_state.gamma,
            cp=inlet_state.cp
        )

        # Power generated (positive for turbine)
        power = mdot_val * cp * (T_in - T_out) * ureg.watt

        return outlet_state, power