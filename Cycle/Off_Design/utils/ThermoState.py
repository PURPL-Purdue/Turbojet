"""
ThermoState class for managing thermodynamic state properties.
Based on the architecture.md specification.
"""

import pint
from typing import Optional, Literal
from utils.units import ureg


class ThermoState:
    """
    Represents a thermodynamic state with total pressure, total temperature,
    mass flow, and gas properties.
    """

    def __init__(
        self,
        P: pint.Quantity,
        T: pint.Quantity,
        mdot: pint.Quantity,
        composition: Literal["air", "combustion_products"] = "air",
        gamma: Optional[float] = None,
        cp: Optional[pint.Quantity] = None
    ):
        """
        Initialize a thermodynamic state.

        Parameters
        ----------
        P : pint.Quantity
            Total pressure (Pa)
        T : pint.Quantity
            Total temperature (K)
        mdot : pint.Quantity
            Mass flow rate (kg/s)
        composition : str, optional
            Gas composition: "air" or "combustion_products", default "air"
        gamma : float, optional
            Specific heat ratio. If None, uses default based on composition.
        cp : pint.Quantity, optional
            Specific heat at constant pressure (J/kg-K). If None, uses default based on composition.
        """
        self.P_total = P.to(ureg.pascal)
        self.T_total = T.to(ureg.kelvin)
        self.mass_flow = mdot.to(ureg.kg / ureg.s)
        self.composition = composition

        # Set gas properties based on composition if not provided
        if composition == "air":
            self.gamma = gamma if gamma is not None else 1.4
            self.cp = cp if cp is not None else 1005 * ureg.joule / (ureg.kg * ureg.kelvin)
        elif composition == "combustion_products":
            self.gamma = gamma if gamma is not None else 1.33
            self.cp = cp if cp is not None else 1148 * ureg.joule / (ureg.kg * ureg.kelvin)
        else:
            raise ValueError(f"Invalid composition: {composition}. Must be 'air' or 'combustion_products'")

    def to_corrected(self, P_ref: pint.Quantity = None, T_ref: pint.Quantity = None):
        """
        Convert to corrected (non-dimensional) parameters.

        Parameters
        ----------
        P_ref : pint.Quantity, optional
            Reference pressure. Default is standard sea level (101325 Pa)
        T_ref : pint.Quantity, optional
            Reference temperature. Default is standard sea level (288.15 K)

        Returns
        -------
        dict
            Dictionary with corrected mass flow, pressure, and temperature
        """
        if P_ref is None:
            P_ref = 101325 * ureg.pascal
        if T_ref is None:
            T_ref = 288.15 * ureg.kelvin

        # Corrected mass flow: mdot * sqrt(T/T_ref) / (P/P_ref)
        mdot_corrected = self.mass_flow * (self.T_total / T_ref)**0.5 / (self.P_total / P_ref)

        # Corrected pressure ratio
        P_corrected = self.P_total / P_ref

        # Corrected temperature ratio
        T_corrected = self.T_total / T_ref

        return {
            'mdot_corrected': mdot_corrected.to_base_units(),
            'P_corrected': P_corrected.to_base_units(),
            'T_corrected': T_corrected.to_base_units()
        }

    @classmethod
    def from_corrected(
        cls,
        mdot_corrected: pint.Quantity,
        P_corrected: float,
        T_corrected: float,
        P_ref: pint.Quantity = None,
        T_ref: pint.Quantity = None,
        composition: Literal["air", "combustion_products"] = "air"
    ):
        """
        Create a ThermoState from corrected parameters.

        Parameters
        ----------
        mdot_corrected : pint.Quantity
            Corrected mass flow
        P_corrected : float
            Corrected pressure ratio (P/P_ref)
        T_corrected : float
            Corrected temperature ratio (T/T_ref)
        P_ref : pint.Quantity, optional
            Reference pressure. Default is standard sea level (101325 Pa)
        T_ref : pint.Quantity, optional
            Reference temperature. Default is standard sea level (288.15 K)
        composition : str, optional
            Gas composition

        Returns
        -------
        ThermoState
            New ThermoState object
        """
        if P_ref is None:
            P_ref = 101325 * ureg.pascal
        if T_ref is None:
            T_ref = 288.15 * ureg.kelvin

        # Convert back to actual values
        P = P_corrected * P_ref
        T = T_corrected * T_ref
        mdot = mdot_corrected * (P / P_ref) / (T / T_ref)**0.5

        return cls(P=P, T=T, mdot=mdot, composition=composition)

    def copy(self):
        """
        Create a deep copy of this ThermoState.

        Returns
        -------
        ThermoState
            New ThermoState object with same properties
        """
        return ThermoState(
            P=self.P_total,
            T=self.T_total,
            mdot=self.mass_flow,
            composition=self.composition,
            gamma=self.gamma,
            cp=self.cp
        )

    def __repr__(self):
        return (f"ThermoState(P={self.P_total:.2f}, T={self.T_total:.2f}, "
                f"mdot={self.mass_flow:.4f}, composition='{self.composition}')")

    def __str__(self):
        return (f"Thermodynamic State:\n"
                f"  P_total: {self.P_total:.2f}\n"
                f"  T_total: {self.T_total:.2f}\n"
                f"  Mass flow: {self.mass_flow:.4f}\n"
                f"  Composition: {self.composition}\n"
                f"  γ: {self.gamma}\n"
                f"  cp: {self.cp:.2f}")
