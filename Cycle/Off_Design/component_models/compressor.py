"""
Compressor performance model using map interpolation.
"""

import pint
from typing import Tuple
import pandas as pd
import numpy as np
from utils.ThermoState import ThermoState

ureg = pint.UnitRegistry()


class CompressorMap:
    """
    Compressor performance model using experimental map data.

    The map provides pressure ratio and efficiency as a function of
    corrected mass flow and corrected RPM.
    """

    def __init__(self, map_file: str):
        """
        Initialize the compressor map.

        Parameters
        ----------
        map_file : str
            Path to the compressor map CSV file
        """
        self.map_file = map_file
        self.map_data = None

        # TODO: Load the map data from CSV
        # Expected columns: corrected_rpm, corrected_mass_flow, pressure_ratio, efficiency
        # self.map_data = pd.read_csv(map_file)

        # TODO: Set up interpolation functions
        # - Pressure ratio as f(corrected_rpm, corrected_mdot)
        # - Efficiency as f(corrected_rpm, corrected_mdot)

    def calculate(
        self,
        inlet_state: ThermoState,
        rpm: pint.Quantity
    ) -> Tuple[ThermoState, pint.Quantity]:
        """
        Calculate compressor outlet state and power consumption.

        Parameters
        ----------
        inlet_state : ThermoState
            Compressor inlet thermodynamic state
        rpm : pint.Quantity
            Shaft rotational speed

        Returns
        -------
        outlet_state : ThermoState
            Compressor outlet thermodynamic state
        power : pint.Quantity
            Compressor power consumption (positive value)
        """

        # TODO: Convert inlet state to corrected parameters
        # corrected = inlet_state.to_corrected()
        # corrected_rpm = rpm * sqrt(T_ref / inlet_state.T_total)

        # TODO: Interpolate map to get pressure ratio and efficiency
        # PR = self._interpolate_pr(corrected_rpm, corrected['mdot_corrected'])
        # eta = self._interpolate_efficiency(corrected_rpm, corrected['mdot_corrected'])

        # TODO: Calculate outlet pressure and temperature
        # P_out = inlet_state.P_total * PR
        # T_out = inlet_state.T_total * (1 + (PR**((gamma-1)/gamma) - 1) / eta)

        # TODO: Create outlet state
        # outlet_state = ThermoState(
        #     P=P_out,
        #     T=T_out,
        #     mdot=inlet_state.mass_flow,
        #     composition=inlet_state.composition
        # )

        # TODO: Calculate power
        # power = inlet_state.mass_flow * inlet_state.cp * (T_out - inlet_state.T_total)

        # Placeholder return
        outlet_state = inlet_state.copy()
        power = 0 * ureg.watt

        return outlet_state, power

    def _interpolate_pr(self, corrected_rpm: float, corrected_mdot: pint.Quantity) -> float:
        """
        Interpolate pressure ratio from map.

        TODO: Implement 2D interpolation
        """
        return 1.5  # Placeholder

    def _interpolate_efficiency(self, corrected_rpm: float, corrected_mdot: pint.Quantity) -> float:
        """
        Interpolate isentropic efficiency from map.

        TODO: Implement 2D interpolation
        """
        return 0.75  # Placeholder
    
    @staticmethod
    def plot_from_file(map_file_path: str):
        """
        Plot compressor map from a tab-separated text file.

        Parameters
        ----------
        map_file_path : str
            Path to the tab-separated text file with columns:
            - RPM (in thousands)
            - Pressure ratio
            - Efficiency
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        # Load map data from tab-separated file
        # No header row, so we assign column names manually
        map_data = pd.read_csv(
            map_file_path,
            sep='\t',
            header=None,
            names=['rpm_thousands', 'pressure_ratio', 'efficiency']
        )

        # Convert RPM from thousands to actual RPM
        map_data['rpm'] = map_data['rpm_thousands'] * 1000

        # Create figure with subplots
        fig = plt.figure(figsize=(14, 5))

        # Plot 1: Pressure ratio vs RPM, colored by efficiency
        ax1 = fig.add_subplot(131)
        scatter1 = ax1.scatter(
            map_data['rpm'],
            map_data['pressure_ratio'],
            c=map_data['efficiency'],
            cmap='viridis',
            marker='o'
        )
        ax1.set_xlabel('RPM')
        ax1.set_ylabel('Pressure Ratio')
        ax1.set_title('Pressure Ratio vs RPM')
        plt.colorbar(scatter1, ax=ax1, label='Efficiency')

        # Plot 2: Efficiency vs RPM
        ax2 = fig.add_subplot(132)
        ax2.scatter(
            map_data['rpm'],
            map_data['efficiency'],
            c='b',
            marker='o'
        )
        ax2.set_xlabel('RPM')
        ax2.set_ylabel('Efficiency')
        ax2.set_title('Efficiency vs RPM')

        # Plot 3: 3D plot
        ax3 = fig.add_subplot(133, projection='3d')
        scatter3 = ax3.scatter(
            map_data['rpm'],
            map_data['pressure_ratio'],
            map_data['efficiency'],
            c=map_data['efficiency'],
            cmap='viridis',
            marker='o'
        )
        ax3.set_xlabel('RPM')
        ax3.set_ylabel('Pressure Ratio')
        ax3.set_zlabel('Efficiency')
        ax3.set_title('3D Map View')

        plt.tight_layout()
        plt.show()

        return map_data
