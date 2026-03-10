import pint
from typing import Tuple, Optional
import pandas as pd
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator
from utils.ThermoState import ThermoState
from utils.units import ureg


class Compressor:

    def __init__(
        self,
        chics_file: str,
        isoeff_file: str,
        design_pr: Optional[float] = None,
        design_efficiency: Optional[float] = None,
        design_rpm: Optional[float] = None,
        design_mass_flow: Optional[float] = None,
    ):
        """
        Initialize compressor with map data files and optional scaling.

        Args:
            chics_file: Path to characteristics CSV (columns: RPM, mass_flow, PR)
            isoeff_file: Path to iso-efficiency CSV (columns: efficiency, mass_flow, PR)
            design_pr: Target design pressure ratio (optional)
            design_efficiency: Target design/max efficiency as decimal (optional)
            design_rpm: Target design RPM [rev/min] (optional)
            design_mass_flow: Target design mass flow [kg/s] (optional)
        """

        # Load raw data
        self._raw_chics = self._load_chics(chics_file)
        self._raw_isoeff = self._load_isoeff(isoeff_file)

        # Find reference design point from raw data
        self.ref_design_point = self._find_reference_design_point()

        # Store target design point (use reference values if not specified)
        self.design_point = {
            'pr': design_pr if design_pr is not None else self.ref_design_point['pr'],
            'efficiency': design_efficiency if design_efficiency is not None else self.ref_design_point['efficiency'],
            'rpm': design_rpm if design_rpm is not None else self.ref_design_point['rpm'],
            'mass_flow': design_mass_flow if design_mass_flow is not None else self.ref_design_point['mass_flow'],
        }

        # Calculate scaling factors
        self.scale_factors = self._calculate_scale_factors()

        # Apply scaling to get working data
        self.chics_data = self._scale_chics(self._raw_chics.copy())
        self.isoeff_data = self._scale_isoeff(self._raw_isoeff.copy())

        # Create interpolators from scaled data
        self._build_interpolators()

    def _load_chics(self, filepath: str) -> pd.DataFrame:
        """Load characteristics (speed lines) data."""
        df = pd.read_csv(filepath)
        df.columns = ['rpm', 'mass_flow', 'pressure_ratio']
        return df

    def _load_isoeff(self, filepath: str) -> pd.DataFrame:
        """Load iso-efficiency contours data."""
        df = pd.read_csv(filepath)
        df.columns = ['efficiency', 'mass_flow', 'pressure_ratio']
        # Convert efficiency from percentage to decimal
        df['efficiency'] = df['efficiency'] / 100.0
        return df

    def _find_reference_design_point(self) -> dict:
        """
        Find the reference design point from raw map data.

        Uses the point of maximum efficiency as the reference design point.
        """
        # Find max efficiency in the isoeff data
        max_eff = self._raw_isoeff['efficiency'].max()

        # Get the mass_flow and PR at max efficiency (average if multiple points)
        max_eff_data = self._raw_isoeff[self._raw_isoeff['efficiency'] == max_eff]
        ref_mass_flow = max_eff_data['mass_flow'].mean()
        ref_pr = max_eff_data['pressure_ratio'].mean()

        # Find the RPM that corresponds to this operating point
        # Use the CHICS data to find the speed line closest to this (mass_flow, PR)
        # Simple approach: find the point in CHICS closest to (ref_mass_flow, ref_pr)
        distances = np.sqrt(
            (self._raw_chics['mass_flow'] - ref_mass_flow)**2 +
            (self._raw_chics['pressure_ratio'] - ref_pr)**2
        )
        closest_idx = distances.idxmin()
        ref_rpm = self._raw_chics.loc[closest_idx, 'rpm']

        return {
            'efficiency': max_eff,
            'mass_flow': ref_mass_flow,
            'pr': ref_pr,
            'rpm': ref_rpm,
        }

    def _calculate_scale_factors(self) -> dict:
        """Calculate scaling factors to transform raw map to target design point."""
        ref = self.ref_design_point
        des = self.design_point

        return {
            'rpm': des['rpm'] / ref['rpm'],
            'mass_flow': des['mass_flow'] / ref['mass_flow'],
            # PR scaling preserves PR=1 as the baseline (no compression)
            'pr': (des['pr'] - 1) / (ref['pr'] - 1),
            # Efficiency scaling: shift to match target max efficiency
            'efficiency_delta': des['efficiency'] - ref['efficiency'],
        }

    def _scale_chics(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply scaling to characteristics data."""
        sf = self.scale_factors
        df['rpm'] = df['rpm'] * sf['rpm']
        df['mass_flow'] = df['mass_flow'] * sf['mass_flow']
        df['pressure_ratio'] = 1 + (df['pressure_ratio'] - 1) * sf['pr']
        return df

    def _scale_isoeff(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply scaling to iso-efficiency data."""
        sf = self.scale_factors
        df['mass_flow'] = df['mass_flow'] * sf['mass_flow']
        df['pressure_ratio'] = 1 + (df['pressure_ratio'] - 1) * sf['pr']
        df['efficiency'] = df['efficiency'] + sf['efficiency_delta']
        # Clamp efficiency to reasonable bounds
        df['efficiency'] = df['efficiency'].clip(0.0, 1.0)
        return df

    def _build_interpolators(self):
        """Build 2D interpolators from the scaled map data."""
        # PR interpolator: (rpm, mass_flow) -> PR
        self.pr_interpolator = CloughTocher2DInterpolator(
            points=np.column_stack([
                self.chics_data['rpm'].values,
                self.chics_data['mass_flow'].values
            ]),
            values=self.chics_data['pressure_ratio'].values
        )

        # Efficiency interpolator: (mass_flow, PR) -> efficiency
        self.eff_interpolator = CloughTocher2DInterpolator(
            points=np.column_stack([
                self.isoeff_data['mass_flow'].values,
                self.isoeff_data['pressure_ratio'].values
            ]),
            values=self.isoeff_data['efficiency'].values
        )

    def query(self, rpm: float, mass_flow: float) -> Tuple[float, float]:
        """
        Query the compressor map for pressure ratio and efficiency.

        Args:
            rpm: Shaft speed [rev/min]
            mass_flow: Mass flow rate [kg/s]

        Returns:
            (pressure_ratio, efficiency) tuple
        """
        pr = float(self.pr_interpolator(rpm, mass_flow))

        # If PR is off-map, find nearest point in chics data
        if np.isnan(pr):
            distances = np.sqrt(
                ((self.chics_data['rpm'] - rpm) / self.chics_data['rpm'].std())**2 +
                ((self.chics_data['mass_flow'] - mass_flow) / self.chics_data['mass_flow'].std())**2
            )
            nearest_idx = distances.idxmin()
            pr = float(self.chics_data.loc[nearest_idx, 'pressure_ratio'])

        eff = float(self.eff_interpolator(mass_flow, pr))
        # eff = 0.7;

        # If efficiency is off-map, find nearest point in isoeff data
        if np.isnan(eff):
            distances = np.sqrt(
                ((self.isoeff_data['mass_flow'] - mass_flow) / self.isoeff_data['mass_flow'].std())**2 +
                ((self.isoeff_data['pressure_ratio'] - pr) / self.isoeff_data['pressure_ratio'].std())**2
            )
            nearest_idx = distances.idxmin()
            eff = float(self.isoeff_data.loc[nearest_idx, 'efficiency'])

        return pr, eff

    def calculate(self, inlet_state: ThermoState, rpm: pint.Quantity) -> Tuple[ThermoState, pint.Quantity]:
        """
        Calculate compressor outlet state and power consumption.

        Args:
            inlet_state: Inlet thermodynamic state
            rpm: Shaft speed with units

        Returns:
            (outlet_state, power) tuple
        """
        rpm_val = rpm.to(ureg.rpm).magnitude
        mdot_val = inlet_state.mass_flow.to(ureg.kg / ureg.s).magnitude

        pr, eff = self.query(rpm_val, mdot_val)

        # # TEMPORARY CLAMP EFFICIENCY/PR:
        eff = 0.70;
        pr = 3.0;

        gamma = inlet_state.gamma
        T_in = inlet_state.T_total.magnitude
        P_in = inlet_state.P_total.magnitude

        # Isentropic temperature rise
        T_out_isentropic = T_in * (pr ** ((gamma - 1) / gamma))

        # Actual temperature rise (accounting for efficiency)
        T_out = T_in + (T_out_isentropic - T_in) / eff

        P_out = P_in * pr

        outlet_state = ThermoState(
            P=P_out * ureg.pascal,
            T=T_out * ureg.kelvin,
            mdot=inlet_state.mass_flow,
            composition=inlet_state.composition,
            gamma=inlet_state.gamma,
            cp=inlet_state.cp
        )

        cp = inlet_state.cp.to(ureg.joule / (ureg.kg * ureg.kelvin)).magnitude
        power = mdot_val * cp * (T_out - T_in) * ureg.watt

        return outlet_state, power

    def plot_map(self):
        """
        Plot the compressor map as a colormap.

        Shows pressure ratio vs mass flow with RPM as speed lines
        and efficiency as color.
        """
        import matplotlib.pyplot as plt

        rpm_min, rpm_max = self.chics_data['rpm'].min(), self.chics_data['rpm'].max()
        mdot_min, mdot_max = self.chics_data['mass_flow'].min(), self.chics_data['mass_flow'].max()

        rpm_grid = np.linspace(rpm_min, rpm_max, 50)
        mdot_grid = np.linspace(mdot_min, mdot_max, 100)
        RPM, MDOT = np.meshgrid(rpm_grid, mdot_grid)

        PR = self.pr_interpolator(RPM, MDOT)
        EFF = self.eff_interpolator(MDOT, PR)

        fig, ax = plt.subplots(figsize=(12, 8))

        mask = ~np.isnan(PR) & ~np.isnan(EFF)

        scatter = ax.scatter(
            MDOT[mask], PR[mask],
            c=EFF[mask],
            cmap='viridis',
            s=5,
            alpha=0.7
        )

        plt.colorbar(scatter, ax=ax, label='Efficiency')

        # Plot speed lines
        for rpm in self.chics_data['rpm'].unique():
            rpm_data = self.chics_data[self.chics_data['rpm'] == rpm].sort_values('mass_flow')
            ax.plot(rpm_data['mass_flow'], rpm_data['pressure_ratio'],
                   'k-', linewidth=1.5, alpha=0.8)
            ax.annotate(f'{int(rpm/1000)}k',
                       xy=(rpm_data['mass_flow'].iloc[-1], rpm_data['pressure_ratio'].iloc[-1]),
                       fontsize=8, alpha=0.7)

        # Mark design point
        ax.scatter(
            [self.design_point['mass_flow']],
            [self.design_point['pr']],
            c='red', s=100, marker='*', zorder=5,
            label=f"Design: PR={self.design_point['pr']:.2f}, η={self.design_point['efficiency']:.1%}"
        )
        ax.legend()

        ax.set_xlabel('Mass Flow [kg/s]')
        ax.set_ylabel('Pressure Ratio')
        ax.set_title('Compressor Map (Scaled)')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

        return fig, ax

    def print_scaling_info(self):
        """Print information about the map scaling."""
        print("Reference Design Point (from raw map):")
        print(f"  RPM:        {self.ref_design_point['rpm']:.0f}")
        print(f"  Mass Flow:  {self.ref_design_point['mass_flow']:.4f} kg/s")
        print(f"  PR:         {self.ref_design_point['pr']:.3f}")
        print(f"  Efficiency: {self.ref_design_point['efficiency']:.1%}")
        print()
        print("Target Design Point:")
        print(f"  RPM:        {self.design_point['rpm']:.0f}")
        print(f"  Mass Flow:  {self.design_point['mass_flow']:.4f} kg/s")
        print(f"  PR:         {self.design_point['pr']:.3f}")
        print(f"  Efficiency: {self.design_point['efficiency']:.1%}")
        print()
        print("Scale Factors:")
        print(f"  RPM:        {self.scale_factors['rpm']:.4f}")
        print(f"  Mass Flow:  {self.scale_factors['mass_flow']:.4f}")
        print(f"  PR:         {self.scale_factors['pr']:.4f}")
        print(f"  Eff Delta:  {self.scale_factors['efficiency_delta']:+.1%}")

    @staticmethod
    def plot_from_files(chics_file: str, isoeff_file: str):
        """Convenience method to load and plot map from files."""
        comp = Compressor(chics_file, isoeff_file)
        return comp.plot_map()
