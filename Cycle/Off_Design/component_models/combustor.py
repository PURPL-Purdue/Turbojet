"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Combustion Startup
Synopsis:
Author: Yuuki 
Input:
Output:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""

from utils.ThermoState import ThermoState
from utils.units import ureg
import numpy as np

class Combustor:

    # Combustor pressure drop:
    combustor_pressure_drop = 0.05  # 5% pressure drop

    def __init__(
        self,
    ):
        pass

    def calculate(self, inlet_state: ThermoState) -> ThermoState:
        T_out = combustorTemp(inlet_state.P_total.to(ureg.atm).magnitude,
                              inlet_state.mass_flow.magnitude,
                              inlet_state.T_total.to(ureg.kelvin).magnitude,
                              targetPhi=1.2) * ureg.kelvin
        outlet_state = ThermoState.copy(inlet_state);
        outlet_state.P_total = inlet_state.P_total * (1-self.combustor_pressure_drop);  # assuming no pressure drop
        outlet_state.T_total = T_out;
        return outlet_state;



def combustorTemp(pressure, massFlow, inlet_temp, targetPhi):
    """
    Combustion Temperature 3D Interpolation Function
    
    Input: 
        pressure (atm)
        massFlow (kg/s)
        inlet_temp (K)
        targetPhi - target equivalence ratio
    Output: 
        temperature (K) using 3D trilinear interpolation
    """
    
    # Constants
    FA_st = 0.064  # stoichiometric fuel-air ratio
    fuelFlow = 0.007  # fuel flow in [kg/s]
    
    # Calculate phi from mass flow
    phi = (fuelFlow / massFlow) / FA_st
    
    # Define dimension arrays
    phi_vals = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])
    press_vals = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
    temp_vals = np.array([270, 300, 350, 400, 450])
    temp_data = np.zeros((9, 5, 5))  # temp_data org = [phi, press, temperature]
    
    # Temperature = 270K
    temp_data[:, :, 0] = np.array([
        [2261.22, 2271.50, 2278.44, 2283.62, 2287.72],
        [2259.54, 2265.83, 2269.72, 2272.43, 2274.46],
        [2195.86, 2198.39, 2199.91, 2200.95, 2201.72],
        [2120.59, 2121.89, 2122.67, 2123.20, 2123.60],
        [2045.30, 2046.06, 2046.51, 2046.82, 2047.05],
        [1971.95, 1972.41, 1972.69, 1972.88, 1973.02],
        [1900.85, 1901.13, 1901.31, 1901.42, 1901.51],
        [1831.94, 1832.12, 1832.23, 1832.30, 1832.35],
        [1765.11, 1765.22, 1765.28, 1765.32, 1765.36]
    ])
    
    # Temperature = 300K
    temp_data[:, :, 1] = np.array([
        [2273.69, 2284.38, 2291.60, 2297.00, 2301.27],
        [2274.08, 2280.97, 2285.26, 2288.27, 2290.53],
        [2212.67, 2215.52, 2217.24, 2218.42, 2219.29],
        [2138.05, 2139.52, 2140.40, 2141.01, 2141.45],
        [2062.91, 2063.77, 2064.28, 2064.63, 2064.90],
        [1989.53, 1990.06, 1990.38, 1990.60, 1990.76],
        [1918.34, 1918.67, 1918.86, 1919.00, 1919.10],
        [1849.30, 1849.51, 1849.63, 1849.71, 1849.77],
        [1782.32, 1782.44, 1782.52, 1782.57, 1782.61]
    ])
    
    # Temperature = 350K
    temp_data[:, :, 2] = np.array([
        [2294.13, 2305.48, 2313.18, 2318.94, 2323.50],
        [2297.61, 2305.53, 2310.53, 2314.07, 2316.74],
        [2240.25, 2243.71, 2245.81, 2247.25, 2248.32],
        [2166.94, 2168.73, 2169.80, 2170.54, 2171.09],
        [2092.14, 2093.19, 2093.82, 2094.26, 2094.58],
        [2018.78, 2019.43, 2019.82, 2020.09, 2020.29],
        [1947.45, 1947.86, 1948.11, 1948.28, 1948.40],
        [1878.22, 1878.48, 1878.63, 1878.74, 1878.82],
        [1811.00, 1811.16, 1811.26, 1811.32, 1811.37]
    ])
    
    # Temperature = 400K
    temp_data[:, :, 3] = np.array([
        [2314.15, 2326.19, 2334.36, 2340.48, 2345.35],
        [2320.34, 2329.31, 2335.04, 2339.13, 2342.23],
        [2267.28, 2271.44, 2273.97, 2275.72, 2277.02],
        [2195.56, 2197.73, 2199.04, 2199.94, 2200.61],
        [2121.24, 2122.52, 2123.29, 2123.82, 2124.22],
        [2047.97, 2048.77, 2049.25, 2049.58, 2049.82],
        [1976.56, 1977.07, 1977.37, 1977.58, 1977.73],
        [1907.16, 1907.48, 1907.67, 1907.80, 1907.90],
        [1839.72, 1839.92, 1840.04, 1840.13, 1840.19]
    ])
    
    # Temperature = 450K
    temp_data[:, :, 4] = np.array([
        [2332.60, 2345.27, 2353.89, 2360.36, 2365.51],
        [2340.97, 2350.93, 2357.35, 2361.97, 2365.50],
        [2292.09, 2297.00, 2299.99, 2302.07, 2303.61],
        [2222.15, 2224.73, 2226.29, 2227.37, 2228.17],
        [2148.42, 2149.96, 2150.88, 2151.52, 2151.99],
        [2075.30, 2076.27, 2076.85, 2077.24, 2077.54],
        [2003.86, 2004.48, 2004.85, 2005.10, 2005.29],
        [1934.34, 1934.73, 1934.97, 1935.13, 1935.25],
        [1866.72, 1866.97, 1867.11, 1867.22, 1867.29]
    ])
    
    if phi < targetPhi:
        # print("Not enough mass flow for target phi... passing compressor temperature.")
        T = inlet_temp
        return T
    
    # Find surrounding indices for phi
    if phi <= phi_vals[0]:
        i1 = 0
        i2 = 0
        w_phi = 0
    elif phi >= phi_vals[-1]:
        i1 = len(phi_vals) - 1
        i2 = len(phi_vals) - 1
        w_phi = 0
    else:
        i1 = np.where(phi_vals <= phi)[0][-1]
        i2 = i1 + 1
        w_phi = (phi - phi_vals[i1]) / (phi_vals[i2] - phi_vals[i1])
    
    # Find surrounding indices for pressure
    if pressure <= press_vals[0]:
        j1 = 0
        j2 = 0
        w_press = 0
    elif pressure >= press_vals[-1]:
        j1 = len(press_vals) - 1
        j2 = len(press_vals) - 1
        w_press = 0
    else:
        j1 = np.where(press_vals <= pressure)[0][-1]
        j2 = j1 + 1
        w_press = (pressure - press_vals[j1]) / (press_vals[j2] - press_vals[j1])
    
    # Find surrounding indices for inlet temperature
    if inlet_temp <= temp_vals[0]:
        k1 = 0
        k2 = 0
        w_temp = 0
    elif inlet_temp >= temp_vals[-1]:
        k1 = len(temp_vals) - 1
        k2 = len(temp_vals) - 1
        w_temp = 0
    else:
        k1 = np.where(temp_vals <= inlet_temp)[0][-1]
        k2 = k1 + 1
        w_temp = (inlet_temp - temp_vals[k1]) / (temp_vals[k2] - temp_vals[k1])
    
    # Trilinear interpolation (8 corner points of the cube)
    T111 = temp_data[i1, j1, k1]
    T112 = temp_data[i1, j1, k2]
    T121 = temp_data[i1, j2, k1]
    T122 = temp_data[i1, j2, k2]
    T211 = temp_data[i2, j1, k1]
    T212 = temp_data[i2, j1, k2]
    T221 = temp_data[i2, j2, k1]
    T222 = temp_data[i2, j2, k2]
    
    # Interpolate along temperature dimension first
    T11 = T111 * (1 - w_temp) + T112 * w_temp
    T12 = T121 * (1 - w_temp) + T122 * w_temp
    T21 = T211 * (1 - w_temp) + T212 * w_temp
    T22 = T221 * (1 - w_temp) + T222 * w_temp
    
    # Interpolate along pressure dimension
    T1 = T11 * (1 - w_press) + T12 * w_press
    T2 = T21 * (1 - w_press) + T22 * w_press
    
    # Interpolate along phi dimension
    T = T1 * (1 - w_phi) + T2 * w_phi
    
    return T