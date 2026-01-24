# Microturbojet Off-Design Cycle Analysis Tool

**This document generated partially by AI**


## Analysis Approach

For each RPM setpoint:
1. Pick an RPM
2. Guess a mass flow rate
3. Use compressor map to determine efficiency and pressure ratio at corrected conditions
4. Calculate compressor power and outlet state
5. Pass through combustor (black box) → returns outlet state
6. Calculate turbine power output (expand to ambient)
7. Add rotordynamic losses (black box)
8. Check power balance; if not balanced, adjust mass flow and iterate
9. Once converged, increment RPM and repeat

**Key assumptions**: Steady-state, ground test (P_inlet = P_ambient), non-dimensionalized using corrected parameters

---

## File Structure

```
microturbojet_cycle/
├── main.py                   # Orchestration
├── config.py                 # Configuration parameters
├── models/
│   ├── compressor.py         # Compressor map interpolation
│   ├── turbine.py            # Turbine performance
│   ├── combustor.py          # Black box interface
│   └── rotordynamics.py      # Black box interface
├── core/
│   ├── thermodynamics.py     # Gas properties, state class
│   ├── solver.py             # Power balance iteration
│   └── nondim.py             # Corrected parameters
├── utils/
│   ├── interpolation.py      # Map utilities
│   └── plotting.py           # Visualization
└── data/
    └── compressor_map.csv    # Map data
```

---

## Core Data Structures

### ThermodynamicState
```
Properties:
- P_total (Pa)
- T_total (K)
- mass_flow (kg/s)
- composition (air vs combustion products)
- gamma
- cp (J/kg-K)

Methods:
- to_corrected()
- from_corrected()
- copy()
```

### OperatingPoint
```
Properties:
- RPM
- mass_flow
- compressor_states (inlet, outlet)
- combustor_states (inlet, outlet)
- turbine_states (inlet, outlet)
- compressor_PR
- compressor_efficiency
- turbine_efficiency
- compressor_power
- turbine_power
- mechanical_losses
- net_power (≈0 at convergence)
- fuel_flow
```

---