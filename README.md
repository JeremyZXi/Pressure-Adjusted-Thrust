# Pressure-Adjusted Thrust
This is the unofficial repo for pressure-adjusted thrust sim 🐜🐜🐜🐜

Run this from the project folder(or just run that file:

```bash
python rocket.py
```

- **`rocket.py`** simulates the flight and saves `result/flight_data.csv`. Change rocket inputs like dry mass here.
- **`visualization.ipynb`** compares the flight runs in `result/`. Open it from the project folder.
- **`engine_model_py/config.py`** holds engine inputs like the starting fuel and nitrous masses.

`flight_data.csv` contains

| Columns | Content                                                                                                   |
|---|-----------------------------------------------------------------------------------------------------------|
| `time`, `altitude`, `velocity`, `acceleration` | Flight motion in s, m, m/s, and m/s²                                                                      |
| `mass`, `CG` | Total rocket mass in kg and center-of-gravity position in m                                               |
| `mach`, `q`, `density` | Mach number, dynamic pressure in Pa, and air density in kg/m³                                             |
| `CD`, `CN`, `CP` | Drag coefficient, normal-force coefficient, and center-of-pressure position from the aero table           |
| `thrust` | Thrust in N, including both tank and ambient pressure effects                                             |
| `tankPressureRatio` | Current nitrous tank pressure divided by its initial pressure. `1.0` means 100% of the starting pressure. |

As for the rest of the files, we don't need to run them separately:

| File | What it does |
| --- | --- |
| `engine_history.py` | Works out how pressure and propellant mass change during the burn. |
| `ox_tank.py` | Calculates the pressure in the nitrous tank. |
| `fuel_tank.py` | Calculates the pressure that pushes fuel out of its tank. |
| `injector.py` | Calculates how much fuel and nitrous flow into the engine. |
| `nozzle.py` | Calculates how much thrust the nozzle produces for a given chamber and air pressure. |
| `integrator.py` | Advances the engine model one time step at a time. |
| `n2o_props.py` and `n2o_saturation.csv` | Look up nitrous properties, such as pressure and density. |
| `Eureka3_Machnum_RASAero2_sim.CSV` | Gives the simulator the drag coefficient at different speeds. |
