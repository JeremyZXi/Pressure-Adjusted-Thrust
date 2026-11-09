# Pressure-Adjusted Thrust
This is the unofficial repo for pressure-adjusted thrust sim 🐜🐜🐜🐜

Run this from the project folder(or just run that file:

```bash
python rocket.py
```

- **`rocket.py`** simulates the flight and saves `result/flight_data.csv`, which contains the thrust curve and other flight data. Change rocket inputs like dry mass here.
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

As for the rest of the files, we don't need to run them separately
