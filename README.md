# Pressure-Adjusted Thrust
This is the unofficial repo for pressure-adjusted thrust sim 🐜🐜🐜🐜

**Run `engine_model_py/rocket.py`.** It simulates the flight and saves the results as `engine_model_py/flight_data_python_engine.csv`.

From the project folder, you can also run:

```bash
python -m engine_model_py.rocket
```

- **`config.py`** — Change engine inputs here, like how much fuel and nitrous you start with.
- **`rocket.py`** — Change rocket inputs here, like dry mass. This is the file to run.
- **`visualization.ipynb`** — Open this to compare the runs in `result/`. Copy your new flight CSV there before rerunning the notebook.

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
