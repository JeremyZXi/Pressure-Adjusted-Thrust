# Pressure-adjusted engine model (N2O / IPA, piston-fed)

Everything needed is in this one folder: copy it anywhere and it runs.

A physics model of the Eureka Spaceshot CF engine that predicts how thrust, chamber
pressure, flow rates and O/F change during a burn as:

- **the N2O tank blows down.** The N2O self-pressurizes, cools as it boils, and loses
  pressure. Lower pressure means less flow, so the burn stretches out.
- **the rocket climbs.** The nozzle's pressure thrust grows as the air thins.

It replaces the old "80/20" linear tank-pressure assumption in the 1D framework with
an equilibrium tank model, injector flow, a chamber pressure balance and nozzle thrust.

> **Status (Nov 2026):** the code is complete and numerically verified. It has **not
> been validated against test data yet**: no full-duration N2O/IPA flow data exists yet.
> Several inputs are placeholders (see [Inputs](#inputs-to-set)). Treat absolute
> numbers (burn time, impulse, apogee) as estimates until those are filled in and
> the model is checked against a hotfire.

## Quick start

```bash
pip3 install numpy matplotlib        # CoolProp is optional (only for verify.py / regenerating the table)
cd engine_model_py

python3 check_cores.py               # self-test of the physics solvers (should say 21/21)
python3 static_fire.py               # sea-level static fire: prints a summary, plots, writes static_fire.csv
python3 static_fire.py --eng         # also writes thrust_curve.eng for OpenRocket
python3 flight.py                    # 1D vertical flight with altitude-adjusted thrust (~15 s to run)
python3 validate.py TEST.csv ...     # compare against measured data (see Validating)
python3 verify.py                    # numerical checks (mass balance, timestep, property table)
```

### What the default run gives (placeholder loads)

| Quantity | Ignition | End of burn |
|---|---|---|
| N2O tank pressure | 600 psi (11.4 °C) | 351 psi (−9.6 °C) |
| Chamber pressure | 345 psi | 240 psi |
| Sea-level thrust | 13.53 kN | 8.71 kN |
| O/F | 4.50 | 4.78 |

At ignition the model reproduces the RPA design point: Pc 345 psi (RPA 350), 6.60 kg/s,
O/F 4.50, 13.53 kN (RPA 13.56 kN). The burn lasts 10.5 s, with 117.6 kN·s total
impulse and an average sea-level Isp of 203 s.

In `flight.py`, the altitude correction raises burnout thrust from 8.7 to 9.9 kN. With
the placeholder vehicle mass, that's +2.3 km of apogee versus using sea-level thrust
the whole way.

## Inputs to set

All inputs live in `config.py` (`EngineConfig`), and `VehicleConfig` in `flight.py`.
Values tagged with a source come from the SEB wiki. **TBD values are placeholders. Get
real numbers from Propulsion.**

| Input | Current value | Source / status |
|---|---|---|
| Fuel feed | piston pushed by N2O gas | Eureka-3 Piston / IPA Tank / Downcomer pages. Confirm the flight vehicle uses this |
| N2O fill pressure | 600 psi (sets T₀ ≈ 11.4 °C) | CF Tank Testing, MEOP 600 psi |
| **N2O mass loaded** | 54 kg | **TBD** (placeholder = 10 s at design flow) |
| **N2O ullage at fill** | 10 % of tank volume | **TBD** (or set `ox_tank_volume` directly) |
| **IPA mass loaded** | 12 kg | **TBD** (placeholder = N2O / 4.5) |
| IPA density | 786 kg/m³ | IPA near 20 °C |
| **Piston seal friction** | 0 psi | **TBD** |
| **Injector C_d·A (N2O, IPA)** | calibrated: 99.2 / 22.8 mm² | **TBD.** Back-solved to hit the RPA design point at ignition. Replace with cold-flow measurements |
| Throat dia., expansion ratio | 72.86 mm, 5.53 | RPA run 8/16/26 |
| c*, efficiency | 1544.8 m/s × 0.974 | RPA run 8/16/26 |
| Nozzle γ, C_F efficiency | 1.25, 0.969 | fit to RPA ideal C_F / delivered sea-level thrust |
| **Vehicle dry mass** (`flight.py`) | 46.7 kg | **TBD** (value from `1DSimExample.cpp`, not the Spaceshot vehicle) |
| Aero table (`flight.py`) | `Eureka3_Machnum_RASAero2_sim.CSV` | in this folder (Eureka-3 RASAero run). Replace with the Spaceshot RASAero table |

`static_fire.py` prints the tank lengths implied by the masses and ullage at the
197.2 mm tank ID, so they can be checked against the hardware.

**About C_d·A:** in this model it lumps together *everything* between the tank and the
chamber: feed lines, valves, manifold and orifices. The Hotfire 7 notes call about
450 psi the typical N2O injector pressure, with tanks filled to about 600 psi, so the
feed system likely takes a large share of the pressure drop. A C_d·A fit from a cold flow (tank PT and flow rate)
captures that automatically.

## Validating

### Done: verification (`verify.py`), meaning the code solves its own equations correctly

- **Mass balance:** propellant leaving the tanks equals integrated injector flow, to within 0.00 %.
- **Timestep:** halving dt from 0.01 s to 0.005 s changes total impulse by 0.04 %.
- **Property table:** interpolated N2O properties match CoolProp to better than 0.002 %.
- **Solver self-tests:** `check_cores.py` checks mass and volume consistency, an
  energy round trip, flow in = flow out, and correct trends (draining or piston
  expansion cools the tank).

### Not done yet: validation, meaning the model matches the real engine

There is no SEB flight data for this vehicle. The two recent N2O/IPA piston-system
hotfires can't validate a multi-second blowdown:

- **Hotfire 6** aborted early with a data blackout.
- **Hotfire 7** (11/08) held nominal Pc of about 360 psi for only about 32 ms before
  the hard start and chamber PT failure.

The model's 345 psi at ignition is in the same range as that 360 psi, but 32 ms of
startup transient is not a validation.

**What would validate it, best first:**

1. **A full-duration N2O cold flow** (tank PT, injector PT, load cell or a tank mass
   trace). This isolates the tank model and gives the real C_d·A.
2. **A full-duration hotfire** with tank, injector and chamber PTs and load cells.
3. **Flight data**, eventually.

**How to run it:** export the test's channels from Grafana as CSV, then:

```bash
python3 validate.py hf.csv --time-col Time \
    --col "P_tank_psi=NOS Tank PT" --col "Pc_psi=Chamber PT" \
    --col "thrust_N=Load Cell Sum" --thrust-scale 9.80665 \
    --t-ignition "2026-11-08 13:41:06.200" \
    --set ox_mass=30 --set ox_fill_pressure=580psi
```

- Use `--set` so the model matches the *test article*: its fill mass, fill pressure
  and any measured C_d·A.
- The script prints bias and RMS error per channel and saves `validation.png`.
- If the tank pressure falls faster or slower than the test, that tests the
  equilibrium assumption directly. A constant offset in Pc usually means C_d·A or c*.

## How the model works

State integrated with 4th-order Runge–Kutta (`integrator.py`):

```
static fire:  s = [M_ox, U_ox, m_fuel]
flight:       s = [h, v, M_ox, U_ox, m_fuel]
```

Each derivative evaluation (`engine.py`) runs this chain:

1. **N2O tank** (`ox_tank.py`): equilibrium model. Liquid and vapor share one
   temperature T, and tank pressure = P_sat(T).
   - `phase_masses` solves the 2×2 system
     `[1 1; 1/ρ_l 1/ρ_v]·[m_l, m_v] = [M, V]`.
   - `solve` bisects for the T where `m_l·u_l + m_v·u_v = U`.
   - A negative m_l means the liquid is gone, and the burn ends.
2. **Piston-fed IPA** (`fuel_tank.py`): the fuel feed pressure is P_tank minus the
   piston friction. As IPA leaves, the N2O volume grows by `dV = dm_f/ρ_f`.
3. **Injectors** (`injector.py`): `ṁ = C_d·A·√(2ρ(P_up − P_c))`, with zero flow if
   P_up ≤ P_c.
4. **Chamber** (`chamber.py`): bisects for the P_c where
   `ṁ_ox(P_c) + ṁ_f(P_c) = P_c·A_t / c*`.
5. **Rates fed back to RK4**:
   - `dM/dt = −ṁ_ox`
   - `dU/dt = −ṁ_ox·h_l − P·dV/dt` (the leaving liquid's enthalpy, plus the work done on the piston)
   - `dm_f/dt = −ṁ_f`
6. **Thrust** (`nozzle.py`): `F = η·C_F(P_c, P_amb)·P_c·A_t`. C_F is an ideal isentropic
   nozzle, so the ambient-pressure term is built in, and `flight.py` passes P_amb from altitude.

N2O properties come from `n2o_saturation.csv`, a 0.1 K table generated with CoolProp by
`n2o_props.py`. The same CSV can be loaded by the C++ framework's `Interpolator1D`.

## Assumptions and limitations

- **Equilibrium tank:** no thermal lag between liquid and vapor, no stratification,
  and no heat from the tank walls. Fast drains can deviate from this. The Zimmerman
  et al. (2013) review compares equilibrium with non-equilibrium models.
- **Shared temperature on the piston side:** the N2O vapor above the piston is assumed
  to be at the same temperature as the main tank.
- **No vapor-only phase:** the vapor-only tail after the liquid runs out isn't
  modeled. The burn ends at liquid depletion.
- **SPI injector:** treats N2O as an incompressible liquid. Real N2O partly flashes
  in the orifice, so SPI over-predicts flow. The Dyer model is the usual upgrade.
- **Constant c\*:** it doesn't change with O/F or P_c. NASA CEA (`rocketcea`) would
  make it vary.
- **Ignored:** hydrostatic and acceleration head, startup and shutdown transients, and
  nozzle flow separation.
- **Flight sim limits:** `flight.py` is a vertical, point-mass flight. Its atmosphere
  model is simplified above 11 km.

## Files

| File | Role |
|---|---|
| `config.py` | All engine inputs, plus the derived starting temperature, tank volume and C_d·A |
| `n2o_props.py`, `n2o_saturation.csv` | N2O saturation properties |
| `ox_tank.py` | Equilibrium N2O tank |
| `fuel_tank.py` | Piston-fed and regulated fuel tanks |
| `injector.py` | SPI injector flow |
| `chamber.py` | Chamber pressure balance |
| `nozzle.py` | Thrust coefficient |
| `integrator.py` | Generic RK4 step |
| `engine.py` | Connects the components into `derivatives(t, s)` and `thrust()` |
| `static_fire.py` | Ground firing: summary, plots, CSV, optional `.eng` |
| `flight.py` | 1D flight with pressure-adjusted thrust |
| `aero.py`, `Eureka3_Machnum_RASAero2_sim.CSV` | Drag table reader and standard atmosphere (copied from the team's 1D framework) |
| `validate.py` | Model vs test-data comparison |
| `verify.py` | Numerical verification |
| `check_cores.py` | Solver self-tests |

## Next steps

1. Replace the TBD inputs with Propulsion's numbers.
2. Get C_d·A from a nitrous cold flow, and validate against the first full-duration hotfire.
3. Port into the C++ 1D framework and 6DOF. The state vector becomes
   `[h, v, M_ox, U_ox, m_fuel]`, and `n2o_saturation.csv` loads with `Interpolator1D`.
4. Model upgrades, if validation calls for them: vapor-only phase, Dyer injector, c*(O/F).

*Model and code: SEB Sims, Pressure-Adjusted Thrust project (Henry Jugan).*
