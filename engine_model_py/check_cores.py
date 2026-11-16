#CREDIT TO CLAUDE FOR THIS FILE!

"""Self-test for the three physics solvers (injector, N2O tank, chamber).

Run this after changing any physics code.

    python3 check_cores.py            # all cores, in order
    python3 check_cores.py injector   # injector flow only
    python3 check_cores.py tank       # N2O tank (phase_masses + solve)
    python3 check_cores.py chamber    # chamber pressure balance

These checks don't compare to a hidden answer key. They test whether the
result obeys the physics.
mass and volume add up, energy round-trips, flow in equals flow out, etc.
If a check fails, the message says which law was broken.
"""
import math
import sys

from n2o_props import SaturationTable

PASS, FAIL = 'PASS', 'FAIL'
results = []


def check(name, condition, hint=''):
    results.append(condition)
    print(f'  [{PASS if condition else FAIL}] {name}' + ('' if condition else f'\n         -> {hint}'))


def close(a, b, rel=1e-6, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def run_section(title, fn):
    print(f'\n{title}')
    try:
        fn()
    except NotImplementedError as e:
        results.append(False)
        print(f'  [TODO] {e}')
    except Exception as e:  # show the error but keep checking other cores
        results.append(False)
        print(f'  [ERROR] {type(e).__name__}: {e}')


def check_injector():
    from injector import spi_mass_flow
    m = spi_mass_flow(1e-4, 800.0, 11e5, 1e5)
    check('magnitude: CdA=1e-4 m^2, rho=800, dP=10 bar gives 4.0 kg/s', close(m, 4.0, 1e-9),
          f'got {m}. Did you use dP = P_up - P_down and a factor of 2 from 1/2 rho v^2?')
    check('no flow when dP = 0', spi_mass_flow(1e-4, 800.0, 5e5, 5e5) == 0.0, 'should return 0.0')
    check('no backflow when P_down > P_up', spi_mass_flow(1e-4, 800.0, 5e5, 6e5) == 0.0,
          'should return 0.0, not an error or a negative number')
    r = spi_mass_flow(1e-4, 800.0, 41e5, 1e5) / spi_mass_flow(1e-4, 800.0, 11e5, 1e5)
    check('4x the pressure drop gives 2x the flow (square-root law)', close(r, 2.0, 1e-9), f'ratio was {r}')
    r = spi_mass_flow(2e-4, 800.0, 11e5, 1e5) / spi_mass_flow(1e-4, 800.0, 11e5, 1e5)
    check('flow is proportional to CdA', close(r, 2.0, 1e-9), f'ratio was {r}')


def check_tank():
    from ox_tank import EquilibriumOxTank
    props = SaturationTable()
    V, M = 0.045, 32.0
    tank = EquilibriumOxTank(V, M, 293.15, props)

    print('  phase_masses():')
    for T in (260.0, 280.0, 293.15, 300.0):
        ml, mv = tank.phase_masses(T, M, V)
        vol = ml / props('rho_l', T) + mv / props('rho_v', T)
        check(f'T={T:.2f} K: masses add to M and volumes add to V',
              close(ml + mv, M, 1e-9) and close(vol, V, 1e-9),
              f'm_l+m_v={ml+mv:.6f} (want {M}), volume={vol:.6f} (want {V}). Check your A matrix rows.')
    ml, mv = tank.phase_masses(293.15, M, V)
    check('at 20 C, most of the mass is liquid (about 31 kg)', 30.0 < ml < 32.0, f'got m_l={ml:.3f} kg')

    print('  solve():')
    for T0 in (265.0, 280.0, 293.15, 300.0):
        ml, mv = tank.phase_masses(T0, M, V)
        U = ml * props('u_l', T0) + mv * props('u_v', T0)
        st = tank.solve(M, U, V)
        check(f'round trip at T0={T0:.2f} K: solve(M, U(T0)) recovers T0',
              abs(st.T - T0) < 1e-3 and close(st.P, props('P', T0), 1e-4)
              and abs(st.m_liq - ml) < 1e-3 and close(st.h_liq, props('h_l', T0), 1e-4),
              f'got T={st.T:.4f} K, P={st.P:.1f} Pa, m_liq={st.m_liq:.4f}. '
              'Is your bracket/sign logic in the bisection right?')
    st0 = tank.solve(tank.M0, tank.U0, V)
    h = st0.h_liq
    st1 = tank.solve(tank.M0 - 1.0, tank.U0 - 1.0 * h, V)
    check('draining 1 kg of liquid cools the tank and drops its pressure',
          st1.T < st0.T and st1.P < st0.P,
          f'T: {st0.T:.3f} -> {st1.T:.3f} K. Removing liquid (carrying h_l) should make the rest boil and cool.')
    print(f'         (1 kg drained: {st0.T-273.15:.2f} C -> {st1.T-273.15:.2f} C, '
          f'{st0.P/1e5:.2f} -> {st1.P/1e5:.2f} bar)')
    dV = 0.001
    st2 = tank.solve(tank.M0, tank.U0 - st0.P * dV, V + dV)
    check('piston sweeps 1 L (N2O expands and does work P*dV): tank cools and pressure drops',
          st2.T < st0.T and st2.P < st0.P and st2.m_vap > st0.m_vap,
          f'T: {st0.T:.3f} -> {st2.T:.3f} K, m_vap: {st0.m_vap:.3f} -> {st2.m_vap:.3f} kg. '
          'Is V passed through to phase_masses?')
    print(f'         (1 L swept: {st0.T-273.15:.2f} C -> {st2.T-273.15:.2f} C, '
          f'{st0.P/1e5:.2f} -> {st2.P/1e5:.2f} bar)')


def check_chamber():
    from chamber import solve_chamber
    from injector import spi_mass_flow
    args = dict(P_ox=50e5, rho_ox=785.0, P_fuel=40e5, rho_fuel=789.0,
                CdA_ox=6.8e-5, CdA_fuel=2.1e-5, A_throat=math.pi * 0.052 ** 2 / 4, c_star=1350.0)
    sol = solve_chamber(**args)
    m_in = (spi_mass_flow(args['CdA_ox'], args['rho_ox'], args['P_ox'], sol.Pc)
            + spi_mass_flow(args['CdA_fuel'], args['rho_fuel'], args['P_fuel'], sol.Pc))
    m_out = sol.Pc * args['A_throat'] / args['c_star']
    check('flow in == flow out at the returned Pc', close(m_in, m_out, 1e-4),
          f'in={m_in:.5f} kg/s, out={m_out:.5f} kg/s at Pc={sol.Pc/1e5:.3f} bar')
    check('Pc is between 0 and the highest feed pressure', 0.0 < sol.Pc < max(args['P_ox'], args['P_fuel']),
          f'Pc={sol.Pc}')
    check('returned flows are the injector flows at that Pc',
          close(sol.mdot_ox, spi_mass_flow(args['CdA_ox'], args['rho_ox'], args['P_ox'], sol.Pc), 1e-9)
          and close(sol.mdot_fuel, spi_mass_flow(args['CdA_fuel'], args['rho_fuel'], args['P_fuel'], sol.Pc), 1e-9),
          'evaluate spi_mass_flow at the final Pc for each side')
    lower = solve_chamber(**{**args, 'P_ox': 40e5})
    check('lower N2O tank pressure -> lower Pc, less N2O flow, MORE fuel flow',
          lower.Pc < sol.Pc and lower.mdot_ox < sol.mdot_ox and lower.mdot_fuel > sol.mdot_fuel,
          'with regulated fuel, a lower Pc means a bigger fuel-side pressure drop')
    empty = solve_chamber(**{**args, 'P_ox': 0.0, 'P_fuel': 0.0})
    check('everything empty -> Pc = 0 and no flow', empty.Pc == 0.0 and empty.mdot == 0.0, f'got {empty}')
    print(f'         (Pc = {sol.Pc/1e5:.2f} bar, O/F = {sol.OF:.2f})')


SECTIONS = {
    'injector': ('Injector: injector.spi_mass_flow', check_injector),
    'tank': ('N2O tank: ox_tank phase_masses + solve', check_tank),
    'chamber': ('Chamber: chamber.solve_chamber', check_chamber),
}

if __name__ == '__main__':
    chosen = sys.argv[1:] or list(SECTIONS)
    for key in chosen:
        run_section(*SECTIONS[key])
    passed = sum(results)
    print(f'\n{passed}/{len(results)} checks passed')
    if passed == len(results):
        print('All good: try  python3 static_fire.py  or  python3 flight.py')
