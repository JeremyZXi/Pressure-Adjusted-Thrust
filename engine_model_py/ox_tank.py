"""Assumes liquid & vapor are at the same temp inside the tank & that both are saturated
this means that tank pressure = vapor pressure P_sat(T)
ALso no heat exchange at the barriers of the tank

    M = total N2O mass in the tank (liquid + vapor)
    U = total internal energy of all N2O in the tank

M and U change because liquid leaves through the injector and the  N2O 
does work pushing the piston down:
    dM/dt = -mdot_ox
    dU/dt = -mdot_ox * h_l  - P * dV/dt
            leaving liquid carries enthalpy h_l = u_l + P/rho_l;
             P*dV/dt is the work done on the piston

T, P, how much is liquid vs vapor is then solved for using M, U & V
^^ this is what what solve() does

"""
from dataclasses import dataclass

import numpy as np


@dataclass
class TankState:
    T: float       # temperature [K]
    P: float       # pressure [Pa]
    m_liq: float   # liquid mass [kg]  (<= 0 means the liquid has run out)
    m_vap: float   # vapor mass [kg]
    rho_liq: float # liquid density [kg/m^3] (what the injector sees)
    h_liq: float   # liquid specific enthalpy [J/kg] (what leaves the tank)

    @property
    def liquid_gone(self):
        return self.m_liq <= 0.0 # if the liquid mass is below zero then there is no liquid left
    #this would indicate the end of the burn


class EquilibriumOxTank:
    def __init__(self, volume, ox_mass, T0, props):
        """volume at fill, ox_mass loaded, initial temperature T0,
        props: n2o_props.SaturationTable"""
        self.V0 = volume
        self.props = props
        density = ox_mass / volume
        if not props('rho_v', T0) < density < props('rho_l', T0):
            raise ValueError(
                f'Fill density {density:.0f} kg/m^3 must be between vapor '
                f'({props("rho_v", T0):.0f}) and liquid ({props("rho_l", T0):.0f}) '
                f'density at T0 = {T0} K, or the tank is not liquid + vapor.')
        m_l, m_v = self.phase_masses(T0, ox_mass, volume)
        self.M0 = ox_mass
        self.U0 = m_l * props('u_l', T0) + m_v * props('u_v', T0)

    def phase_masses(self, T, M, V):
        #m_l/rho_l + m_v/rho_v   = V
        #this gives the linear equation 1/rho * m = V
        # m_l + m_v = M
        # this gives linear equation 1 * m = M
        #this gives 2 equations for the unknowns that can be solved as a 2x2 linear system
        #follows the form Ax = b
        #this then gives m_l & m_v (x in the lin)
        #ezpz

        rho_l = self.props('rho_l', T)
        rho_v = self.props('rho_v', T)

        A = np.array([[1, 1],
              [1 / rho_l, 1 / rho_v]])

        b = np.array([M, V])

        x = np.linalg.solve(A, b)

        m_l = x[0]
        m_v = x[1]
        return m_l, m_v

    def solve(self, M, U, V):

        #How this owrks is basically we guess a temp T & phase_masses(T, M, V)
        #tells us how much the mass would split at that T
        #then we can compute teh internal energy at that T
        #we subtract this calculate internal energy from the actual internal energy
        #U is the current internal energy tracked by the integrator 
        # it starts at the fill value U0
        #U(T) - U(known) = 0 tells us what the proper temp for the tank state is
        #U(T) is the internal energy derived from teh guessed temp
        def f(T):
            m_l, m_v = self.phase_masses(T, M, V)
            u_l = self.props('u_l', T)
            u_v = self.props('u_v', T)
            return (m_l * u_l) + (m_v * u_v) - U   # if this = 0 we have the correct T
        lo, hi = self.props.T_min, self.props.T_max
        while hi - lo > 1e-6:
            mid = 0.5 * (lo + hi)
            if f(mid) > 0: #the root must be in the lower half 
                hi = mid
            else: 
                lo = mid
        T = mid
        m_l, m_v = self.phase_masses(T, M, V)
        
        return TankState(T, self.props('P', T), m_l, m_v, self.props('rho_l', T), self.props('h_l', T))
        

