"""Reduced equations of state for fluids.

Reduced equations of state are general in the sense that they
do not involve material parameters (e.g. physico-chemical constants)
which characterize a particular, real fluid.
"""
import numpy as np
from abc import ABCMeta, abstractproperty, abstractmethod


class ReducedEquationOfState(object):

    """Abstract base class for the reduced equations of state."""

    __metaclass__ = ABCMeta

    @abstractproperty
    def Zc(self):
        """Compressibility factor at the critical point."""

    @abstractproperty
    def VrmolMin(self):
        """Lower limit for the reduced molar volume."""

    @abstractmethod
    def Pr(self, Vrmol, Tr):
        """Equation of state for the reduced pressure.

        Parameters
        ----------
        Vrmol : float
            Reduced molar volume > VrmolMin;
            Note, Vrmol = 1/reduced density
        Tr : float
            Reduced temperature > 0

        Returns
        -------
        Reduced pressure : float
        """


class ReducedIdealGas(ReducedEquationOfState):

    """Reduced ideal gas equation of state."""

    @property
    def Zc(self):
        """Compressibility factor at the critical point."""
        return 1.0

    @property
    def VrmolMin(self):
        """Lower limit for the reduced molar volume."""
        return 0.0

    def Pr(self, Vrmol, Tr):
        """Equation of state for the reduced pressure.

        Parameters
        ----------
        Vrmol : float
            Reduced molar volume > VrmolMin;
            Note, Vrmol = 1/reduced density
        Tr : float
            Reduced temperature > 0

        Returns
        -------
        Reduced pressure : float
        """
        return Tr/Vrmol


class ReducedVanDerWaals(ReducedEquationOfState):

    """Reduced Van der Waals equation of state."""

    @property
    def Zc(self):
        """Compressibility factor at the critical point."""
        return 3.0/8.0

    @property
    def VrmolMin(self):
        """Lower limit for the reduced molar volume."""
        return 1.0/3.0

    def Pr(self, Vrmol, Tr):
        """Equation of state for the reduced pressure.

        Parameters
        ----------
        Vrmol : float
            Reduced molar volume > VrmolMin;
            Note, Vrmol = 1/reduced density
        Tr : float
            Reduced temperature > 0

        Returns
        -------
        Reduced pressure : float
        """
        rep = 8*Tr/(3*Vrmol - 1)
        att = 3/(Vrmol*Vrmol)

        return rep - att


class ReducedRedlichKwong(ReducedEquationOfState):

    """Reduced Redlich-Kwong equation of state."""

    def __init__(self):
        """Precompute coefficients of the equation of state."""
        self._cff = np.power(2, 1.0/3.0) - 1

    @property
    def Zc(self):
        """Compressibility factor at the critical point."""
        return 1.0/3.0

    @property
    def VrmolMin(self):
        """Lower limit for the reduced molar volume."""
        return self._cff

    def Pr(self, Vrmol, Tr):
        """Equation of state for the reduced pressure.

        Parameters
        ----------
        Vrmol : float
            Reduced molar volume > VrmolMin;
            Note, Vrmol = 1/reduced density
        Tr : float
            Reduced temperature > 0

        Returns
        -------
        Reduced pressure : float
        """
        rep = 3*Tr/(Vrmol - self._cff)
        att = 1/(self._cff*np.sqrt(Tr)*Vrmol*(Vrmol + self._cff))
        return rep - att


class ReducedSoave(ReducedEquationOfState):

    """Reduced Soave equation of state.

    Attributes
    ----------
    Omega: float
        Acentric factor (default value for water)
    """

    def __init__(self, Omega=0.344):
        """Precompute coefficients of the equation of state."""
        self._Omega = Omega
        self._cff = np.power(2, 1.0/3.0) - 1
        self._alpha_cff = 0.48508 + 1.55171*Omega - 0.15613*Omega*Omega

    @property
    def Zc(self):
        """Compressibility factor at the critical point."""
        return 1.0/3.0

    @property
    def Omega(self):
        """Acentric factor"""
        return self._Omega

    @property
    def VrmolMin(self):
        """Lower limit for the reduced molar volume."""
        return self._cff

    def Pr(self, Vrmol, Tr):
        """Equation of state for the reduced pressure.

        Parameters
        ----------
        Vrmol : float
            Reduced molar volume > VrmolMin;
            Note, Vrmol = 1/reduced density
        Tr : float
            Reduced temperature > 0

        Returns
        -------
        Reduced pressure : float
        """
        rep = 3*Tr/(Vrmol - self._cff)

        alpha1 = 1 + self._alpha_cff*(1 - np.sqrt(Tr))
        alpha = alpha1*alpha1

        att = alpha/(self._cff*Vrmol*(Vrmol + self._cff))
        return rep - att


class ReducedPengRobinson(ReducedEquationOfState):

    """Reduced Peng-Robinson equation of state.

    Attributes
    ----------
    Omega : float
        Acentric factor (default value for water)
    """

    def __init__(self, Omega=0.344):
        """Precompute coefficients of the equation of state."""
        self._Omega = Omega
        self._cff = 0.2534
        self._alpha_cff = 0.37464 + 1.54226*Omega - 0.26992*Omega*Omega

    @property
    def Zc(self):
        """Compressibility factor at the critical point."""
        return 0.307

    @property
    def Omega(self):
        """Acentric factor"""
        return self._Omega

    @property
    def VrmolMin(self):
        """Lower limit for the reduced molar volume."""
        return self._cff

    def Pr(self, Vrmol, Tr):
        """Equation of state for the reduced pressure.

        Parameters
        ----------
        Vrmol : float
            Reduced molar volume > VrmolMin;
            Note, Vrmol = 1/reduced density
        Tr : float
            Reduced temperature > 0

        Returns
        -------
        Reduced pressure : float
        """
        rep = 3.2573*Tr/(Vrmol - self._cff)

        alpha1 = 1 + self._alpha_cff*(1 - np.sqrt(Tr))
        alpha = 4.8514*alpha1*alpha1

        att = alpha/(Vrmol*Vrmol + 2*self._cff*Vrmol - self._cff*self._cff)
        return rep - att


class ReducedCarnahanStarling(ReducedEquationOfState):

    """Reduced Carnahan-Starling equation of state."""

    def __init__(self):
        """Precompute coefficients of the equation of state."""
        self._cff1 = 2.785855166
        self._cff2 = 0.1304438842
        self._cff3 = 3.852462257

    @property
    def Zc(self):
        """Compressibility factor at the critical point."""
        return 0.3589562

    @property
    def VrmolMin(self):
        """Lower limit for the reduced molar volume."""
        return self._cff2

    def Pr(self, Vrmol, Tr):
        """Equation of state for the reduced pressure.

        Parameters
        ----------
        Vrmol : float
            Reduced molar volume > VrmolMin;
            Note, Vrmol = 1/reduced density
        Tr : float
            Reduced temperature > 0

        Returns
        -------
        Reduced pressure : float
        """
        aux1 = self._cff2/Vrmol
        aux2 = aux1*aux1

        aux3 = 1 - aux1
        aux4 = aux3*aux3*aux3

        rep1 = self._cff1*Tr/Vrmol
        rep2 = (1 + aux1 + aux2 - aux1*aux2)/aux4
        rep = rep1*rep2

        att = self._cff3/(Vrmol*Vrmol)
        return rep - att
