#!/usr/local/bin/python
# -*- coding: utf-8 -*-

######################################################################################################################################################
#								IMPORT PACKAGES
######################################################################################################################################################
import sys
import moments
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import numpy


######################################################################################################################################################
#							STRICT ISOLATION MODELS - 2 POPULATIONS
######################################################################################################################################################
def SI(params, ns):
    nu1, nu2, Ts, O = params
    n1, n2 = ns
    """
    Strict Isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: Duration of divergence in isolation.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Equilibrium ancestral population
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phi = moments.Spectrum(sts)
    # Split event
    phi = moments.Manips.split_1D_to_2D(phi, n1, n2)
    # Integrate by setting divergence time to Ts, the population sizes to nu1 and nu2, and the migration rates to zero
    phi.integrate([nu1, nu2], Ts)
    # Calculate the well-oriented spectrum
    fsO = phi
    # Calculate the mis-oriented spectrum
    fsM = moments.Numerics.reverse_array(fsO)

    # Sum all spectra
    fs = O*fsO + (1-O)*fsM
    return fs


def SIexp(params, ns):
    nu1e, nu1, nu2e, nu2, Ts, Ti, O = params
    n1, n2 = ns
    """
    Strict Isolation with exponential size changes.

    nu1e: Size of population 1 after split.
    nu1: Size of population 1 after exponential size change.
    nu2e: Size of population 2 after split.
    nu2: Size of population 2 after exponential size change.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence in isolation with pop1: exponential change / pop2: exponential change.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phi = moments.Spectrum(sts)
    phi = moments.Manips.split_1D_to_2D(phi, n1, n2)
    phi.integrate([nu1e, nu2e], Ts)
    nu1_func = lambda t: nu1e*(nu1/nu1e) ** (t/Ti)
    nu2_func = lambda t: nu2e*(nu2/nu2e) ** (t/Ti)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    phi.integrate(nu_func, Ti)
    fsO = phi
    fsM = moments.Numerics.reverse_array(fsO)

    # Sum all spectra
    fs = O*fsO + (1-O)*fsM
    return fs


def SI2N(params, ns):
    nu1, nu2, Ts, nr, bf, O = params
    n1, n2 = ns
    """
    Strict Isolation with heterogeneous Ne (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: Duration of divergence in isolation.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phinr = moments.Spectrum(sts)
    phinr = moments.Manips.split_1D_to_2D(phinr, n1, n2)
    phinr.integrate([nu1*bf, nu2*bf], Ts)
    fsnrO = phinr
    fsnrM = moments.Numerics.reverse_array(fsnrO)
    
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phir = moments.Spectrum(sts)
    phir = moments.Manips.split_1D_to_2D(phir, n1, n2)
    phir.integrate([nu1, nu2], Ts)
    fsrO = phir
    fsrM = moments.Numerics.reverse_array(fsrO)

    # Sum all spectra
    fs = O*( nr*fsnrO + (1-nr)*fsrO ) + (1-O)*( nr*fsnrM + (1-nr)*fsrM )
    return fs


######################################################################################################################################################
#							ISOLATION WITH MIGRATION MODELS - 2 POPULATIONS
######################################################################################################################################################
def IM(params, ns):
    nu1, nu2, m12, m21, Ts, O = params
    n1, n2 = ns
    """
    Isolation with Migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence with migration.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phi = moments.Spectrum(sts)
    phi = moments.Manips.split_1D_to_2D(phi, n1, n2)
    phi.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsO = phi
    fsM = moments.Numerics.reverse_array(fsO)

    # Sum all spectra
    fs = O*fsO + (1-O)*fsM
    return fs


def IMexp(params, ns):
    nu1e, nu1, nu2e, nu2, m12, m21, Ts, Ti, O = params
    n1, n2 = ns
    """
    Isolation with Migration with exponential size changes.

    nu1e: Size of population 1 after split.
    nu1: Size of population 1 after exponential size change.
    nu2e: Size of population 2 after split.
    nu2: Size of population 2 after exponential size change.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence with migration.
    Ti: Duration of divergence with migration with pop1: exponential change / pop2: exponential change.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phi = moments.Spectrum(sts)
    phi = moments.Manips.split_1D_to_2D(phi, n1, n2)
    phi.integrate([nu1e, nu2e], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    nu1_func = lambda t: nu1e*(nu1/nu1e) ** (t/Ti)
    nu2_func = lambda t: nu2e*(nu2/nu2e) ** (t/Ti)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    phi.integrate(nu_func, Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsO = phi
    fsM = moments.Numerics.reverse_array(fsO)

    # Sum all spectra
    fs = O*fsO + (1-O)*fsM
    return fs


def IM2N(params, ns):
    nu1, nu2, m12, m21, Ts, nr, bf, O = params
    n1, n2 = ns
    """
    Isolation with Migration with heterogeneous Ne (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence with migration.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phinr = moments.Spectrum(sts)
    phinr = moments.Manips.split_1D_to_2D(phinr, n1, n2)
    phinr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsnrO = phinr
    fsnrM = moments.Numerics.reverse_array(fsnrO)
    
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phir = moments.Spectrum(sts)
    phir = moments.Manips.split_1D_to_2D(phir, n1, n2)
    phir.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsrO = phir
    fsrM = moments.Numerics.reverse_array(fsrO)

    # Sum all spectra
    fs = O*( nr*fsnrO + (1-nr)*fsrO ) + (1-O)*( nr*fsnrM + (1-nr)*fsrM )
    return fs 


def IM2M(params, ns):
    nu1, nu2, m12, m21, Ts, P, O = params
    n1, n2 = ns
    """
    Isolation with Migration with heterogeneous Me (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence with migration.
    P: Fraction of the genome evolving neutrally.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN = moments.Spectrum(sts)
    phiN = moments.Manips.split_1D_to_2D(phiN, n1, n2)
    phiN.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsNO = phiN
    fsNM = moments.Numerics.reverse_array(fsNO)

    # Barrier spectrum
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI = moments.Spectrum(sts)
    phiI = moments.Manips.split_1D_to_2D(phiI, n1, n2)
    phiI.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    fsIO = phiI
    fsIM = moments.Numerics.reverse_array(fsIO)

    # Sum all spectra
    fs = O*( P*fsNO + (1-P)*fsIO ) + (1-O)*( P*fsNM + (1-P)*fsIM )
    return fs
    

def IM2M2P(params, ns):
    nu1, nu2, m12, m21, Ts, P1, P2, O = params
    n1, n2 = ns
    """
    Isolation with Migration with heterogeneous Me (population-specific).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence with migration.
    P1: Fraction of the genome evolving neutrally in population 1.
    P2: Fraction of the genome evolving neutrally in population 2.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum in population 1 and 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2 = moments.Spectrum(sts)
    phiN1N2 = moments.Manips.split_1D_to_2D(phiN1N2, n1, n2)
    phiN1N2.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsN1N2O = phiN1N2
    fsN1N2M = moments.Numerics.reverse_array(fsN1N2O)

    # Barrier spectrum in population 1 and 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2 = moments.Spectrum(sts)
    phiI1I2 = moments.Manips.split_1D_to_2D(phiI1I2, n1, n2)
    phiI1I2.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2O = phiI1I2
    fsI1I2M = moments.Numerics.reverse_array(fsI1I2O)

    # Neutral spectrum in population 1, barrier spectrum in population 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2 = moments.Spectrum(sts)
    phiN1I2 = moments.Manips.split_1D_to_2D(phiN1I2, n1, n2)
    phiN1I2.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[0, 0]]))
    fsN1I2O = phiN1I2
    fsN1I2M = moments.Numerics.reverse_array(fsN1I2O)

    # Neutral spectrum in population 2, barrier spectrum in population 1
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2 = moments.Spectrum(sts)
    phiI1N2 = moments.Manips.split_1D_to_2D(phiI1N2, n1, n2)
    phiI1N2.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[m21, 0]]))
    fsI1N2O = phiI1N2
    fsI1N2M = moments.Numerics.reverse_array(fsI1N2O)

    # Sum all spectra
    fs = O*( P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O ) + (1-O)*( P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M )
    return fs


def IM2N2M(params, ns):
    nu1, nu2, m12, m21, Ts, nr, bf, P, O = params
    n1, n2 = ns
    """
    Isolation with Migration with heterogeneous Ne (shared between populations) and heterogeneous Me (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence with migration.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    P: Fraction of the genome evolving neutrally.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiNnr = moments.Spectrum(sts)
    phiNnr = moments.Manips.split_1D_to_2D(phiNnr, n1, n2)
    phiNnr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsNOnr = phiNnr
    fsNMnr = moments.Numerics.reverse_array(fsNOnr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiNr = moments.Spectrum(sts)
    phiNr = moments.Manips.split_1D_to_2D(phiNr, n1, n2)
    phiNr.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsNOr = phiNr
    fsNMr = moments.Numerics.reverse_array(fsNOr)

    # Barrier spectrum
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiInr = moments.Spectrum(sts)
    phiInr = moments.Manips.split_1D_to_2D(phiInr, n1, n2)
    phiInr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    fsIOnr = phiInr
    fsIMnr = moments.Numerics.reverse_array(fsIOnr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiIr = moments.Spectrum(sts)
    phiIr = moments.Manips.split_1D_to_2D(phiIr, n1, n2)
    phiIr.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    fsIOr = phiIr
    fsIMr = moments.Numerics.reverse_array(fsIOr)

    # Sum all spectra
    fs = O*( P*nr*fsNOnr + P*(1-nr)*fsNOr + (1-P)*nr*fsIOnr + (1-P)*(1-nr)*fsIOr ) + (1-O)*( P*nr*fsNMnr + P*(1-nr)*fsNMr + (1-P)*nr*fsIMnr + (1-P)*(1-nr)*fsIMr )
    return fs


def IM2N2M2P(params, ns):
    nu1, nu2, m12, m21, Ts, nr, bf, P1, P2, O = params
    n1, n2 = ns
    """
    Isolation with Migration with heterogeneous Ne (shared between populations) and heterogeneous Me (population-specific).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence with migration.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    P1: Fraction of the genome evolving neutrally in population 1.
    P2: Fraction of the genome evolving neutrally in population 2.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum in population 1 and 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2nr = moments.Spectrum(sts)
    phiN1N2nr = moments.Manips.split_1D_to_2D(phiN1N2nr, n1, n2)
    phiN1N2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsN1N2Onr = phiN1N2nr
    fsN1N2Mnr = moments.Numerics.reverse_array(fsN1N2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2r = moments.Spectrum(sts)
    phiN1N2r = moments.Manips.split_1D_to_2D(phiN1N2r, n1, n2)
    phiN1N2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    fsN1N2Or = phiN1N2r
    fsN1N2Mr = moments.Numerics.reverse_array(fsN1N2Or)

    # Barrier spectrum in population 1 and 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2nr = moments.Spectrum(sts)
    phiI1I2nr = moments.Manips.split_1D_to_2D(phiI1I2nr, n1, n2)
    phiI1I2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2Onr = phiI1I2nr
    fsI1I2Mnr = moments.Numerics.reverse_array(fsI1I2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2r = moments.Spectrum(sts)
    phiI1I2r = moments.Manips.split_1D_to_2D(phiI1I2r, n1, n2)
    phiI1I2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2Or = phiI1I2r
    fsI1I2Mr = moments.Numerics.reverse_array(fsI1I2Or)

    # Neutral spectrum in population 1, barrier spectrum in population 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2nr = moments.Spectrum(sts)
    phiN1I2nr = moments.Manips.split_1D_to_2D(phiN1I2nr, n1, n2)
    phiN1I2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, m12],[0, 0]]))
    fsN1I2Onr = phiN1I2nr
    fsN1I2Mnr = moments.Numerics.reverse_array(fsN1I2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2r = moments.Spectrum(sts)
    phiN1I2r = moments.Manips.split_1D_to_2D(phiN1I2r, n1, n2)
    phiN1I2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[0, 0]]))
    fsN1I2Or = phiN1I2r
    fsN1I2Mr = moments.Numerics.reverse_array(fsN1I2Or)

    # Neutral spectrum in population 2, barrier spectrum in population 1
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2nr = moments.Spectrum(sts)
    phiI1N2nr = moments.Manips.split_1D_to_2D(phiI1N2nr, n1, n2)
    phiI1N2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[m21, 0]]))
    fsI1N2Onr = phiI1N2nr
    fsI1N2Mnr = moments.Numerics.reverse_array(fsI1N2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2r = moments.Spectrum(sts)
    phiI1N2r = moments.Manips.split_1D_to_2D(phiI1N2r, n1, n2)
    phiI1N2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[m21, 0]]))
    fsI1N2Or = phiI1N2r
    fsI1N2Mr = moments.Numerics.reverse_array(fsI1N2Or)

    # Sum all spectra
    fs = O*( P1*P2*nr*fsN1N2Onr + P1*P2*(1-nr)*fsN1N2Or + (1-P1)*(1-P2)*nr*fsI1I2Onr + (1-P1)*(1-P2)*(1-nr)*fsI1I2Or + P1*(1-P2)*nr*fsN1I2Onr + P1*(1-P2)*(1-nr)*fsN1I2Or + (1-P1)*P2*nr*fsI1N2Onr + (1-P1)*P2*(1-nr)*fsI1N2Or ) + (1-O)*( P1*P2*nr*fsN1N2Mnr + P1*P2*(1-nr)*fsN1N2Mr + (1-P1)*(1-P2)*nr*fsI1I2Mnr + (1-P1)*(1-P2)*(1-nr)*fsI1I2Mr + P1*(1-P2)*nr*fsN1I2Mnr + P1*(1-P2)*(1-nr)*fsN1I2Mr + (1-P1)*P2*nr*fsI1N2Mnr + (1-P1)*P2*(1-nr)*fsI1N2Mr )
    return fs


######################################################################################################################################################
#							SECONDARY CONTACT MODELS - 2 POPULATIONS
######################################################################################################################################################
def SC(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, O = params
    n1, n2 = ns
    """
    Secondary Contact.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence with migration.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phi = moments.Spectrum(sts)
    phi = moments.Manips.split_1D_to_2D(phi, n1, n2)
    phi.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phi.integrate([nu1, nu2], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsO = phi
    fsM = moments.Numerics.reverse_array(fsO)

    # Sum all spectra
    fs = O*fsO + (1-O)*fsM
    return fs


def SCexp(params, ns):
    nu1e, nu1, nu2e, nu2, m12, m21, Ts, Ti, O = params
    n1, n2 = ns
    """
    Secondary Contact with exponential size changes.

    nu1e: Size of population 1 after split.
    nu1: Size of population 1 after exponential size change.
    nu2e: Size of population 2 after split.
    nu2: Size of population 2 after exponential size change.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence with migration with pop1: exponential change / pop2: exponential change.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phi = moments.Spectrum(sts)
    phi = moments.Manips.split_1D_to_2D(phi, n1, n2)
    phi.integrate([nu1e, nu2e], Ts, m=numpy.array([[0, 0],[0, 0]]))
    nu1_func = lambda t: nu1e*(nu1/nu1e) ** (t/Ti)
    nu2_func = lambda t: nu2e*(nu2/nu2e) ** (t/Ti)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    phi.integrate(nu_func, Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsO = phi
    fsM = moments.Numerics.reverse_array(fsO)

    # Sum all spectra
    fs = O*fsO + (1-O)*fsM
    return fs


def SC2N(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, nr, bf, O = params
    n1, n2 = ns
    """
    Secondary Contact with heterogeneous Ne (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence with migration.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phinr = moments.Spectrum(sts)
    phinr = moments.Manips.split_1D_to_2D(phinr, n1, n2)
    phinr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phinr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsnrO = phinr
    fsnrM = moments.Numerics.reverse_array(fsnrO)
    
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phir = moments.Spectrum(sts)
    phir = moments.Manips.split_1D_to_2D(phir, n1, n2)
    phir.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phir.integrate([nu1, nu2], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsrO = phir
    fsrM = moments.Numerics.reverse_array(fsrO)

    # Sum all spectra
    fs = O*( nr*fsnrO + (1-nr)*fsrO ) + (1-O)*( nr*fsnrM + (1-nr)*fsrM )
    return fs


def SC2M(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, P, O = params
    n1, n2 = ns
    """
    Secondary Contact with heterogeneous Me (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence with migration.
    P: Fraction of the genome evolving neutrally.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN = moments.Spectrum(sts)
    phiN = moments.Manips.split_1D_to_2D(phiN, n1, n2)
    phiN.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiN.integrate([nu1, nu2], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsNO = phiN
    fsNM = moments.Numerics.reverse_array(fsNO)

    # Barrier spectrum
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI = moments.Spectrum(sts)
    phiI = moments.Manips.split_1D_to_2D(phiI, n1, n2)
    phiI.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsIO = phiI
    fsIM = moments.Numerics.reverse_array(fsIO)

    # Sum all spectra
    fs = O*( P*fsNO + (1-P)*fsIO ) + (1-O)*( P*fsNM + (1-P)*fsIM )
    return fs


def SC2M2P(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, P1, P2, O = params
    n1, n2 = ns
    """
    Secondary Contact with heterogeneous Me (population-specific).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence with migration.
    P1: Fraction of the genome evolving neutrally in population 1.
    P2: Fraction of the genome evolving neutrally in population 2.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum in population 1 and 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2 = moments.Spectrum(sts)
    phiN1N2 = moments.Manips.split_1D_to_2D(phiN1N2, n1, n2)
    phiN1N2.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiN1N2.integrate([nu1, nu2], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsN1N2O = phiN1N2
    fsN1N2M = moments.Numerics.reverse_array(fsN1N2O)

    # Barrier spectrum in population 1 and 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2 = moments.Spectrum(sts)
    phiI1I2 = moments.Manips.split_1D_to_2D(phiI1I2, n1, n2)
    phiI1I2.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1I2.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2O = phiI1I2
    fsI1I2M = moments.Numerics.reverse_array(fsI1I2O)

    # Neutral spectrum in population 1, barrier spectrum in population 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2 = moments.Spectrum(sts)
    phiN1I2 = moments.Manips.split_1D_to_2D(phiN1I2, n1, n2)
    phiN1I2.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiN1I2.integrate([nu1, nu2], Ti, m=numpy.array([[0, m12],[0, 0]]))
    fsN1I2O = phiN1I2
    fsN1I2M = moments.Numerics.reverse_array(fsN1I2O)

    # Neutral spectrum in population 2, barrier spectrum in population 1
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2 = moments.Spectrum(sts)
    phiI1N2 = moments.Manips.split_1D_to_2D(phiI1N2, n1, n2)
    phiI1N2.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1N2.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[m21, 0]]))
    fsI1N2O = phiI1N2
    fsI1N2M = moments.Numerics.reverse_array(fsI1N2O)

    # Sum all spectra
    fs = O*( P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O ) + (1-O)*( P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M )
    return fs


def SC2N2M(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, nr, bf, P, O = params
    n1, n2 = ns
    """
    Secondary Contact with heterogeneous Ne (shared between populations) and heterogeneous Me (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence with migration.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    P: Fraction of the genome evolving neutrally.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiNnr = moments.Spectrum(sts)
    phiNnr = moments.Manips.split_1D_to_2D(phiNnr, n1, n2)
    phiNnr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiNnr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsNOnr = phiNnr
    fsNMnr = moments.Numerics.reverse_array(fsNOnr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiNr = moments.Spectrum(sts)
    phiNr = moments.Manips.split_1D_to_2D(phiNr, n1, n2)
    phiNr.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiNr.integrate([nu1, nu2], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsNOr = phiNr
    fsNMr = moments.Numerics.reverse_array(fsNOr)

    # Barrier spectrum
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiInr = moments.Spectrum(sts)
    phiInr = moments.Manips.split_1D_to_2D(phiInr, n1, n2)
    phiInr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiInr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsIOnr = phiInr
    fsIMnr = moments.Numerics.reverse_array(fsIOnr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiIr = moments.Spectrum(sts)
    phiIr = moments.Manips.split_1D_to_2D(phiIr, n1, n2)
    phiIr.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiIr.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsIOr = phiIr
    fsIMr = moments.Numerics.reverse_array(fsIOr)

    # Sum all spectra
    fs = O*( P*nr*fsNOnr + P*(1-nr)*fsNOr + (1-P)*nr*fsIOnr + (1-P)*(1-nr)*fsIOr ) + (1-O)*( P*nr*fsNMnr + P*(1-nr)*fsNMr + (1-P)*nr*fsIMnr + (1-P)*(1-nr)*fsIMr )
    return fs


def SC2N2M2P(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, nr, bf, P1, P2, O = params
    n1, n2 = ns
    """
    Secondary Contact with heterogeneous Ne (shared between populations) and heterogeneous Me (population-specific).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence with migration.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    P1: Fraction of the genome evolving neutrally in population 1.
    P2: Fraction of the genome evolving neutrally in population 2.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum in population 1 and 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2nr = moments.Spectrum(sts)
    phiN1N2nr = moments.Manips.split_1D_to_2D(phiN1N2nr, n1, n2)
    phiN1N2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiN1N2nr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsN1N2Onr = phiN1N2nr
    fsN1N2Mnr = moments.Numerics.reverse_array(fsN1N2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2r = moments.Spectrum(sts)
    phiN1N2r = moments.Manips.split_1D_to_2D(phiN1N2r, n1, n2)
    phiN1N2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiN1N2r.integrate([nu1, nu2], Ti, m=numpy.array([[0, m12],[m21, 0]]))
    fsN1N2Or = phiN1N2r
    fsN1N2Mr = moments.Numerics.reverse_array(fsN1N2Or)

    # Barrier spectrum in population 1 and 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2nr = moments.Spectrum(sts)
    phiI1I2nr = moments.Manips.split_1D_to_2D(phiI1I2nr, n1, n2)
    phiI1I2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1I2nr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2Onr = phiI1I2nr
    fsI1I2Mnr = moments.Numerics.reverse_array(fsI1I2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2r = moments.Spectrum(sts)
    phiI1I2r = moments.Manips.split_1D_to_2D(phiI1I2r, n1, n2)
    phiI1I2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1I2r.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2Or = phiI1I2r
    fsI1I2Mr = moments.Numerics.reverse_array(fsI1I2Or)

    # Neutral spectrum in population 1, barrier spectrum in population 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2nr = moments.Spectrum(sts)
    phiN1I2nr = moments.Manips.split_1D_to_2D(phiN1I2nr, n1, n2)
    phiN1I2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiN1I2nr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, m12],[0, 0]]))
    fsN1I2Onr = phiN1I2nr
    fsN1I2Mnr = moments.Numerics.reverse_array(fsN1I2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2r = moments.Spectrum(sts)
    phiN1I2r = moments.Manips.split_1D_to_2D(phiN1I2r, n1, n2)
    phiN1I2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiN1I2r.integrate([nu1, nu2], Ti, m=numpy.array([[0, m12],[0, 0]]))
    fsN1I2Or = phiN1I2r
    fsN1I2Mr = moments.Numerics.reverse_array(fsN1I2Or)

    # Neutral spectrum in population 2, barrier spectrum in population 1
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2nr = moments.Spectrum(sts)
    phiI1N2nr = moments.Manips.split_1D_to_2D(phiI1N2nr, n1, n2)
    phiI1N2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1N2nr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[m21, 0]]))
    fsI1N2Onr = phiI1N2nr
    fsI1N2Mnr = moments.Numerics.reverse_array(fsI1N2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2r = moments.Spectrum(sts)
    phiI1N2r = moments.Manips.split_1D_to_2D(phiI1N2r, n1, n2)
    phiI1N2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1N2r.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[m21, 0]]))
    fsI1N2Or = phiI1N2r
    fsI1N2Mr = moments.Numerics.reverse_array(fsI1N2Or)

    # Sum all spectra
    fs = O*( P1*P2*nr*fsN1N2Onr + P1*P2*(1-nr)*fsN1N2Or + (1-P1)*(1-P2)*nr*fsI1I2Onr + (1-P1)*(1-P2)*(1-nr)*fsI1I2Or + P1*(1-P2)*nr*fsN1I2Onr + P1*(1-P2)*(1-nr)*fsN1I2Or + (1-P1)*P2*nr*fsI1N2Onr + (1-P1)*P2*(1-nr)*fsI1N2Or ) + (1-O)*( P1*P2*nr*fsN1N2Mnr + P1*P2*(1-nr)*fsN1N2Mr + (1-P1)*(1-P2)*nr*fsI1I2Mnr + (1-P1)*(1-P2)*(1-nr)*fsI1I2Mr + P1*(1-P2)*nr*fsN1I2Mnr + P1*(1-P2)*(1-nr)*fsN1I2Mr + (1-P1)*P2*nr*fsI1N2Mnr + (1-P1)*P2*(1-nr)*fsI1N2Mr )
    return fs


######################################################################################################################################################
#							ANCIENT MIGRATION MODELS - 2 POPULATIONS
######################################################################################################################################################
def AM(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, O = params
    n1, n2 = ns
    """
    Ancient Migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence with migration.
    Ti: Duration of divergence in isolation.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phi = moments.Spectrum(sts)
    phi = moments.Manips.split_1D_to_2D(phi, n1, n2)
    phi.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phi.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsO = phi
    fsM = moments.Numerics.reverse_array(fsO)

    # Sum all spectra
    fs = O*fsO + (1-O)*fsM
    return fs


def AMexp(params, ns):
    nu1e, nu1, nu2e, nu2, m12, m21, Ts, Ti, O = params
    n1, n2 = ns
    """
    Ancient Migration with exponential size changes.

    nu1e: Size of population 1 after split.
    nu1: Size of population 1 after exponential size change.
    nu2e: Size of population 2 after split.
    nu2: Size of population 2 after exponential size change.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence with migration.
    Ti: Duration of divergence in isolation with pop1: exponential change / pop2: exponential change.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phi = moments.Spectrum(sts)
    phi = moments.Manips.split_1D_to_2D(phi, n1, n2)
    phi.integrate([nu1e, nu2e], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    nu1_func = lambda t: nu1e*(nu1/nu1e) ** (t/Ti)
    nu2_func = lambda t: nu2e*(nu2/nu2e) ** (t/Ti)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    phi.integrate(nu_func, Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsO = phi
    fsM = moments.Numerics.reverse_array(fsO)

    # Sum all spectra
    fs = O*fsO + (1-O)*fsM
    return fs


def AM2N(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, nr, bf, O = params
    n1, n2 = ns
    """
    Ancient Migration with heterogeneous Ne (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence with migration.
    Ti: Duration of divergence in isolation.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phinr = moments.Spectrum(sts)
    phinr = moments.Manips.split_1D_to_2D(phinr, n1, n2)
    phinr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phinr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsnrO = phinr
    fsnrM = moments.Numerics.reverse_array(fsnrO)
    
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phir = moments.Spectrum(sts)
    phir = moments.Manips.split_1D_to_2D(phir, n1, n2)
    phir.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phir.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsrO = phir
    fsrM = moments.Numerics.reverse_array(fsrO)

    # Sum all spectra
    fs = O*( nr*fsnrO + (1-nr)*fsrO ) + (1-O)*( nr*fsnrM + (1-nr)*fsrM )
    return fs


def AM2M(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, P, O = params
    n1, n2 = ns
    """
    Ancient Migration with heterogeneous Me (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence with migration.
    Ti: Duration of divergence in isolation.
    P: Fraction of the genome evolving neutrally.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN = moments.Spectrum(sts)
    phiN = moments.Manips.split_1D_to_2D(phiN, n1, n2)
    phiN.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phiN.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsNO = phiN
    fsNM = moments.Numerics.reverse_array(fsNO)

    # Barrier spectrum
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI = moments.Spectrum(sts)
    phiI = moments.Manips.split_1D_to_2D(phiI, n1, n2)
    phiI.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsIO = phiI
    fsIM = moments.Numerics.reverse_array(fsIO)

    # Sum all spectra
    fs = O*( P*fsNO + (1-P)*fsIO ) + (1-O)*( P*fsNM + (1-P)*fsIM )
    return fs


def AM2M2P(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, P1, P2, O = params
    n1, n2 = ns
    """
    Ancient Migration with heterogeneous Me (population-specific).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence with migration.
    Ti: Duration of divergence in isolation.
    P1: Fraction of the genome evolving neutrally in population 1.
    P2: Fraction of the genome evolving neutrally in population 2.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum in population 1 and 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2 = moments.Spectrum(sts)
    phiN1N2 = moments.Manips.split_1D_to_2D(phiN1N2, n1, n2)
    phiN1N2.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phiN1N2.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsN1N2O = phiN1N2
    fsN1N2M = moments.Numerics.reverse_array(fsN1N2O)

    # Barrier spectrum in population 1 and 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2 = moments.Spectrum(sts)
    phiI1I2 = moments.Manips.split_1D_to_2D(phiI1I2, n1, n2)
    phiI1I2.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1I2.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2O = phiI1I2
    fsI1I2M = moments.Numerics.reverse_array(fsI1I2O)

    # Neutral spectrum in population 1, barrier spectrum in population 2
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2 = moments.Spectrum(sts)
    phiN1I2 = moments.Manips.split_1D_to_2D(phiN1I2, n1, n2)
    phiN1I2.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[0, 0]]))
    phiN1I2.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsN1I2O = phiN1I2
    fsN1I2M = moments.Numerics.reverse_array(fsN1I2O)

    # Neutral spectrum in population 2, barrier spectrum in population 1
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2 = moments.Spectrum(sts)
    phiI1N2 = moments.Manips.split_1D_to_2D(phiI1N2, n1, n2)
    phiI1N2.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[m21, 0]]))
    phiI1N2.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1N2O = phiI1N2
    fsI1N2M = moments.Numerics.reverse_array(fsI1N2O)

    # Sum all spectra
    fs = O*( P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O ) + (1-O)*( P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M )
    return fs


def AM2N2M(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, nr, bf, P, O = params
    n1, n2 = ns
    """
    Ancient Migration with heterogeneous Ne (shared between populations) and heterogeneous Me (shared between populations).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence with migration.
    Ti: Duration of divergence in isolation.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    P: Fraction of the genome evolving neutrally.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiNnr = moments.Spectrum(sts)
    phiNnr = moments.Manips.split_1D_to_2D(phiNnr, n1, n2)
    phiNnr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phiNnr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsNOnr = phiNnr
    fsNMnr = moments.Numerics.reverse_array(fsNOnr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiNr = moments.Spectrum(sts)
    phiNr = moments.Manips.split_1D_to_2D(phiNr, n1, n2)
    phiNr.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phiNr.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsNOr = phiNr
    fsNMr = moments.Numerics.reverse_array(fsNOr)

    # Barrier spectrum
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiInr = moments.Spectrum(sts)
    phiInr = moments.Manips.split_1D_to_2D(phiInr, n1, n2)
    phiInr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiInr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsIOnr = phiInr
    fsIMnr = moments.Numerics.reverse_array(fsIOnr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiIr = moments.Spectrum(sts)
    phiIr = moments.Manips.split_1D_to_2D(phiIr, n1, n2)
    phiIr.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiIr.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsIOr = phiIr
    fsIMr = moments.Numerics.reverse_array(fsIOr)

    # Sum all spectra
    fs = O*( P*nr*fsNOnr + P*(1-nr)*fsNOr + (1-P)*nr*fsIOnr + (1-P)*(1-nr)*fsIOr ) + (1-O)*( P*nr*fsNMnr + P*(1-nr)*fsNMr + (1-P)*nr*fsIMnr + (1-P)*(1-nr)*fsIMr )
    return fs


def AM2N2M2P(params, ns):
    nu1, nu2, m12, m21, Ts, Ti, nr, bf, P1, P2, O = params
    n1, n2 = ns
    """
    Ancient Migration with heterogeneous Ne (shared between populations) and heterogeneous Me (population-specific).

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 in neutral regions. In barrier regions, m12=0.
    m21: Migration from pop 1 to pop 2 in neutral regions. In barrier regions, m21=0.
    Ts: Duration of divergence with migration.
    Ti: Duration of divergence in isolation.
    nr: Fraction of non-recombining regions, i.e. with reduced Ne.
    bf: Background factor that defines to which extent Ne is reduced.
    P1: Fraction of the genome evolving neutrally in population 1.
    P2: Fraction of the genome evolving neutrally in population 2.
    O: Fraction of SNPs accurately orientated.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Neutral spectrum in population 1 and 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2nr = moments.Spectrum(sts)
    phiN1N2nr = moments.Manips.split_1D_to_2D(phiN1N2nr, n1, n2)
    phiN1N2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phiN1N2nr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsN1N2Onr = phiN1N2nr
    fsN1N2Mnr = moments.Numerics.reverse_array(fsN1N2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1N2r = moments.Spectrum(sts)
    phiN1N2r = moments.Manips.split_1D_to_2D(phiN1N2r, n1, n2)
    phiN1N2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[m21, 0]]))
    phiN1N2r.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsN1N2Or = phiN1N2r
    fsN1N2Mr = moments.Numerics.reverse_array(fsN1N2Or)

    # Barrier spectrum in population 1 and 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2nr = moments.Spectrum(sts)
    phiI1I2nr = moments.Manips.split_1D_to_2D(phiI1I2nr, n1, n2)
    phiI1I2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1I2nr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2Onr = phiI1I2nr
    fsI1I2Mnr = moments.Numerics.reverse_array(fsI1I2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1I2r = moments.Spectrum(sts)
    phiI1I2r = moments.Manips.split_1D_to_2D(phiI1I2r, n1, n2)
    phiI1I2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[0, 0]]))
    phiI1I2r.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1I2Or = phiI1I2r
    fsI1I2Mr = moments.Numerics.reverse_array(fsI1I2Or)

    # Neutral spectrum in population 1, barrier spectrum in population 2
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2nr = moments.Spectrum(sts)
    phiN1I2nr = moments.Manips.split_1D_to_2D(phiN1I2nr, n1, n2)
    phiN1I2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, m12],[0, 0]]))
    phiN1I2nr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsN1I2Onr = phiN1I2nr
    fsN1I2Mnr = moments.Numerics.reverse_array(fsN1I2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiN1I2r = moments.Spectrum(sts)
    phiN1I2r = moments.Manips.split_1D_to_2D(phiN1I2r, n1, n2)
    phiN1I2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, m12],[0, 0]]))
    phiN1I2r.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsN1I2Or = phiN1I2r
    fsN1I2Mr = moments.Numerics.reverse_array(fsN1I2Or)

    # Neutral spectrum in population 2, barrier spectrum in population 1
    # Spectrum of non-recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2nr = moments.Spectrum(sts)
    phiI1N2nr = moments.Manips.split_1D_to_2D(phiI1N2nr, n1, n2)
    phiI1N2nr.integrate([nu1*bf, nu2*bf], Ts, m=numpy.array([[0, 0],[m21, 0]]))
    phiI1N2nr.integrate([nu1*bf, nu2*bf], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1N2Onr = phiI1N2nr
    fsI1N2Mnr = moments.Numerics.reverse_array(fsI1N2Onr)
    # Spectrum of recombining regions
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2)
    phiI1N2r = moments.Spectrum(sts)
    phiI1N2r = moments.Manips.split_1D_to_2D(phiI1N2r, n1, n2)
    phiI1N2r.integrate([nu1, nu2], Ts, m=numpy.array([[0, 0],[m21, 0]]))
    phiI1N2r.integrate([nu1, nu2], Ti, m=numpy.array([[0, 0],[0, 0]]))
    fsI1N2Or = phiI1N2r
    fsI1N2Mr = moments.Numerics.reverse_array(fsI1N2Or)

    # Sum all spectra
    fs = O*( P1*P2*nr*fsN1N2Onr + P1*P2*(1-nr)*fsN1N2Or + (1-P1)*(1-P2)*nr*fsI1I2Onr + (1-P1)*(1-P2)*(1-nr)*fsI1I2Or + P1*(1-P2)*nr*fsN1I2Onr + P1*(1-P2)*(1-nr)*fsN1I2Or + (1-P1)*P2*nr*fsI1N2Onr + (1-P1)*P2*(1-nr)*fsI1N2Or ) + (1-O)*( P1*P2*nr*fsN1N2Mnr + P1*P2*(1-nr)*fsN1N2Mr + (1-P1)*(1-P2)*nr*fsI1I2Mnr + (1-P1)*(1-P2)*(1-nr)*fsI1I2Mr + P1*(1-P2)*nr*fsN1I2Mnr + P1*(1-P2)*(1-nr)*fsN1I2Mr + (1-P1)*P2*nr*fsI1N2Mnr + (1-P1)*P2*(1-nr)*fsI1N2Mr )
    return fs
