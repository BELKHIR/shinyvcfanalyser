import numpy as np
import dadi
# Define custom demography models
def prior_onegrow_mig(params, ns, pts):
    """
    Model with growth, split, bottleneck in pop2, exp recovery, migration

    params list is
    nu1F: The ancestral population size after growth. (Its initial size is
          defined to be 1.)
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    m: The scaled migration rate
    Tp: The scaled time between ancestral population growth and the split.
    T: The time between the split and present

    ns = (n1,n2): Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu1F, nu2B, nu2F, m, Tp, T = params
    n1,n2 = ns
    # Define the grid we'll use
    xx = yy = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nu1F)

    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so.
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nu1F, nu2=nu2_func, 
                                    m12=m, m21=m)

    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def SC(params, ns, pts):
    nu1, nu2, m12, m21, Ts, Ti= params
    n1,n2 = ns
    """
    Secondary Contact.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1.
    m21: Migration from pop 1 to pop 2.
    Ts: Duration of divergence in isolation.
    Ti: Duration of divergence with migration.
   
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phi = dadi.Integration.two_pops(phi, xx, T=Ti, nu1=nu1, nu2=nu2, m12=m12, m21=m21)
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return sfs

def SI(params, ns, pts):
    nu1, nu2, Ts = params
    n1,n2 = ns
    """
    Strict Isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: Duration of divergence in isolation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # Equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Split event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # Integrate by setting divergence time to Ts, the population sizes to nu1 and nu2, and the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, T=Ts, nu1=nu1, nu2=nu2, m12=0, m21=0)
    # Calculate the well-oriented spectrum
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))

    # Sum all spectra
    return sfs