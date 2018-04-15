# coding=utf-8

import numpy as np
from scipy.special import jvp, jv, jnp_zeros, jn_zeros
from const import eps0, mue0, C0

pi = np.pi


#2D cart2pol transformation
def cart2pol(x, y):
    th = np.arctan2(y, x)
    r = np.hypot(x, y)
    return th, r


def calc_te(n, m, freq, a):
    """Calculat mode in circular waveguide
    :param n: TE n
    :param m: TE m
    :param freq: frequency in Hz
    :param a: radius of circle waveguide in m
    :return: Er, Ephi, Hr, Hphi
    """
    wavelength = C0 / freq
    ds = wavelength / 3
    k = 2 * pi / wavelength

    # here I multiply 1.1 only to calculate more space.
    N = 3 * round(1.1 * a * 2 / ds) + 1
    x0 = np.linspace(-1.1*a, 1.1*a, N)
    x, y = np.meshgrid(x0, x0)

    pnm = jnp_zeros(n, m)
    pnm = pnm[m-1]

    kc = pnm / a
    fc = C0 * kc / 2 / pi

    if(freq<fc):
        raise Exception("This mode cannot spread!")
    
    beta = np.sqrt(k**2 - kc**2)
    w = 2 * pi * freq
    epsr = eps0 * 1

    A = 1
    B = 1

    phi, r = cart2pol(x, y)
    r[np.where(r==0)] = ds / 100
    z = 0

    #Transverse Electric field
    Er = -1j * w * mue0 * n / kc / kc / r * (A * np.cos(n*phi) - B * np.sin(n*phi)) * jv(n, kc*r) * np.exp(-1j*beta*z)
    Ephi = 1j * w * mue0 / kc * (A * np.sin(n*phi) + B * np.cos(n*phi)) * jvp(n, kc*r) * np.exp(-1j*beta*z)

    #Transverse Magnetic field
    Hr = -1j * beta / kc * (A * np.sin(n*phi) + B * np.cos(n*phi)) * jvp(n, kc*r) * np.exp(-1j*beta*z)
    Hphi = -1j * beta * n / kc / kc / r * (A * np.cos(n*phi) - B * np.sin(n*phi)) * jv(n, kc*r) * np.exp(-1j*beta*z)

    E_abs = np.hypot(np.abs(Er), np.abs(Ephi))
    E_abs[np.where(r>a)] = 1e-20

    # filename = 'TE{n},{m}.dat'.format(n=n, m=m)
    # E_abs.tofile(filename)
    print("Calculation Complete!")
    return Er, Ephi, Hr, Hphi


def calc_tm(n, m, freq, a):
    """Calculate mode in circular waveguide
    :param n: TE n
    :param m: TE m
    :param freq: frequency in Hz
    :param a: radius of circle waveguide in m
    :return: Er, Ephi, Hr, Hphi
    """
    wavelength = C0 / freq
    ds = wavelength / 3
    k = 2 * pi / wavelength

    N = 8 * round(1.1 * a * 2 / ds) + 1
    x0 = np.linspace(-1.1*a, 1.1*a, N)
    x, y = np.meshgrid(x0, x0)

    pnm = jn_zeros(n, m)
    pnm = pnm[m-1]

    kc = pnm / a
    fc = C0 * kc / 2 / pi

    if freq < fc:
        raise Exception("This mode cannot spread!")
    
    beta = np.sqrt(k**2 - kc**2)
    w = 2 * pi * freq
    epsr = eps0 * 1

    A = 1
    B = 1

    phi, r = cart2pol(x, y)
    r[np.where(r==0)] = ds / 100
    z = 0

    #Transverse Electric field
    Er = -1j * beta / kc * (A * np.sin(n*phi) + B * np.cos(n*phi)) * jvp(n, kc*r) * np.exp(-1j*beta*z)
    Ephi = 1j * beta * n / kc / kc / r * (A * np.cos(n*phi) - B * np.sin(n*phi)) * jv(n, kc*r) * np.exp(-1j*beta*z)

    #Transverse Magnetic field
    Hr = -1j * beta / kc * (A * np.sin(n*phi) + B * np.cos(n*phi)) * jv(n, kc*r) * np.exp(-1j*beta*z)
    Hphi = -1j * beta * n / kc / kc / r * (A * np.cos(n*phi) - B * np.sin(n*phi)) * jvp(n, kc*r) * np.exp(-1j*beta*z)

    E_abs = np.hypot(np.abs(Er), np.abs(Ephi))
    E_abs[np.where(r>a)] = 1e-20

    # filename = 'TM{n},{m}.dat'.format(n=n, m=m)
    # E_abs.tofile(filename)
    print("Calculation Complete!")
    return Er, Ephi, Hr, Hphi
