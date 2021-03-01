import numpy as np


def gaussian(center, width, height=1.0):
    """Compute gaussian function."""
    w = width / ( 2.0 * np.sqrt( 2.0 * np.log( 2.0 )) )
    def compute_for(pos):
        return height / (w * (2.0 * np.pi) ** 0.5) * (
            np.exp( -(pos - center)**2 /
                     (2.0 * w**2)))
    return compute_for


def lorentzian(center, width, height=1.0):
    """Compute lorentzian function."""
    def compute_for(pos):
        return height / np.pi * (
            0.5 * width / ((pos - center)**2 + (0.5 * width)**2))
    return compute_for


def broaden(positions, amplitudes,
            width=8.0, xmin=0, xmax=4000, dx=1.0, fun=lorentzian):
    """Perform Lorentzian or Gaussian broadening.

    Parameters
    ----------
    positions, amplitudes : lists
        Arrays of positions and intensities.
    width : float
        Full-width at half-maximum.
    xmin, xmax : int
        Extrema for the calculation.
    dx : float
        Resolution of the resulting broadening.
    fun : {gaussian, lorentzian}


    Returns
    -------
    skx, spky : tuple(x-values, y-values)

    """
    nmin, nmax = int(round(xmin/dx)), int(round(xmax/dx))
    spkx = np.array(range(nmin, nmax)) * dx
    spky = np.zeros(len(spkx))
    for pos, intensity in zip(positions, amplitudes):
        broaden = fun(pos, width, intensity)
        spky += broaden(spkx)
    return spkx, spky


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    ar = [[2.0, 4.0, 12.0], [10.0, 4.0, 4.0]]
    spkxl, spkyl = broaden(*ar, width=0.5, xmax=20, dx=0.1, fun=lorentzian)
    spkxg, spkyg = broaden(*ar, width=0.5, xmax=20, dx=0.1, fun=gaussian)
    plt.plot(spkxl, spkyl, spkxg, spkyg)
    plt.show()

