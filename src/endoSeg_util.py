# Core utility functions for endoSeg


################################################## SUPPORT FUNCTIONS ###################################################

def convert_to_Npy(filename):
    """Converts a CZI file to np array format."""
    with cz.CziFile(filename) as czi:
        imageArrays = czi.asarray()

    return imageArrays


def normalize(x, factor=256):
    """Normalize an image to go from 0 to factor"""
    y = (x - x.min()) / (x.max() - x.min()) * factor
    return y


def invert(x):
    """Invert the constrast of an image."""
    y = - (x - x.max())
    return y


def rolling_average(a, n):
    """Calculates the rolling average of an array a with a step size n."""
    outAv = np.cumsum(a, dtype=float)
    outAv[n:] = out_av[n:] - out_av[:-n]
    return outAv[n-1:] / n


def savitzky_golay(y, windowSize, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques."""
    try:
      windowSize = np.abs(np.int(windowSize))
      order = np.abs(np.int(order))
    except ValueError:
      raise ValueError("windowSize and order have to be of type int")

    if windowSize % 2 != 1 or windowSize < 1:
      raise TypeError("windowSize size must be a positive odd number")

    if windowSize < order + 2:
      raise TypeError("windowSize is too small for the polynomials order")

    orderRange = range(order+1)
    halfWindow = (windowSize -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in orderRange] for k in range(-halfWindow, halfWindow+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:halfWindow+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-halfWindow-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))

    return np.convolve( m[::-1], y, mode='valid')


def smooth_by_convolve(y, kernalSize, window):
    """Calculates a smooth version of a function by convolving with a specified
       kernal type ('flat', 'hanning', 'hamming', 'bartlett', 'blackman'). Kernal
       size should be odd."""
    kern = np.r_[y[kernalSize-1:0:-1],y,y[-2:-kernalSize-1:-1]]
    print(len(kern))
    w=eval('np.'+window+'(kernalSize)')
    print(len(w))
    ySmooth = np.convolve(w/w.sum(), kern, mode='same')
    return ySmooth


def smooth_by_convolve_linear(y, kernalSize):
    """Calculates a smooth version of a function by convolving with a linear
       kernal"""
    kern = np.ones(kernalSize)/kernalSize
    ySmooth = np.convolve(y, kern, mode='same')
    return ySmooth


def median_filter(sigIn):
  """Subtracts the median value of the input signal and sets residual negative values to zero"""

  sigOut = sigIn - np.median(sigIn[sigIn > 0])
  sigOut[sigOut < 0] = 0

  return sigOut


def residual(params, x, y_hat, epsData):
    """Residual for calculating an assymetric loss function."""
    amp = params['amp']
    decay = params['decay']
    a = params['a']
    b = params['b']
    c = params['c']

    model = (amp * np.exp(x * decay)) + a * x ** 2 + b * x + c
    dy = (model - y_hat)
    lossFunc = (dy ** 2) * ((np.sign(dy) + 0.8) ** 4) / epsData

    return lossFunc

