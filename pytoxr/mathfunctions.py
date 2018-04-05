def residuals(constants, function, x, y):
    """
    This is a value function (Kostenfunktion) that is used to optimise the fit of the curve to the data.
    It simply calculates the distance between y-value from real data and y-value from the function (sigmoid/sine/etc).
    """
    return y - function(constants, x)

def sine(sine_constants, x):
    ''' Sine equation used to fit sine curve to data.
    f (x) = a*sin (bx + c) + d
    see http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?reg_classic_dr_variable.htm
    '''
    import numpy as np
    a, b, c, d = sine_constants
    y = a * np.sin (b * x + c) + d
    return y

def sine_perfect_helix(sine_constants_cd, x):
    ''' Sine equation which is forced to retain alpha helical periodicity, with 3.6 residues per turn.
    f (x) = a*sin (bx + c) + d
    a = 0.2
    b = 1.745
    Why is a fixed to 0.2?
        This is arbitrary, resulting in a curve that is 0.4 in height
    Why is b fixed to 1.745?
        Since periodicity = 2*np.pi/a, for a perfect alpha helix b = 2*np.pi/3.6 = 1.745
    '''
    import numpy as np
    a = 0.2
    b = 1.745
    c, d = sine_constants_cd
    y = a * np.sin (b * x + c) + d
    return y