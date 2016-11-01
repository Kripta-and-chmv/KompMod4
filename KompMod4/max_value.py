from scipy import exp as exp

def pdf(v, u, x):
    coef = -((x - u) / v)
    return exp(coef - exp(coef)) / v

def cdf(v, u, x):
    coef = -((x - u) / v)
    return exp(-coef) 