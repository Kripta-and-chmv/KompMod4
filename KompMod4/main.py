import math

def max_value_distribution(x, u, q):
    t = (-x - u) / q
    return math.exp(t - math.exp(t))

def inverse_transform_sampling(u, q):
    

u1, u2, u3 = 0, 8, 8
q1, q2, q3 = 1, 4, 0.4