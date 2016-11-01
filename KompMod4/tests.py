import scipy
import sys
import numpy
import max_value
from math import sqrt

def Kolmagorov(seq, v, u, alpha):
    def calc_s_star(seq, v, u):
        def calc_dn(seq, v, u):
            def calc_d_plus(seq, v, u):
                d = []
                lng = len(seq)
                for i, x in zip(range(1, lng+1), seq):
                    el = i/lng - max_value.distribution(v, u, x)
                    d.append(el)
                return max(d)
            def calc_d_minus(seq, v, u):
                d = []
                lng = len(seq)
                for i, x in zip(range(1, lng+1), seq):
                    el = max_value.distribution(v, u, x) - (i - 1)/lng
                    d.append(el)
                return max(d)
            
            d_min = calc_d_minus(seq, v, u)
            d_plus = calc_d_plus(seq, v, u)
            return max(d_min, d_plus)
        
        dn = calc_dn(seq, v, u)
        lng = len(seq)
        s_star = (6 * lng * dn + 1) / (6 * sqrt(lng))
        return s_star

    def calc_prob_s_grtr_sstr(s_star):
        i = 0
        k = 0
        for i in range(-10000, 10000):
            k += (-1)**i * scipy.exp(-2 * i**2 * s_star**2)
        return 1 - k
    
    seq.sort()
    
    s_star = calc_s_star(seq, v, u)
    p = calc_prob_s_grtr_sstr(s_star)
    
    prob_s_sstar = p > alpha
    return prob_s_sstar