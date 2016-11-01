import scipy
import scipy.special
import sys
import numpy
import max_value
from math import sqrt

def kolmagorov(seq, v, u, alpha):
    def calc_d_plus(seq, v, u):
        d = []
        lng = len(seq)
        for i, x in zip(range(1, lng+1), seq):
            el = i/lng - max_value.cdf(v, u, x)
            d.append(el)
        return max(d)

    def calc_d_minus(seq, v, u):
        d = []
        lng = len(seq)
        for i, x in zip(range(1, lng+1), seq):
            el = max_value.cdf(v, u, x) - (i - 1)/lng
            d.append(el)
        return max(d)

    def calc_dn(seq, v, u):
        d_min = calc_d_minus(seq, v, u)
        d_plus = calc_d_plus(seq, v, u)
        return max(d_min, d_plus)

    def calc_s_star(seq, v, u):
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

    print("Тест Колмогорова:")

    seq.sort()
    
    s_star = calc_s_star(seq, v, u)
    print("\tЗначение статистики - {}".format(s_star))

    prob_s = calc_prob_s_grtr_sstr(s_star)
    print("\tP(S* > S) - {}".format(prob_s))

    
    hit = prob_s > alpha
    print("\tРезультат прохождения теста - {}\n".format(hit))
    return hit

def chisqr_test(sequence, alpha, v, u):
    """Тест Хи-квадрат"""
    print("Тест хи квадрат:")

    mod = max(sequence)
    len_seq = len(sequence)
    # разбиваем отрезок от 0 до mod на интервалы
    intervals_amount = int(5 * scipy.log10(len_seq))
    K = intervals_amount
    lngth = mod/K   
    intervals = [x * lngth for x in range(0, K+1)]
    
    #определяем количество попаданий в интервалы
    hits_amount = []    
    for a, b in zip(intervals[:-1], intervals[1:]):
            count = sum([a <= x < b for x in sequence])
            hits_amount.append(count)

    # Вычисляется вероятность попадания слчайной величины в заданные
    # интервалы при равномерном распределении.
    def calc_probs(intervals):
        return [max_value.cdf(v, u, x) - max_value.cdf(v, u, y) for x, y in zip(intervals[1:], intervals[:-1])]

    probabils = calc_probs(intervals)

    # вычисляется статистика
    addition = 0
    for hits, probs in zip(hits_amount, probabils):
        if probs == 0: continue
        addition += (hits / len_seq - probs)**2 / probs

    s_star = len(sequence) * addition
    print("\tЗначение статистики - {}".format(s_star))

    # вычисляется P(S*>S)
    r = intervals_amount - 1
    print("\tКоличество степеней свободы - {}".format(r))

    def integrand(x, r):
        return x ** (r / 2 - 1) * scipy.exp(-x / 2)

    prob_s = scipy.integrate.quad(integrand, s_star, numpy.inf, args = (r))
    prob_s = prob_s[0] / (2 ** (r / 2) * scipy.special.gamma(r / 2))

    print("\tP(S*>S) - {}".format(prob_s))
    print("\tПрохождение теста хи квадрат - {}\n".format(prob_s > alpha))


    
