import scipy as sc
import random

def pdf(v, u, x):
    """Плотность распределения максимального значения"""
    coef = sc.exp(-(x - u) / v)
    return (coef * sc.exp(-coef)) / v

def cdf(v, u, x):
    """Функция распределения максимального значения"""
    coef = sc.exp(-(x - u) / v)
    return sc.exp(-coef)

def reverse_function_method(v, u):
    """Обратная функция распределения"""
    a = random.random()
    x = u - v * (sc.log(-sc.log(a)))
    return x