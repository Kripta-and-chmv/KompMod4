import max_value
import tests
import scipy as sc
import random
from IPython.display import display

def integ(v, u):
    a = random.random()
    x = v * (sc.log(-sc.log(a))) + u
    return x

def get_arguments():
    with open("arguments.txt", "r") as f:
        file_str = f.read()
        args = file_str.split(" ")
        v, u, alpha = float(args[0]), float(args[1]), float(args[2])
        return s, p, alpha

u1, u2, u3 = 0, 8, 8
q1, q2, q3 = 1, 4, 0.4

t = max_value.pdf(q1, u1, 8)
print(t)

seq = [x for x in range(100)]

tests.Kolmagorov(seq, q1, u1)

integ(q1, u1)