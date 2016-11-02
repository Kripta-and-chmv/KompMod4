import max_value
import numpy as np
import time
import tests
import scipy as sc
from IPython.display import display
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys

def get_arguments():
    with open("arguments.txt", "r") as f:
        file_str = f.read()
        args = file_str.split(" ")
        v, u, alpha = float(args[0]), float(args[1]), float(args[2])
        return v, u, alpha

def show_prob_density_and_function(v, u):

    absc = [x for x in range(-10, 11)]
    ordin_pdf = [max_value.pdf(v, u, x) for x in absc]
    ordin_cdf = [max_value.cdf(v, u, x) for x in absc]

    plt.plot(absc, ordin_pdf, label='pdf')
    plt.plot(absc, ordin_cdf, 'r', label='cdf')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)

    plt.title('PDF and CDF')
    plt.xticks(absc)
    plt.show()

def testing_seq(v, u, alpha, N):
    def emperic_graph(seq, v, u):
        seq.sort() 
        teor_y = [max_value.pdf(v, u, x) for x in seq]
        yvals=np.arange(len(seq)) / float(len(seq)) 
        plt.plot(seq, yvals) 
        #plt.plot(seq, teor_y, 'r') 

        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Emperical denisty function')
        plt.show()
            
    print("Размер выборки: {}".format(N))
    print("Парметры выборки: v = {}, u = {}".format(v, u))
    start = time.time()
    seq = [max_value.reverse_function_method(v, u) for x in range(N)]
    building_time = time.time() - start

    emperic_graph(seq, v, u)
    with open("max_val{}_{}_{}.txt".format(len(seq), v, u), "w") as f:
        f.write(str(seq))

    print("Время генерации выборки: {}".format(building_time))
    
    tests.chisqr_test(seq, alpha, v, u)
    tests.kolmagorov(seq, v, u, alpha)


def main():
    sys.stdout = open("output.txt", "w+")

    v, u, alpha = get_arguments()
    testing_seq(v, u, alpha, 50)
    testing_seq(v, u, alpha, 200)
    testing_seq(v, u, alpha, 1000)

    show_prob_density_and_function(v, u)

if __name__ == "__main__":
    main()