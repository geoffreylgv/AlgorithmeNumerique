#!/usr/bin/python
# -*- coding:utf-8 -*-


def bisection(f, x_L, x_R, eps, return_x_list=False):
    f_L = f(x_L)
    if f_L*f(x_R) > 0:
        print "Erreur! Impossible de résoudre l'équation dans cet intervalle. Rappel"
        sys.exit(1)
    x_M = float(x_L + x_R)/2.0
    f_M = f(x_M)
    iteration_counter = 1
    if return_x_list:
        x_list = []

    while abs(f_M) > eps:
        if f_L*f_M > 0:   # meme signe
            x_L = x_M
            f_L = f_M
        else:
            x_R = x_M
        x_M = float(x_L + x_R)/2
        f_M = f(x_M)
        iteration_counter += 1
        if return_x_list:
            x_list.append(x_M)
    if return_x_list:
        return x_list, iteration_counter
    else:
        return x_M, iteration_counter

def f(x):
    return x**2 +x - 1

a = -1;   b = 1

solution, no_iterations = bisection(f, a, b, eps=1.0e-10)

print "Nombre d'itérration: %d" % (1 + 2*no_iterations)
print "La Solution est égal à: %f" % (solution)
