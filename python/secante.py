#!/usr/bin/python
# -*- coding:utf-8 -*-


def secant(f, x0, x1, eps):
    f_x0 = f(x0)
    f_x1 = f(x1)
    iteration_counter = 0
    while abs(f_x1) > eps and iteration_counter < 100:
        try:
            denominator = float(f_x1 - f_x0)/(x1 - x0)
            x = x1 - float(f_x1)/denominator
        except ZeroDivisionError:
            print "Error! - denominator zero for x = ", x
            sys.exit(1)     # sorti du code avec erreur
        x0 = x1
        x1 = x
        f_x0 = f_x1
        f_x1 = f(x1)
        iteration_counter += 1
    # Une solution ou plusieur iterration
    if abs(f_x1) > eps:
        iteration_counter = -1
    return x, iteration_counter

def f(x):
    return x**2 +x -1

x0 = -1;   x1 = 1

solution, no_iterations = secant(f, x0, x1, eps=1.0e-10)

if no_iterations > 0:    # La Solution
    print "Nombre d'iterration: %d" % (2 + no_iterations)
    print "La solution est: %f" % (solution)
else:
    print "Pas de solution!"
