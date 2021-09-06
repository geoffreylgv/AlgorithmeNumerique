#!/usr/bin/python
# -*- coding:utf-8 -*-


def Newton(f, dfdx, x, eps):
    f_value = f(x)
    iteration_counter = 0
    while abs(f_value) > eps and iteration_counter < 100:
        try:
            x = x - float(f_value)/dfdx(x)
        except ZeroDivisionError:
            print "Oups! - derivative zero for x = ", x
            sys.exit(1)     # Abort with error

        f_value = f(x)
        iteration_counter += 1

    # Une sollution, ou iterration depassee
    if abs(f_value) > eps:
        iteration_counter = -1
    return x, iteration_counter

def f(x):
    #return x**2 +x - 1
    return x**3-1

def dfdx(x):
    return 3*x**2

solution, no_iterations = Newton(f, dfdx, x=1000, eps=1.0e-10)

if no_iterations > 0:    # La Solution 
    print "Nombre d'iterration: %d" % (1 + 2*no_iterations)
    print "La Solution est : %f" % (solution)
else:
    print "Aahh ! Pas de solution trouvée"
