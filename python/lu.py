#!/usr/bin/python
# -*- coding:utf-8 -*-

"""
    DUPIN Léa - Aéro 2 classe F2
    Ma 223 - Tp 1 : Méthode de Gauss pour la résolution de systèmes linéaires.
    Programme de résolution d'un système linéaire avec la méthode de la décomposition LU.
    Institut Polytechnique des Sciences Avancées - IPSA Paris
"""

# ---------------- Import des modules nécessaires
import numpy as np
import math
import matplotlib.pyplot as plt
from time import process_time


# ---------------- Partie décomposition LU
def DecompositionLU(A):
    # On récupère la taille de A
    n, m = np.shape(A)
    # On créé une copie de A
    U = np.copy(A)
    # On créé une matrice identité de taille n
    L = np.eye(n)
    if m != n:
        print("La matrice n'est pas carrée.")
        return(None)
    # On calcule les éléments de L et U
    for j in range(n):
        for i in range(j + 1, n):
            pivot = U[i, j] / U[j, j]
            U[i, :] = U[i, :] - pivot * U[j, :]
            L[i, j] = pivot
    # On renvoie les matrices L et U
    return(L, U)

def ResolutionLU(L, U, B):
    Y = []
    # On récupère la taille de B
    n, m = np.shape(B)
    # On résoud le sytème grâce à la décomposition L U
    for i in range(n):
        Y.append(B[i])
        for j in range(i):
            Y[i] = Y[i] - (L[i, j] * Y[j])
        Y[i] = Y[i] / L[i, i]
    X = np.zeros(n)
    for i in range(n, 0, - 1):
        X[i - 1] = (Y[i - 1] - np.dot(U[i - 1, i:], X[i:])) / U[i - 1, i - 1]
    # On renvoie la matrice solution X
    return(X)

def LU(A, B):
    # On décompose la matrice A en deux matrices L et U
    L, U = DecompositionLU(A)
    # On trouve la solution du système
    solution = ResolutionLU(L, U, B)
    # On renvoie la matrice solution du système
    return(solution)


# ---------------- Calculs de précision
def precision(A, X, B):
    # On récupère la taille de B (nombre de colonnes)
    n = len(B)
    # On la redimensionne
    B = np.reshape(B, (1, n))
    # on récupère les éléments de X et B dans une matrice 1-D
    X = np.ravel(X)
    B = np.ravel(B)
    # On calcule l'erreur
    a = np.dot(A, X) - B
    # On renvoie la norme ||A X - B||
    return(np.linalg.norm(a))


# ---------------- Graphiques
def graph():
    # Mise en page pour mettre les 3 graphiques
    plt.gcf().subplots_adjust(wspace = 0.5, hspace = 0.5)
    plt.subplot(3, 1, 1)

    # Demande de la taille de la matrice maximale à calculer
    taille = int(input("Taille max de la matrice souhaitée ? \n"))

    # ------ Temps de calcul
    print("\n------------LU------------")
    time_list_LU = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        LU(A, B)
        t = process_time()
        time_list_LU.append(t)

    # Calculs de temps & création de la liste les contenant tous
    T_list_LU = []
    for i in range(len(time_list_LU)):
        if i == len(time_list_LU)-1:
            T = time_list_LU[-1] - time_list_LU[0]
        elif i == len(time_list_LU):
            None
        else:
            T = time_list_LU[i + 1] - time_list_LU[i]
        T_list_LU.append(T)

    # Affichage des temps en console
    if taille > 100:
        for i in range(0, len(T_list_LU) - 1, 100):
            print("Le temps de calcul pour une matrice de taille ", i, "est de :", T_list_LU[i], "secondes.")

    print("Le temps de calcul pour une matrice de taille ", taille, "est de :", T_list_LU[-2], "secondes.")

    print("\nLe temps de calcul total est de", T_list_LU[-1], "secondes")
    minutes = int(T_list_LU[-1]//60)
    secondes = int(T_list_LU[-1] % 60)
    if minutes == 1:
        print("Soit environ", minutes, "minute et", secondes, "secondes.")
    elif minutes > 1:
        print("Soit environ", minutes, "minutes et", secondes, "secondes.")

    # On supprime le temps total afin de pouvoir afficher les temps de calcul
    del(T_list_LU[- 1])

    abscisse = []
    for i in range(0, taille):
        abscisse.append(i)

    #  -- Création de la courbe
    plt.plot(abscisse, T_list_LU, color = 'r', label = 'Méthode de LU')

    # -- Affichage de la courbe
    # Graphique 1 : Temps / taille
    plt.title("Temps de calcul en fonction de la taille de la matrice")
    plt.ylabel('Temps de calcul en secondes')
    plt.xlabel('Taille de la matrice')
    plt.grid(True)
    plt.legend(loc = 'best')
    # Graphique 2 : Temps en échelle logarithmique / taille
    plt.subplot(3, 1, 2)
    plt.plot(abscisse, T_list_LU, color = 'r', label = 'Méthode de LU')
    plt.title("Temps de calcul en fonction de la taille de la matrice \n Echelle logarithmique ")
    plt.ylabel('Temps de calcul en secondes (log)')
    plt.xlabel('Taille de la matrice')
    plt.yscale('log')
    plt.grid(True)
    plt.legend(loc = 'best')

    # ------ Erreur en fonction de la taille
    plt.subplot(3, 1, 3)
    print("\n------------LU------------")
    size = []
    erreur = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        X = LU(A, B)
        result = precision(A, X, B)
        size.append(i)
        erreur.append(result)
    # Graphique 3 : Erreur / taille
    #  -- Création des courbes : LU
    plt.plot(size, erreur, color = 'r', label = 'Méthode de LU')
    # -- Affichage des courbes
    plt.xlabel('Taille de la matrice')
    plt.ylabel('Erreur ||A X - B||')
    plt.title('Erreur en fonction de la taille de la matrice')
    plt.grid(True)
    plt.legend(loc = 'best')
    plt.show()

graph()
