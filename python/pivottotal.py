#!/usr/bin/python
# -*- coding:utf-8 -*-

"""
    DUPIN Léa - Aéro 2 classe F2
    Ma 223 - Tp 1 : Méthode de Gauss pour la résolution de systèmes linéaires.
    Programme de résolution d'un système linéaire avec la méthode du pivot total.
    Institut Polytechnique des Sciences Avancées - IPSA Paris
"""

# ---------------- Import des modules nécessaires
import numpy as np
import math
import matplotlib.pyplot as plt
from time import process_time


# ---------------- Partie pivot de Gauss

def ResolutionSystTriSup(Taug):
    # On créé une copie de Taug
    A = np.copy(Taug)
    # On récupère la taille de A
    n, m = np.shape(A)
    # On créé une matrice colonne X remplie de 0
    X = np.zeros(n)
    # On calcule les termes solutions de cette matrice X
    for k in range(n - 1, -1, -1):
        S = 0
        for j in range(k + 1, n):
            S = S + A[k, j] * X[j]
        X[k] = (A[k, -1] - S) / A[k, k]
    # On renvoie la matrice solution X
    return(X)


# ---------------- Partie pivot total

def ReducGauss_PivotTotal(Aaug):
    # On créé une copie de Aaug
    A = np.copy(Aaug)
    # On récupère la taille de A
    n, m = np.shape(A)
    # On calcule les nouveaux éléments de la matrice A
    # grâce au meilleur pivot possible
    for k in range(0, n - 1):
        value_max = 0
        for i in range(k, n):
            for j in range(k, n):
                if(abs(Aaug[i][j]) > value_max):
                    value_max = abs(Aaug[i][j])
                    ligne_value_max = i
                    colonne_value_max = j
        for j in range(k, n):
            value = A[k][j]
            A[k][j] = A[ligne_value_max][j]
            A[ligne_value_max][j] = value

        for i in range(k, n):
            value = A[i][k]
            A[i][k] = A[i][colonne_value_max]
            A[i][colonne_value_max] = value
        pivot = A[k, k]
        if (pivot == 0):
            print("Le pivot est nul")

        elif (pivot != 0):
            for i in range(k + 1, n):
                A[i, :] = A[i, :] - (A[i, k] / pivot) * A[k, :]
    # On renvoie la matrice A, matrice réduite de Aaug
    return(A)

def Gauss_PivotTotal(A, B):
    # On créé Aaug en ajoutant B à A
    stack = np.column_stack([A, B])
    # On réduit la matrice
    reduc = ReducGauss_PivotTotal(stack)
    # On trouve la solution du système
    solution = ResolutionSystTriSup(reduc)
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
    print("\n------------Pivot total------------")
    time_list_pivot_total = []

    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        Gauss_PivotTotal(A, B)
        t = process_time()
        time_list_pivot_total.append(t)

    # Calculs de temps & création de la liste les contenant tous
    T_list_pivot_total = []
    for i in range(len(time_list_pivot_total)):
        if i == len(time_list_pivot_total)-1:
            T = time_list_pivot_total[-1] - time_list_pivot_total[0]
        elif i == len(time_list_pivot_total):
            None
        else:
            T = time_list_pivot_total[i + 1] - time_list_pivot_total[i]
        T_list_pivot_total.append(T)

    # Affichage des temps en console
    if taille > 100:
        for i in range(0, len(T_list_pivot_total) - 1, 100):
            print("Le temps de calcul pour une matrice de taille ", i, "est de :", T_list_pivot_total[i], "secondes.")

    print("Le temps de calcul pour une matrice de taille ", taille, "est de :", T_list_pivot_total[-2], "secondes.")

    print("\nLe temps de calcul total est de", T_list_pivot_total[-1], "secondes")
    minutes = int(T_list_pivot_total[-1]//60)
    secondes = int(T_list_pivot_total[-1] % 60)
    if minutes == 1:
        print("Soit environ", minutes, "minute et", secondes, "secondes.")
    elif minutes > 1:
        print("Soit environ", minutes, "minutes et", secondes, "secondes.")

    # On supprime le temps total afin de pouvoir afficher les temps de calcul
    del(T_list_pivot_total[- 1])

    abscisse = []
    for i in range(0, taille):
        abscisse.append(i)

    #  -- Création de la courbe
    plt.plot(abscisse, T_list_pivot_total, color = 'g', label = 'Méthode du pivot total')
    # -- Affichage de la courbe
    # Graphique 1 : Temps / taille    plt.title("Temps de calcul en fonction de la taille de la matrice")
    plt.title("Temps de calcul en fonction de la taille de la matrice")
    plt.ylabel('Temps de calcul en secondes')
    plt.xlabel('Taille de la matrice')
    plt.grid(True)
    plt.legend(loc = 'best')
    # Graphique 2 : Temps en échelle logarithmique / taille
    plt.subplot(3, 1, 2)
    plt.plot(abscisse, T_list_pivot_total, color = 'g', label = 'Méthode du pivot total')
    plt.title("Temps de calcul en fonction de la taille de la matrice \n Echelle logarithmique ")
    plt.ylabel('Temps de calcul en secondes (log)')
    plt.xlabel('Taille de la matrice')
    plt.yscale('log')
    plt.grid(True)
    plt.legend(loc = 'best')

    # ------ Erreur en fonction de la taille
    plt.subplot(3, 1, 3)
    print("\n------------Pivot Total------------")
    size = []
    erreur = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        X = Gauss_PivotTotal(A, B)
        result = precision(A, X, B)
        size.append(i)
        erreur.append(result)
    # Graphique 3 : Erreur / taille
    #  -- Création de la courbe
    plt.plot(size, erreur, color = 'g', label = 'Méthode du pivot total')
    # -- Affichage de la courbe
    plt.xlabel('Taille de la matrice')
    plt.ylabel('Erreur ||A X - B||')
    plt.title('Erreur en fonction de la taille de la matrice')
    plt.grid(True)
    plt.legend(loc = 'best')
    plt.show()

graph()
