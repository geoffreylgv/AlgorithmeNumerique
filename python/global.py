#!/usr/bin/python
# -*- coding:utf-8 -*-


"""
    DUPIN Léa - Aéro 2 classe F2
    Ma 223 - Tp 1 : Méthode de Gauss pour la résolution de systèmes linéaires.
    Programme fournissant les graphiques d'analyse avec toutes les méthodes.
    Institut Polytechnique des Sciences Avancées - IPSA Paris
"""

# ---------------- Import des modules nécessaires
import numpy as np
import math
import matplotlib.pyplot as plt
from time import process_time


# ---------------- Partie pivot de Gauss

def ReductionGauss(Aaug):
    # On créé une copie de Aaug
    A = np.copy(Aaug)
    # On récupère la taille de A
    n, m = np.shape(A)
    # On calcule les nouveaux éléments de la matrice A
    for j in range(1, n):
        for i in range(j, n):
            pivot = A[i, j-1] / A[j-1, j-1]
            A[i, :] = A[i, :] - pivot * A[j-1, :]
    # On renvoie la matrice A, matrice réduite de Aaug
    return(A)

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

def Gauss(A, B):
    # On créé Aaug en ajoutant B à A
    stack = np.column_stack([A, B])
    # On réduit la matrice
    reduc = ReductionGauss(stack)
    # On trouve la solution du système
    solution = ResolutionSystTriSup(reduc)
    # On renvoie la matrice solution du système
    return(solution)


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


# ---------------- Partie variantes (pivot partiel et pivot total)
# ---------------- Partie pivot partiel
def PivotPartiel(A, k, n):
    # On récupère la meilleure valeur possible du pivot
    value = k
    for p in range(k, n):
        if abs(A[p][k]) > abs(A[value][k]):
            value = p
    return(value)

def Transposition(A, i, p):
    # On calcule la matrice transposée
    copy = A[i].copy()
    A[i] = A[p]
    A[p] = copy
    return(A)

def ReducGauss_PivotPartiel(Aaug):
    # On créé une copie de Aaug
    A = np.copy(Aaug)
    # On récupère la taille de A
    n, m = np.shape(A)
    # On calcule les nouveaux éléments de la matrice A
    # grâce au meilleur pivot possible
    for k in range(0, n - 1):
        p = PivotPartiel(A, k, n)
        pivot = A[p, k]
        A = Transposition(A, k, p)
        if pivot == 0:
            print("Le pivot est nul")
        else:
            for i in range(k + 1, n):
                G = A[i, k] / pivot
                A[i, :] = A[i, :] - G * A[k, :]
    # On renvoie la matrice A, matrice réduite de Aaug
    return(A)

def Gauss_PivotPartiel(A, B):
    # On créé Aaug en ajoutant B à A
    stack = np.column_stack([A, B])
    # On réduit la matrice
    reduc = ReducGauss_PivotPartiel(stack)
    # On trouve la solution du système
    solution = ResolutionSystTriSup(reduc)
    # On renvoie la matrice solution du système
    return(solution)


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
    plt.subplot(2, 2, 1)

    # Demande de la taille de la matrice maximale à calculer
    taille = int(input("Taille max de la matrice souhaitée ? \n"))

    # ------ Temps de calcul
    print("------------Gauss------------")
    time_list_gauss = []

    for i in range(0, taille + 1, 50):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        Gauss(A, B)
        t = process_time()
        time_list_gauss.append(t)

    # Calculs de temps & création de la liste les contenant tous
    T_list_gauss = []
    for i in range(len(time_list_gauss)):
        if i == len(time_list_gauss)-1:
            T = time_list_gauss[-1] - time_list_gauss[0]
        elif i == len(time_list_gauss):
            None
        else:
            T = time_list_gauss[i + 1] - time_list_gauss[i]
        T_list_gauss.append(T)

    # Affichage des temps en console
    if taille > 50:
        for i in range(0, len(T_list_gauss) -1):
            print("Le temps de calcul pour une matrice de taille ", i*50, "est de :", T_list_gauss[i], "secondes.")
    print("Le temps de calcul pour une matrice de taille ", taille, "est de :", T_list_gauss[-2], "secondes.")

    print("\nLe temps de calcul total est de", T_list_gauss[-1], "secondes")
    minutes = int(T_list_gauss[-1]//60)
    secondes = int(T_list_gauss[-1] % 60)
    if minutes == 1:
        print("Soit environ", minutes, "minute et", secondes, "secondes.")
    elif minutes > 1:
        print("Soit environ", minutes, "minutes et", secondes, "secondes.")  

    # On supprime le temps total afin de pouvoir afficher les temps de calcul
    del(T_list_gauss[- 1])

    print("\n------------LU------------")
    time_list_LU = []
    for i in range(0, taille + 1, 50):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        LU(A, B)
        t = process_time()
        time_list_LU.append(t)
    T_list_LU = []
    for i in range(len(time_list_LU)):
        if i == len(time_list_LU)-1:
            T = time_list_LU[-1] - time_list_LU[0]
        elif i == len(time_list_LU):
            None
        else:
            T = time_list_LU[i + 1] - time_list_LU[i]
        T_list_LU.append(T)
    if taille > 50:
        for i in range(0, len(T_list_LU) -1):
            print("Le temps de calcul pour une matrice de taille ", i*50, "est de :", T_list_LU[i], "secondes.")
    print("Le temps de calcul pour une matrice de taille ", taille, "est de :", T_list_LU[-2], "secondes.")
    print("\nLe temps de calcul total est de", T_list_LU[-1], "secondes")
    minutes = int(T_list_LU[-1]//60)
    secondes = int(T_list_LU[-1] % 60)
    if minutes == 1:
        print("Soit environ", minutes, "minute et", secondes, "secondes.")
    elif minutes > 1:
        print("Soit environ", minutes, "minutes et", secondes, "secondes.") 
    del(T_list_LU[- 1])

    print("\n------------Pivot partiel------------")
    time_list_pivot_partiel = []
    for i in range(0, taille + 1, 50):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        Gauss_PivotPartiel(A, B)
        t = process_time()
        time_list_pivot_partiel.append(t)
    T_list_pivot_partiel = []
    for i in range(len(time_list_pivot_partiel)):
        if i == len(time_list_pivot_partiel)-1:
            T = time_list_pivot_partiel[-1] - time_list_pivot_partiel[0]
        elif i == len(time_list_pivot_partiel):
            None
        else:
            T = time_list_pivot_partiel[i + 1] - time_list_pivot_partiel[i]
        T_list_pivot_partiel.append(T)
    if taille > 50:
        for i in range(0, len(T_list_pivot_partiel) -1):
            print("Le temps de calcul pour une matrice de taille ", i*50, "est de :", T_list_pivot_partiel[i], "secondes.")
    print("Le temps de calcul pour une matrice de taille ", taille, "est de :", T_list_pivot_partiel[-2], "secondes.")
    print("\nLe temps de calcul total est de", T_list_pivot_partiel[-1], "secondes")
    minutes = int(T_list_pivot_partiel[-1]//60)
    secondes = int(T_list_pivot_partiel[-1] % 60)
    if minutes == 1:
        print("Soit environ", minutes, "minute et", secondes, "secondes.")
    elif minutes > 1:
        print("Soit environ", minutes, "minutes et", secondes, "secondes.") 
    del(T_list_pivot_partiel[- 1])

    print("\n------------Pivot total------------")
    """time_list_pivot_total = []
    for i in range(0, taille + 1, 50):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        Gauss_PivotTotal(A, B)
        t = process_time()
        time_list_pivot_total.append(t)
    T_list_pivot_total = []
    for i in range(len(time_list_pivot_total)):
        if i == len(time_list_pivot_total)-1:
            T = time_list_pivot_total[-1] - time_list_pivot_total[0]
        elif i == len(time_list_pivot_total):
            None
        else:
            T = time_list_pivot_total[i + 1] - time_list_pivot_total[i]
        T_list_pivot_total.append(T)
    if taille > 50:
        for i in range(0, len(T_list_pivot_total) -1):
            print("Le temps de calcul pour une matrice de taille ", i*50, "est de :", T_list_pivot_total[i], "secondes.")
    print("Le temps de calcul pour une matrice de taille ", taille, "est de :", T_list_pivot_total[-2], "secondes.")
    print("\nLe temps de calcul total est de", T_list_pivot_total[-1], "secondes")
    minutes = int(T_list_pivot_total[-1]//60)
    secondes = int(T_list_pivot_total[-1] % 60)
    if minutes == 1:
        print("Soit environ", minutes, "minute et", secondes, "secondes.")
    elif minutes > 1:
        print("Soit environ", minutes, "minutes et", secondes, "secondes.")   
    del(T_list_pivot_total[- 1])"""

    print("\n------------Solveur python------------")
    time_list_solveur = []

    for i in range(0, taille + 1, 50):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        np.linalg.solve(A, B)
        t = process_time()
        time_list_solveur.append(t)

    # Calculs de temps & création de la liste les contenant tous
    T_list_solveur = [] 
    for i in range(len(time_list_solveur)):
        if i == len(time_list_solveur)-1:
            T = time_list_solveur[-1] - time_list_solveur[0]
        elif i == len(time_list_solveur):
            None
        else:
            T = time_list_solveur[i + 1] - time_list_solveur[i]
        T_list_solveur.append(T)

    # Affichage des temps en console
    if taille > 50:
        for i in range(0, len(time_list_solveur) -1):
            print("Le temps de calcul pour une matrice de taille ", i*50, "est de :", T_list_solveur[i], "secondes.")

    print("Le temps de calcul pour une matrice de taille ", taille, "est de :", T_list_solveur[-2], "secondes.")

    print("\nLe temps de calcul total est de", T_list_solveur[-1], "secondes")
    minutes = int(T_list_solveur[-1]//60)
    secondes = int(T_list_solveur[-1] % 60)
    if minutes == 1:
        print("Soit environ", minutes, "minute et", secondes, "secondes.")
    elif minutes > 1:
        print("Soit environ", minutes, "minutes et", secondes, "secondes.")  

    # On supprime le temps total afin de pouvoir afficher les temps de calcul
    del(T_list_solveur[- 1])

    abscisse = []
    for i in range(0, taille, 50):
        abscisse.append(i)  

    #  -- Création des courbes : Gauss
    plt.plot(abscisse, T_list_gauss, color = 'b', label = 'Méthode de Gauss')
    #  -- Création des courbes : LU
    plt.plot(abscisse, T_list_LU, color = 'r', label = 'Méthode de LU')
    #  -- Création des courbes : Pivot partiel
    plt.plot(abscisse, T_list_pivot_partiel, color = 'm', label = 'Méthode du pivot partiel')
    #  -- Création des courbes : Pivot total
    # plt.plot(abscisse, T_list_pivot_total, color = 'g', label = 'Méthode du pivot total')
    #  -- Création des courbes : Solveur python
    plt.plot(abscisse, T_list_solveur, color = 'tab:orange', label = 'Solveur python')
    # -- Affichage des courbes
    # Graphique 1 : Temps / taille
    plt.title("Temps de calcul en fonction de la taille de la matrice")
    plt.title("Temps de calcul en fonction de la taille de la matrice")
    plt.ylabel('Temps de calcul en secondes')
    plt.xlabel('Taille de la matrice')
    plt.legend(loc = 'best')
    # Graphique 2 : Temps en échelle logarithmique / taille
    plt.subplot(2, 2, 2)
    plt.plot(abscisse, T_list_gauss, color = 'b', label = 'Méthode de Gauss')
    plt.plot(abscisse, T_list_LU, color = 'r', label = 'Méthode de LU')
    plt.plot(abscisse, T_list_pivot_partiel, color = 'm', label = 'Méthode du pivot partiel')
    # plt.plot(abscisse, T_list_pivot_total, color = 'g', label = 'Méthode du pivot total')
    plt.plot(abscisse, T_list_solveur, color = 'tab:orange', label = 'Solveur python')
    plt.title("Temps de calcul en fonction de la taille de la matrice \n Echelle logarithmique ")
    plt.ylabel('Temps de calcul en secondes (log)')
    plt.xlabel('Taille de la matrice')
    plt.yscale('log')
    plt.legend(loc = 'best')

    # ------ Erreur en fonction de la taille
    # Graphique 3 : Erreur / taille
    plt.subplot(2, 2, 3)
    print("\n------------Gauss------------")
    size = []
    erreur = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        X = Gauss(A, B)
        result = precision(A, X, B)
        size.append(i)
        erreur.append(result)
    #  -- Création des courbes : Gauss
    plt.plot(size, erreur, color = 'b', label = 'Méthode de Gauss')
    plt.subplot(2, 2, 4)
    plt.plot(size, erreur, color = 'b', label = 'Méthode de Gauss')
    plt.subplot(2, 2, 3)

 
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
    #  -- Création des courbes : LU
    plt.plot(size, erreur, color = 'r', label = 'Méthode de LU') 
    plt.subplot(2, 2, 4)
    plt.plot(size, erreur, color = 'r', label = 'Méthode de LU') 
    plt.subplot(2, 2, 3)

    print("\n------------Pivot Total------------")
    """size = []
    erreur = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        X = Gauss_PivotTotal(A, B)
        result = precision(A, X, B)
        size.append(i)
        erreur.append(result)
    #  -- Création des courbes : Pivot total
    # plt.plot(size, erreur, color = 'g', label = 'Méthode du pivot total')"""

    print("\n------------Pivot Partiel------------")
    size = []
    erreur = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        X = Gauss_PivotPartiel(A, B)
        result = precision(A, X, B)
        size.append(i)
        erreur.append(result)
    #  -- Création des courbes : Pivot partiel
    plt.plot(size, erreur, color = 'm', label = 'Méthode du pivot partiel')
    plt.subplot(2, 2, 4)
    plt.plot(size, erreur, color = 'm', label = 'Méthode du pivot partiel')
    plt.subplot(2, 2, 3)

    print("\n------------Solveur python------------")
    size = []
    erreur = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        X = np.linalg.solve(A, B)
        result = precision(A, X, B)
        size.append(i)
        erreur.append(result)
    # Graphique 3 : Erreur / taille
    #  -- Création des courbes : LU
    plt.plot(size, erreur, color = 'tab:orange', label = 'Solveur python')
    plt.subplot(2, 2, 4)
    plt.plot(size, erreur, color = 'tab:orange', label = 'Solveur python')
    plt.subplot(2, 2, 3)

    # -- Affichage des courbes
    plt.xlabel('Taille de la matrice')
    plt.ylabel('Erreur ||A X - B||')
    plt.title('Erreur en fonction de la taille de la matrice')
    plt.grid(True)
    plt.legend(loc = 'best')

    plt.subplot(2, 2, 4)
    plt.xlabel('Taille de la matrice')
    plt.ylabel('Erreur ||A X - B||')
    plt.title('Erreur en fonction de la taille de la matrice')
    plt.yscale('log')
    plt.grid(True)
    plt.legend(loc = 'best')

    plt.suptitle('Courbes TP 1 - Ma223')
    plt.show()

graph()
