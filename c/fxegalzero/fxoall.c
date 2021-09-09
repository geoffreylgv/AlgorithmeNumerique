/**
*@author geoffrey logovi geoffreylogovi2@gmail.com code écrit avec soin
*@name fxoall.c method  [f(x)=0]
*@description  algorithm for all Matrice Ax=B 
    (the most methods : gauss, gauss partiel, pivot and gauss jordan, cholesky, 
    lu crout, lu doolittle, jacobie, gauss-seidel)
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



double fonction(float x)
{
    return 10*x-5;//pow((double)x,2)-2*x + 1;
}
double fonctDeriver(float x)
{
    return 10;//2*x-2;
}

// test function (ici tu ecris ta fonction 
// i.e la fonction que tu veux utiliser pour faire le test)
double fonctionFixe(float x)
{
    return (x*x +1) / 2;
}

// fonction de contrôlll 
float controlEntier(char *message)
{
    float n;
    int retour;
    /// controle de la saisie
    do
    {
        printf("\n%s",message);
        retour= scanf("%f",&n);
        fflush(stdin);
    }
    while(!retour);
    return n;
}

int controlEntierN(char *message)
{
    int n;
    int retour;
    /// controle de la saisie
    do
    {
        printf("\n%s",message);
        retour= scanf("%d",&n);
        fflush(stdin);
    }
    while(!retour);
    return n;
}

