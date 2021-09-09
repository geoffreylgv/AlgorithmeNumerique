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


/* @dichotomie la méthode de dichotomieaffiche tout ce qui a été saisie par l'utulisateur en faite le système matriciel
*/
void dichotomie(float a, float b, int n)
{
    double erreur = fabs(n);
    erreur = 1/pow(10,n);
    float m,inf,sup;
    int i=1000, trouve=0, it=0;
    if(!fonction(a) && !fonction(b))
    {
        if (!fonction(a))
            printf("\nLa solution est l'intervalle [%f ; %f]",a,b);
        else
            printf("\nL'équation n'admet pas de solution");
    }
    else
    {
        if(fonction(a)*fonction(b)<0)
        {
            m= (fonction(a) + fonction(b))/2; /// formule de Ck Dichotomie
            inf = a;
            sup=b;
            while(trouve== 0 && i>0 && erreur< (sup-inf))
            {
                it++;

                if (fonction(inf)*fonction(m)<0)
                {
                    sup=m;
                    m= (fonction(inf) + fonction(sup))/2; /// formule de Ck Dichotomie
                    printf("\nL'iteration N° %d => [%fl, %fl] ", it, inf, sup);
                }
                else if (fonction(inf)*fonction(m)>0)
                {
                    inf= m;
                    m=(fonction(inf) + fonction(sup))/2; /// formule de Ck Dichotomie
                    printf("\nL'iteration N° %d => [%fl, %fl] ", it, inf, sup);
                }
                else
                {
                    printf("\nL'iteration N° %d ", it);
                    if (!fonction(inf))
                        printf("\nLa solution est x= %f",inf);
                    if(!fonction(m))
                        printf("\nLa solution est x= %f",m);
                    trouve=1;
                }
                i--;
            }
            if(!trouve)
                printf("\nLa solution :  x = %f ",sup);
        }
        else if (fonction(a)*fonction(b)>0)
        {
            ///printf("\n Pas de solution dans cet intervalle");
            balayage(a,b,n);
        }
        else
        {
            if(fonction(a)==0)
                printf("\nLa solution est x= %.0f",a);
            if (fonction(b)==0)
                printf("\nLa solution est x= %.0f",b);
        }
    }

}

void lagrange(float a, float b, int n)
{
    double erreur;
    n = fabs(n);
    erreur = 1/pow(10,n);
    float m,inf,sup;
    int i=1000, it=0, trouve=0;
    if(!fonction(a) && !fonction(b))
    {
        if (!fonction(a))
            printf("\nL'intervalle [%f ; %f] est solution",a,b);
        else
            printf("\nL'équation n'admet pas de solution");
    }
    else
    {
        if(fonction(a)*fonction(b)<0)
        {
            m= a - (fonction(a)*(b-a)/(fonction(b)-fonction(a))); /// formule de Ck Lagrange
            inf = a;
            sup=b;
            while(trouve== 0 && i>0 && erreur< (sup-inf))
            {
                if (fonction(inf)*fonction(m)<0)
                {
                    sup=m;
                    m=inf - (fonction(inf)*(sup - inf))/(fonction(sup)-fonction(inf));/// formule de Ck Lagrange
                    printf("\nL'iteration N° %d => [%fl, %fl] ", it, inf, sup);
                }
                else if (fonction(inf)*fonction(m)>0)
                {
                    inf= m;
                    m=inf - (fonction(inf)*(sup - inf))/(fonction(sup)-fonction(inf));/// formule de Ck Lagrange
                    printf("\nL'iteration N° %d => [%fl, %fl] ", it, inf, sup);
                }
                else
                {
                    if (!fonction(inf))
                        printf("\nLa solution est x= %f",inf);
                    if(!fonction(m))
                        printf("\nLa solution est x= %f",m);
                    trouve=1;
                }
                i--;
            }
            if(!trouve)
                printf("\nLa solution :  x = %f ",sup);
        }
        else if (fonction(a)*fonction(b)>0)
        {
            balayage(a,b,n);
            ///printf("\n Pas de solution dans cet intervalle");
        }
        else
        {
            if(fonction(a)==0)
                printf("\nLa solution est x= %.0f",a);
            if (fonction(b)==0)
                printf("\nLa solution est x= %.0f",b);
        }
    }
}


