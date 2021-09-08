/**
*@author geoffrey logovi geoffreylogovi2@gmail.com
*@name axball.c method  [Ax=B]
*@description  algorithm for all Matrice Ax=B 
    (the most methods : gauss, gauss partiel, pivot and gauss jordan,
    cholesky, lu crout, lu doolittle, jacobie, gauss-seidel)
*@authorinit web and the arabe code. ;)
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/* @name fillmatrice : on remplie par l'utulisateur la matrice définie
*/
void fillmatrice(float a[19][19],float b[19],int n){
    int i,j;
    printf("\n  La matrice A  de gauche vers la droite:\n\n");
    for (i=0; i<n; i++)
    {
        printf(" ligne %i \n", i+1);
        for (j=0; j<n; j++)
            //printf(" (%i%d) ",[i][j]);
            scanf("%f",&a[i][j]);
    }
    printf("\n * Le vecteur B:\n\n");
    for(i=0; i<n; i++)
    {
        printf(" b%d : ", i+1);
        scanf("%f",&b[i]);
    }
}

/* @shome_iwrote affiche tout ce qui a été saisie par l'utulisateur en faite le système matriciel
*/
void showme_iwrote(float a[19][19],float b[19],int n){
    int i,j;
    printf("\n\n");
    for (i=0; i<n; i++)
    {
        printf(" [");
        for (j=0; j<n; j++)
        {
            printf(" %.2f ",a[i][j]);
        }
        printf("] [ %.2f ]",b[i]);
        printf("\n");
    }
}

/**
 * @rule : vérification
 */
float rule(float x[19],int n){
    float ref;
    int i;
    ref=0;
    for(i=0; i<n; i++)
        if (x[i]>ref)
            ref=x[i];
    return(ref);
}
