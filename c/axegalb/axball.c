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



/**
 * @zerofunction remplissage par zero les éléments dans la résolution détèrminée
*/ 

void zerofuntion(float a[19][19],float b[19],int n){
    int i,j;
    float eps=1e-4;
    for(i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            if (fabs(a[i][j])<eps)
                a[i][j]=0;
        if (fabs(b[i])<eps)
            b[i]=0;
    }
}

/**
 * @showme_matrice afficahe de la matrice
 */

void showme_matrice(float a[19][19],int n){
    int i,j;
    printf("\n\n");
    for (i=0; i<n; i++)
    {
        printf(" [");
        for (j=0; j<n; j++)
        {
            printf(" %.3f ",a[i][j]);
        }
        printf("]\n");
    }
}


/**
 * @comatrice : calcul de la comatrice
 */
void comatrices(float a[19][19],float c[19][19],int i,int j,int n){
    int l,k;
    for(l=0; l<n; l++)
        for(k=0; k<n; k++){
            if ((l<i)&&(k<j))
                c[l][k]=a[l][k];
            if ((l>i)&&(k<j))
                c[l-1][k]=a[l][k];
            if ((l<i)&&(k>j))
                c[l][k-1]=a[l][k];
            if ((l>i)&&(k>j))
                c[l-1][k-1]=a[l][k];
        }
}


/**
 * @det : determinant de matrice
*/
float det(float a[19][19],int n){
    int k,j;
    float c[19][19],s;

    k=n-1;

    if(n==0)
        return(1);

    s=0;
    for(j=0; j<n; j++)
    {
        comatrices(a,c,k,j,n);
        s=s+pow(-1,k+j)*a[k][j]*det(c,k);
    }
    return(s);
}