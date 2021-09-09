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


/* ****************
* Méthode de Gauss
* *****************
*/ 

void gauss(float a[19][19],float b[19],int n){
    float x[19],p,s;
    int i,j,k;

    for(k=0; k<n-1; k++){
        if (a[k][k]==0){
            printf("\n\n Oups pivot nul !! \n");
            printf("La méthode de Gauss ne peut pas être utilisée\n\n");
        }


        for(i=k+1; i<n; i++){//réduisons la matrice
            p=a[i][k]/a[k][k];
            for (j=k; j<n; j++)
                a[i][j]=a[i][j]-p*a[k][j];
            b[i]=b[i]-p*b[k];
        }
    }


    for(i=n-1; i>=0; i--){//résolvons avec le for
        s=0;
        for(j=i+1; j<n; j++)
            s=s+a[i][j]*x[j];
        x[i]=(b[i]-s)/a[i][i];
    }

    zerofuntion(a,b,n);
    printf("\n */* Méthode de Gauss */* \n");
    printf("\n Matrice pécédemment saisie dévient :");
    showme_iwrote(a,b,n);
    printf("\n On a la solution suivante :\n\n");
    for (i=0; i<n; i++)
        printf(" X_%d = %.2f ;\n",i+1,x[i]);
}


/* *****************************
* Méthode de Gauss Pivot partiel
* ******************************
*/
void gauss_partielle(float a[19][19],float b[19],int n){
    float x[19],p,s,ref,temp;
    int i,j,k,ligne;

    for(k=0; k<n-1; k++)
    {
        //*************recherche du maximum pour réaliser le pivot*/
        ref=0;
        for(i=k; i<n; i++){
            if(fabs(a[i][k])>ref){
                ref=fabs(a[i][k]);
                ligne=i;
            }
        }
        //***************pivotation */
        for(j=k; j<n; j++){//for pivotation
            temp=a[k][j];
            a[k][j]=a[ligne][j] ;
            a[ligne][j]=temp;
        }

        temp=b[k];
        b[k]=b[ligne];
        b[ligne]=temp;

        if (a[k][k]==0){
            printf("\n\n Impossible d'appliquer la méthode de Gauss pivot partiel (le pivot n'est pas différent de zero)\n\n");
        }

        /******redct */
        for(i=k+1; i<n; i++){//la matrice devient
            p=a[i][k]/a[k][k];
            for (j=k; j<n; j++)
                a[i][j]=a[i][j]-p*a[k][j];
            b[i]=b[i]-p*b[k];
        }
    }

    /*********re */
    for(i=n-1; i>=0; i--){//résolvons
        s=0;
        for(j=i+1; j<n; j++)
            s=s+a[i][j]*x[j];
        x[i]=(b[i]-s)/a[i][i];
    }
    zerofuntion(a,b,n);
    printf("\n */* Méthode de Gauss avec Stratégie de Pivot */* \n");
    printf("\n Matrice devient :");
    showme_iwrote(a,b,n);
    printf("\n On a la solution suivante :\n\n");
    for (i=0; i<n; i++)
        printf(" X_%d = %f ;\n",i+1,x[i]);
    printf("\n");
}

/* *****************************
* Méthode de Gauss Pivot total
* ******************************
*/
void gauss_total(float a[19][19],float b[19],int n){
    float x[19],p,s,ref,temp;
    int i,j,k,ligne,colonne,pivot_sol[19],temps;

    // vecteur de pivotation des solutions
    for(i=0; i<n; i++)
        pivot_sol[i]=i;

    for(k=0; k<n-1; k++){
        // recherche du maximum (pivot total)
        ref=0;
        for(i=k; i<n; i++){
            for (j=k; j<n; j++) {
                if(fabs(a[i][j])>ref)
                {
                    ref=fabs(a[i][j]);
                    ligne=i;
                    colonne=j;
                }
            }
        }
        // pivotations
        for(j=k; j<n; j++){ //for pivot
            temp=a[k][j];
            a[k][j]=a[ligne][j] ;
            a[ligne][j]=temp;
        }

        temp=b[k];
        b[k]=b[ligne];
        b[ligne]=temp;

        for(i=0; i<n; i++){//résolution
            temp=a[i][k];
            a[i][k]=a[i][colonne] ;
            a[i][colonne]=temp;
        }

        // fillmatrice la fonction qui remplie la matrice : ici elle remplie le vecteur b qui permettra la pivotations
        temps=pivot_sol[k];
        pivot_sol[k]=pivot_sol[colonne];
        pivot_sol[colonne]=temps;

        if (a[k][k]==0){
            printf("\n\n La méthode de Gauss avec stratégie total ne peut pas etre oppérée, pivot non différent de zero \n\n");
        }

        //rd
        for(i=k+1; i<n; i++){//reduisons
            p=a[i][k]/a[k][k];
            for (j=k; j<n; j++)
                a[i][j]=a[i][j]-p*a[k][j];
            b[i]=b[i]-p*b[k];
        }
    }

    //rs
    for(i=n-1; i>=0; i--){//résolution
        s=0;
        for(j=i+1; j<n; j++)
            s=s+a[i][j]*b[j];
        b[i]=(b[i]-s)/a[i][i];
    }

    
    for(i=0; i<n; i++){//pivot total chgm pv
        x[pivot_sol[i]]=b[i];
    }
    
    zerofuntion(a,b,n);
    printf("\n *********Méthode de Gauss avec stratégie de pivot total \n");
    printf("\n la matrice écrite devient :");
    showme_iwrote(a,b,n);
    printf("\n La solution est :\n\n");
    for (i=0; i<n; i++)
        printf(" X_%d = %f ;\n",i+1,x[i]);
    printf("\n");
}

/* ***********************
* Méthode de Gauss Jordan
* ************************
*/
void gauss_jordan(float a[19][19],float b[19],int n){
    float p;
    int i,j,k;

    for(k=0; k<n; k++){
        if (a[k][k]==0){
            printf("\n\n Jordan impossible, pivot non different de zero !!\n\n");
        }

        p=a[k][k];

        //la de normalisation  de la matrice
        for (j=k; j<n; j++)
            a[k][j]=a[k][j]/p;
        b[k]=b[k]/p;

        //rdt
        for(i=0; i<n; i++){//redusons
            if (i!=k){
                p=a[i][k];
                for (j=k; j<n; j++){
                    a[i][j]=a[i][j]-p*a[k][j];
                }
                b[i]=b[i]-p*b[k];
            }
        }
    }

    zerofuntion(a,b,n);
    printf("\n *********** La méthode Gauss Jordan ********** \n");
    printf("\n La matrice devient :");
    showme_iwrote(a,b,n);
    printf("\n solution est :\n\n");
    for(i=0; i<n; i++)
        printf(" X_%d = %f ;\n",i+1,b[i]);
    printf("\n");
}



/* *****************************
* Méthode de LU Doolitlttle
* ******************************
*/
void LU_doolittle(float a[19][19],float b[19],int n){
    float L[19][19],U[19][19],x[19],y[19],s;
    int i,j,k,m;

    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {
            if(i==j)
                L[i][j]=1;
            else
                L[i][j]=0;
            U[i][j]=0;
        }

    for (m=0; m<n; m++){
        for (j=m; j<n; j++){//decomposition U
            s=0;
            for (k=0; k<m; k++)
                s=s+L[m][k]*U[k][j];
            U[m][j]=a[m][j]-s;
        }

        if (U[k][k]==0){
            printf("\n\n Nous ne pouvons pas appliquer LU, diagonal nul \n\n");
        }

        for (i=m+1; i<n; i++){//decomposition L
            s=0;
            for (k=0; k<m; k++)
                s=s+L[i][k]*U[k][m];
            L[i][m]=(a[i][m]-s)/U[m][m];
        }
    }

    // rslv
    for(i=0; i<n; i++){//résolution proprement dite
        s=0;
        for(j=0; j<i; j++)
            s=s+L[i][j]*y[j];
        y[i]=(b[i]-s)/L[i][i];
    }

    for(i=n-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<n; j++)
            s=s+U[i][j]*x[j];
        x[i]=(y[i]-s)/U[i][i];
    }

    printf("\n ***** Méthode Doolittle (LU) ****\n");
    printf("\n Elle s'écrit : A = L * U \n Nous avons les deux décompositions L et U suivante: \n");
    printf("\n Décomposition L :");
    showme_matrice(L,n);
    printf("\n Décomposition U :");
    showme_matrice(U,n);
    printf("\n Solution est :\n\n");
    for (i=0; i<n; i++)
        printf(" X_%d = %f ;\n",i+1,x[i]);
}

/* *****************************
* Méthode de LU Crout
* ******************************
*/

void LU_crout(float a[19][19],float b[19],int n){
    float L[19][19],U[19][19],x[19],y[19],s;
    int i,j,k,m;

    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {
            if(i==j)
                U[i][j]=1;
            else
                U[i][j]=0;
            L[i][j]=0;
        }

    for (m=0; m<n; m++)
    {
        for (i=m; i<n; i++)
        {
            s=0;
            for (k=0; k<m; k++)
                s=s+L[i][k]*U[k][m];
            L[i][m]=a[i][m]-s;
        }

        if (L[k][k]==0)
        {
            printf("\n\n Crout ne peut pas etre utilisée, diagonal nul\n\n");
        }

        for (j=m+1; j<n; j++)
        {
            s=0;
            for (k=0; k<m; k++)
                s=s+L[m][k]*U[k][j];
            U[m][j]=(a[m][j]-s)/L[m][m];
        }
    }

    // resolution par crout
    for(i=0; i<n; i++){
        s=0;
        for(j=0; j<i; j++)
            s=s+L[i][j]*y[j];
        y[i]=(b[i]-s)/L[i][i];
    }

    for(i=n-1; i>=0; i--){
        s=0;
        for(j=i+1; j<n; j++)
            s=s+U[i][j]*x[j];
        x[i]=(y[i]-s)/U[i][i];
    }

    printf("\n ******* La méthode LU Crout : Solution **\n");
    printf("\n Explicite :  A = L * U  décomposition L et U \n");
    printf("\n Décomposition L :");
    showme_matrice(L,n);
    printf("\n Décomposition UU :");
    showme_matrice(U,n);
    printf("\n La solution est :\n\n");
    for (i=0; i<n; i++){
        printf(" X_%d = %f ;\n",i+1,x[i]);
    }
}


/* *****************************
* Méthode de CHOLESKYY :)
* ******************************
*/
void _cholesky(float a[19][19],float b[19],int n){
    float L[19][19],Lt[19][19],x[19],y[19],s,p;
    int i,j,k;

    // contrôl: vérification si M est symétrique
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if (a[i][j]!=a[j][i]){
                printf("\n\n La matrice n'est pas symétrique");
                printf("Impossible donc d'implémenter la méthode de Cholesky \n\n");
            }
        }
    }

    //contrôl : vérification si M est définie positive
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            L[i][j]=0;

    for (i=0; i<n; i++)
    {
        s=0;
        for (k=0; k<i; k++)
            s=s+pow(L[i][k],2);
        p=a[i][i]-s;

        if (p<=0)
        {
            printf("\n\n La matrice n'est pas définie positive ");
            printf("Méthode de Cholesky non implémentable \n\n");
        }

        L[i][i]=sqrt(p);

        for(j=i+1; j<n; j++)
        {
            s=0;
            for (k=0; k<i; k++)
                s=s+L[i][k]*L[j][k];
            L[j][i]=(a[j][i]-s)/L[i][i];
        }
    }

    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            Lt[i][j]=L[j][i];

    // redt
    for(i=0; i<n; i++){//réduisons M
        s=0;
        for(j=0; j<i; j++)
            s=s+L[i][j]*y[j];
        y[i]=(b[i]-s)/L[i][i];
    }

    for(i=n-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<n; j++)
            s=s+Lt[i][j]*x[j];
        x[i]=(y[i]-s)/Lt[i][i];
    }

    printf("\n- ***** Affichage solution avec Cholesky **** \n");
    printf("\n La décomp est la suivante :  A = L * L't \n");
    printf("\n La Décomp L est :");
    showme_matrice(L,n);
    printf("\n La Décomp L't est :");
    showme_matrice(Lt,n);
    printf("\n * La resolution donne :\n\n");
    for (i=0; i<n; i++)
        printf(" X_%d = %f ;\n",i+1,x[i]);
}
void cholesky(float a[19][19], float b[19], int taille)
{
    float L[19][19],Lt[19][19],x[19],y[19],s,p;
    int i,j,k;

// véification de le symérie
    for (i=0; i<taille; i++)
        for (j=0; j<taille; j++)
            if (a[i][j]!=a[j][i])
            {
                printf("\n\n Non symetrique! On ne peut appliquer la methode de Cholesky\n\n");
                exit(-1);
            }

    for (i=0; i<taille; i++)
        for (j=0; j<taille; j++)
            L[i][j]=0;

    for (i=0; i<taille; i++)
    {
        s=0;
        for (k=0; k<i; k++)
            s=s+pow(L[i][k],2);
        p=a[i][i]-s;

        if (p<=0)
        {
            printf("\n\n Non definie positive! On ne peut pas appliquer la methode de Cholesky\n\n");
            exit(-1);
        }

        L[i][i]=sqrt(p);

        for(j=i+1; j<taille; j++)
        {
            s=0;
            for (k=0; k<i; k++)
            {
                s=s+L[i][k]*L[j][k];
            }
            L[j][i]=(a[j][i]-s)/L[i][i];
        }
    }

    for (i=0; i<taille; i++)
        for (j=0; j<taille; j++)
            Lt[i][j]=L[j][i];

/// resolution
    for(i=0; i<taille; i++)
    {
        s=0;
        for(j=0; j<i; j++)
            s=s+L[i][j]*y[j];
        y[i]=(b[i]-s)/L[i][i];
    }

    for(i=taille-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<taille; j++)
            s=s+Lt[i][j]*x[j];
        x[i]=(y[i]-s)/Lt[i][i];
    }

    printf("\n - A = L * Lt \n");
    printf("\n - La matrice L :");
    showme_matrice(L,taille);
    printf("\n - La matrice Lt :");
    showme_matrice(Lt,taille);
    printf("\n - La resolution donne :\n\n");
    for (i=0; i<taille; i++)
        printf(" X_%d = %f ;\n",i+1,x[i]);
}


/* *****************************
* Méthode de Jacobienne
* ******************************
*/
void jacobie(float a[19][19],float b[19],int n){
    float x[19],x1[19],x2[19],s,eps=1e-4;
    int i,j,k,iter=0;

    //Vecteur initial Xo
    printf("\n Innitialison le vecteur solution (Xo) : \n"); 
    printf("\n Saisissez  (Xo) : \n\n");
    for (i=0; i<n; i++){
        printf(" X(0)[%d]= ",i+1);
        scanf("%f",&x1[i]);
    }

    do{
        for(i=0; i<n; i++){
            s=0;
            for (j=0; j<n; j++)
                if (i!=j)
                    s=s+a[i][j]*x1[j];
            x2[i]=(b[i]-s)/a[i][i];
        }
        for (k=0; k<n; k++){
            x[k]=fabs(x1[k]-x2[k]);
            x1[k]=x2[k];
        }

        iter++;
    }while (rule(x,n)>eps) ;

    printf("\n ***********Solution avec la méthode Jacoobienne *************\n");
    printf("\n solution égal :\n\n");
    for (i=0; i<n; i++)
        printf(" X_%d = %f ;\n",i+1,x2[i]);
    printf("\n Indice de précision (espillon) 10^-10. \n");
    printf("\n Itérration:  %d\n",iter);
}

/* *****************************
* Méthode de Gauss-Seidel
* ******************************
*/
void gauss_seidel(float a[19][19],float b[19],int n){
    float x[19],x1[19],x2[19],s,p,eps=1e-4;
    int i,j,k,iter=0;

    //vecteur initial Xo
   //Vecteur initial Xo
    printf("\n Innitialison le vecteur solution (Xo) : \n"); 
    printf("\n Saisissez  (Xo) : \n\n");
    for (i=0; i<n; i++){
        printf(" X(0)[%d]= ",i+1);
        scanf("%f",&x1[i]);
    }

    do{
        for(i=0; i<n; i++){
            s=0;
            p=0;
            for (j=i+1; j<n; j++)
                s=s+a[i][j]*x1[j];
            for (j=0; j<i; j++)
                p=p+a[i][j]*x2[j];
            x2[i]=(b[i]-s-p)/a[i][i];
        }
        for (k=0; k<n; k++){
            x[k]=fabs(x1[k]-x2[k]);
            x1[k]=x2[k];
        }

        iter++;
    }while (rule(x,n)>eps) ;

    printf("\n *******Seidel ** Méthode de Gauss Seidel *********\n");
    printf("\n  La resolution donne :\n\n");
    for (i=0; i<n; i++){
        printf(" X_%d = %f ;\n",i+1,x2[i]);
    }
    //printf("\n Indice de précision 10^-10. \n");
    printf("\n Intérration compté : %d \n", iter);
}


