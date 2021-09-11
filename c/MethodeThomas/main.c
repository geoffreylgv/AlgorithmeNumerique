#include <stdio.h>
#include <stdlib.h>

int main()
{

    void saisirMatrice(float a[100][100], int n);
    void affiche(float a[100][100], int n);
    void afficheColumn(float a[100], int n);
    void saisirMatriceSecond(float r[100], int n);
    void methodeThomas(float a[100][100], float r[100], float x[100], int n);



    float a[100][100], x[100], r[100], b[100], c[100], d[100];
    int dim;
    printf("\n\t=============Methode de Thomas===================\n");

    printf("\n\n\tVeuillez saisir la dimension de la matrice : ");
    scanf("%d",&dim);

    //saisirMatrice(a,dim);
    saisirMatriceTridiagonale(a,dim);
    affiche(a,dim);

    printf("\n\tValeur a la diagonale principale \n");
    afficheColumn(d,dim);

    printf("\n\tValeur a la diagonale superieur \n");
    afficheColumn(c,dim-1);

    printf("\n\tValeur a la diagonale inferieur \n");
    afficheColumn(d,dim-1);


    //saisirMatriceSecond(r,dim);

    return 0;
}


void saisirMatrice(float a[100][100], int n){
    int i,j;

    printf("\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            printf("\tA<%d,%d>",i,j);
            scanf("%f",&a[i][j]);
        }
    }
}

void saisirMatriceSecond(float r[100], int n){
    int i;

    printf("\n\tSaisir matrice second : \n");
    for(i=0; i<n; i++){
        printf("\t");
        printf("R<%d>",i);
        scanf("%f",&r[i]);
    }
}

void affiche(float a[100][100], int n){

     int i,j;


    for(i=0; i<n; i++){
        printf("\n");
        for(j=0; j<n; j++){
           printf("\t");
           printf("%.2f",a[i][j]);
        }
    }
    printf("\n\n");
}

void afficheColumn(float a[100], int n){

    int i;

    printf("\n");
    for(i=0; i<n; i++){
        printf("\n\tx%d = ",i+1);
        printf("%.2f",a[i]);
    }
}

void methodeThomas(float a[100][100], float r[100], float x[100], int n){

    int i,j;
}

void saisirMatriceTridiagonale(float a[100][100], int n){

    int i,j;

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i == j){
                printf("\t");
                printf("A<%d,%d>",i,j);
                scanf("%f",&a[i][j]);
            }
            else if(j == i+1){
                printf("\t");
                printf("A<%d,%d>",i,j);
                scanf("%f",&a[i][j]);
            }else if(j == i- 1){
                printf("\t");
                printf("A<%d,%d>",i,j);
                scanf("%f",&a[i][j]);
            }
            else {
                a[i][j] = 0;
            }
        }
    }
}



void  normaliation(float a[100][100], float c[100], float b[100], float d[100], int n){

    int i,j;



    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i == j){
                d[i] = a[i][j];
            }
            else if(j == i+1){
               c[i] = a[i][j];
            }else if(j == i- 1){
                b[i] = a[i][j];
            }
        }
    }


}


//void mainDiagonale(float a[100][100], float d[100][100], int n){
//    int i,j;
//
//    for(i=0; i<n; i++){
//        for(j=0; j<n; j++){
//            if(i == j){
//                d[i][j] = a[i][j];
//            }
//            else {
//                d[i][j] = 0;
//            }
//        }
//    }
//
//}
//
//
//void LowerDiagonale(float a[100][100], float b[100][100], int n){
//
//    int i,j;
//
//    for(i=0; i<n; i++){
//        for(j=0; j<n; j++){
//            if(j == i-1){
//                b[i][j] = a[i][j];
//            }
//            else {
//                b[i][j] = 0;
//            }
//        }
//    }
//
//}
//
//
//void upperDiagonale(float a[100][100], float c[100][100], int n){
//
//    int i,j;
//
//    for(i=0; i<n; i++){
//        for(j=0; j<n; j++){
//            if(j == i+1){
//                c[i][j] = a[i][j];
//            }
//            else {
//                c[i][j] = 0;
//            }
//        }
//    }
//}
