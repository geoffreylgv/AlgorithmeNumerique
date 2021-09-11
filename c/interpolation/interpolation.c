#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void menu();

float* AllocatVect(int n);
float **AllocateMat(int n);
void saisir(float *x,float *y,int n);
float *calcul_de_phi(float *x,float *y,int n);
void Lagrange(float *x,float *y,int n,float *coeff_phi);
void Newton(float *x,float *y,int n,float **T);
void cholesky(float **A, float *b, float *x,  int n);
void Moindre_carre(int n,float *x,float *y);

int main(){
	int i,choix,n;
	char reponse='o';

	float *x = (float *) malloc(sizeof(float));
	float *y = (float *) malloc(sizeof(float));
	float **tab;
	
	printf("\t\t\t------------------------------------\n");
	printf("\t\t\tPROGRAMME D'INTERPOLATION LINEAIRE\n");
	printf("\t\t\t------------------------------------\n");
	
	do{
		printf("\n\t Saisir le nombre points qu'on va interpoler:");
		scanf("%d",&n);
		saisir(x,y,n);
	    menu();

		printf("\n\t\t-----------------------------------------------------\n");
		printf("\n\tMethode de Lagrange");
		Lagrange(x,y,n,calcul_de_phi(x,y,n));
		printf("\n\t\t-----------------------------------------------------\n");	
		printf("\n\tMethode de Newton");
		tab=AllocateMat(n);
		Newton(x,y,n,tab);
		printf("\n\t\t-----------------------------------------------------\n");
		printf("\n\tMethode des Moindres carres");
		Moindre_carre(n,x,y);
			
	/*	printf("\n Au revoir !!!");
		exit(-1);*/
		
		printf("\n Refaire le test:");
		scanf("%c",&reponse);
		fflush(stdin);
	}while(reponse!='n' && reponse!='N');
	
	return 0;
}

void menu(){
	printf("1-) Methode de Lagrange\n");
	printf("2-) Methode de Newton\n");
	printf("3-) Methode des moindres carrees\n\n");
	/*printf("Appuyer sur n'importe qu'elle touche pour quitter.");
	printf("\nVotre choix:");*/
}


float* AllocatVect(int n)
{
    float* Vect = (float*) malloc(sizeof(float)*n);
    if(Vect==NULL)
    {
        printf("\n\t\tMEMOIRE INSUFFISANTE!\n");
        exit(EXIT_FAILURE);
    }
    return Vect;
}

float **AllocateMat(int n)
{
    int i;
    float** Mat;
    Mat = malloc(sizeof(float*)*n);
    for(i=0; i<=n; i++)
    {
        Mat[i] = (float*)malloc(sizeof(float)*n);
        if(Mat[i]==NULL)
        {
            printf("\n\t\tMEMOIRE INSUFFISANTE !\n");
            exit(EXIT_FAILURE);
        }
    }
    return Mat;
}


void saisir(float *x,float *y,int n){
	int i;
    for(i=0; i<n; i++)
    {
        printf("\n\t\tEntrer x[%d]valeur :\t",i);
        scanf("%f",&x[i]);
        printf("\n\t\tEntrer y[%d]valeur :\t",i);
        scanf("%f",&y[i]);
        printf("\n");
    }
    //Affichage des données entrés
    printf("\n\t\t Les donnees saisies :\n");
    for(i=0; i<n; i++)
    {
        printf("\n\tx[%i]==> %f || y[%i]==> %f\n", i, x[i], i, y[i]);
    }
}

float *calcul_de_phi(float *x,float *y,int n){
	float *t=(float *) malloc(sizeof(float));
    int i,j;
    for(i=0; i<n; i++)
    {
        t[i]=1;
        for(j=0; j<n; j++)
        {
            if (i!=j)
            {
                t[i] *= (x[i] - x[j]);
            }
        }
        t[i] = y[i]/t[i];
    }
    return t;
}

void Lagrange(float *x,float *y,int n,float *coeff_phi){
	  int i,j;
	  float *coeff=(float *) malloc(sizeof(float));
    for(i=0; i<n; i++)
    {  
    		 printf("\n\n\tPhi[%d](x)==> %f.",i, coeff_phi[i]);
        	 for(j=0; j<n; j++)
        	 {
                if(i!=j)
                {
                printf("( x - %f)", x[j]);
                }
             }
       
        printf("\n\n");
    }
    printf("\t\tp(x)=");
    for (i=0; i<n; i++)
    {
       // if(coeff_phi[i]!=0.00000){
        	 printf("%f", coeff_phi[i]);
             for (j=0; j<n; j++)
             {
                if (i!=j)
                printf("(x -%.3f)",x[j]);
             }
             if (i<n-1)
               printf("+");
	//	}
       
    }
    printf("\n\n");
}

void Newton(float *x,float *y,int n,float **T){
	int i,j;
    i=0;
    while (i!=n)
    {
        T[i][0]=y[i];
        i++;
    }
    for(j=1; j<n; j++)
    {
        for(i=0; i<n-j; i++)
        {
            T[i][j]= ((T[i][j-1]-T[i+1][j-1])/(x[i]-x[i+j]));
        }
    }
    printf("\n\n\tN.0(x) =  %f\n", T[0][0]);
    for(i=1; i<n; i++)
    {
        printf("\n\n\tN.%d(x) = ",i);
        for(j=0; j<i; j++)
        {
            printf("( x - %f)",x[j]);
        }
        printf("\n\n");
    }
    printf("\n\n\tp(x)=");
    for (i=0; i<n; i++)
    {
        printf("%f", T[0][i]);
        for (j=0; j<i; j++)
        {
            printf("(x-%.3f)", x[j]);
        }
        if (i<n-1)
            printf("+");
    }
    printf("\n\n");
}

void cholesky(float **A, float *b, float *x,  int n)
{
    float **L ;         
    float **T ;         
    float *Y ;          
    int i,j,k  ;        
    double aide,p  ;    

    L=AllocateMat(n);
    T=AllocateMat(n);
    Y=AllocatVect(n);

    // véification de la symérie
    for (i=0; i<n ; i++)
    {
        for (j=0; j<n ; j++)
        {
            if (A[i][j]!=A[j][i])
            {

                printf("\n\n\t La matrice n'est pas symetrique donc on ne peut pas appliquer la methode de cholesky!!!!!!!!!\n\t\t");
                system("PAUSE");
            }
        }

    }

    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            L[i][j]=0;

    for (i=0; i<n; i++)
    {
        aide=0;
        for (k=0; k<i; k++) aide=aide+pow(L[i][k],2);
        p=A[i][i]-aide;

        if (p<=0)
        {
            printf("\n\n\tLa matrice n'est pas positive donc on ne peut pas appliquer la methode de cholesky!!!!!!!!!\n\t");
            system("PAUSE");
        }

        L[i][i]=sqrt(p);

        for(j=i+1; j<n; j++)
        {
            aide=0;
            for (k=0; k<i; k++) aide=aide+L[i][k]*L[j][k];
            L[j][i]=(A[j][i]-aide)/L[i][i];
        }
    }

    for (i=0; i<n; i++) for (j=0; j<n; j++) T[i][j]=L[j][i];

// resolution
    for(i=0; i<n; i++)
    {
        aide=0;
        for(j=0; j<i; j++) aide=aide+L[i][j]*Y[j];
        Y[i]=(b[i]-aide)/L[i][i];
    }

    for(i=n-1; i>=0; i--)
    {
        aide=0;
        for(j=i+1; j<n; j++) aide=aide+T[i][j]*x[j];
        if (T[i][i]==0)
        {
            printf("\n\n\t\terreur,division par zero impossible \n\t\t");
            system("pause");
            exit(-1);
        }
        x[i]=(Y[i]-aide)/T[i][i];
    }
    free(L);
    free(T);
    free(Y);
}

void Moindre_carre(int n,float *x,float *y)
{
    int p;               
    int i, j, k;         
    int aide1, aide2;    
    float **P;           
    float *B;          
    float *A;           

    do
    {
        printf("\n\n\t\tLe degres du polynome : n= %i :\n", n);
        scanf("%i", &p);
    }
    while(p>n || p<0);

    P=AllocateMat(p+1);
    B=AllocatVect(p+1);
    A=AllocatVect(p+1);
    aide2=2*p;

    //Calcul de la matrice
    for (i=0 ; i<p+1 ; i++)
    {
        aide1=aide2;
        for (j=0 ; j<p+1 ; j++)
        {
            P[i][j]=0;
            for (k=0; k<n ; k++)
            {
                P[i][j]+=pow(x[k], aide1);
            }
            aide1--;
        }
        aide2--;
    }


    //Calcul du vecteur B
    aide1=p;
    for (i=0 ; i<p+1 ; i++)
    {
        B[i]=0;
        for (k=0 ; k<n ; k++)
        {
            B[i]+=(pow(x[k], aide1)*y[k]);
        }
        aide1--;
    }
    cholesky(P, B, A, p+1);

    printf("\n\n\tLa matrice P : \n\n\n");

    for(i =0 ; i<p+1 ; i++)
    {
        for (j=0 ; j<p+1 ; j++)
        {
            printf("\tX%d   =",i+1);
            printf("\t%.6f",P[i][j]);
            printf("\n");
        }

    }
    for (j=0 ; j<p+1 ; j++)
        {
            printf("X%d  =",i+1);
            printf("\t%.6f",B[j]);
            printf("\n");
        }
    //Affichage de la solution
    aide2=p;
    printf("\n\t Affichage de la fonction P");
    printf("\n\n\tp(x)= ");
    for (i=0; i<p+1 ; i++)
    {
        printf(" %fx^%i ", A[i], aide2);
        aide2--;
        if (i<p)
            printf("+");
    }
    printf("\n");
}

