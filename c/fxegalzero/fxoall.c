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


/* @dichotomie la méthode de dichotomie
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


/* @lagrange la méthode de lagrange
*/
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

/* @pointfixe la méthode de pointfixe
*/
void pointFixe(float a, int n)
{
	int i = 10000, trouver = 0;
	int temp = 1;
	double erreur, calcErreur;
	double xn, xn1;
	n = fabs(n);
	erreur = 1 / pow(10, n);
	xn = a;

	while(temp < i && calcErreur > erreur && trouver != 1)
	{
		xn1 = fonctionFixe(xn);
		calcErreur = fabs(xn1 - xn) / fabs(xn1);
		printf("\nItération n°%d Xn+1 = %f ", temp, xn1);
		if(calcErreur < erreur)
		{
			printf("\n\n La convergence est atteinte en %d itérations \n Et le point fixe est Xn+1 = %f", temp, xn1);
			trouver = 1;
		}
		else if(!(temp < i))
		{
			printf("\n La convergence n'est pas atteinte en %d itérations", i);
			trouver = 1;
		}

		xn = xn1;
		temp++;
	}

}

/* @pointFixe la méthode de pointFixe dexiemmmmm (petit probleme)
*/
/*
void pointFixe(float a, int n)
{
    double erreur;
    double temp;
    int i=1000,trouve=0;
    double compteur;
    n = fabs(n);
    erreur = 1/(pow(10,n));
    temp=0;
    while(temp<i && fabs(fonction(a)-a/fonction(a))>erreur && trouve==0)
    {
        compteur=a;
        a = fonction(a);
        if(fonction(compteur)==0)
        {
            printf("\n Le resultat est : %lf",compteur);
            trouve=1;
        }
        temp++;
    }

    if(trouve==0)
    {
        printf("\n Erreur de convergence de la methode.");
    }
}
*/

/* @secante la méthode de secante
*/
void secante(float a, float b, int n)
{
    double erreur;
    n = fabs(n);
    erreur = 1/pow(10,n);
    float t=0,x2= 0, f1=0,f2=0 ;
    int it=1;

    do
    {
        f1=fonction(a);
        f2=fonction(b);
        x2 = b - ((f2*(b-a))/(f2-f1));
        a=b;
        b=x2;
        if(f2<0)
            t=fabs(f2);
        else
            t=f2;
        printf("\nItération N° %d : X = %f \n",it,x2);
        it++;
    }
    while(t>erreur);
    if( isfinite(x2) == 0)
        printf("\nune division par zéro est faite à ce niveau.\n");
    else
        printf("\nLa solution est %f",x2);
}

/* @newton la méthode de newton
*/
void newton(float x0, int n)
{
    double prec;
    n = fabs(n);
    prec = 1/pow(10,n);
    float y, dy;
    double iprec;
    int i=10000, iter = 1;
    if (fonctDeriver(x0))
    {
        y= fonction(x0);
        dy= fonctDeriver(x0);
        iprec= 1/prec;
        while (fabs(y)>prec && dy<iprec && i>0)
        {
            printf("\nItération N°: %d, f(%f) = %f, df(%f) = %f",iter,x0,y,x0,dy);
            x0 -= y/dy;
            y= fonction(x0);
            dy= fonctDeriver(x0);
            i--;
            iter++;
        }

        if(i)
        {
            printf("\nLa solution est X= %f",x0);
        }
        else
            printf("\nLa méthode ne converge pas après %d itérations", iter);
    }
    else
        printf("\nLa méthode ne converge pas car df(%f)= 0",x0);
}


/* @corde1 la méthode de corde1
*/
void corde1(float a, float b, float x0, int n)
{
    double prec;
    n = fabs(n);
    prec = 1/pow(10,n);
    float y, qk;
    float x1 = x0;
    int i=10000, iter = 1;
    if (fonctDeriver(x0))
    {
        y= fonction(x0);
        qk= fonctDeriver(x0);
        while (fabs(y)>prec && i>0)
        {
            printf("\nItération N°: %d, f(%f) = %f, df(%f) = %f",iter,x0,y,x1,qk);
            x0 -= y/qk;
            y= fonction(x0);
            i--;
            iter++;
        }

        if(i)
        {
            printf("\nLa solution est X= %f",x0);
        }
        else
            printf("\nLa méthode ne converge pas après %d itérations", iter);
    }
    else
        printf("\nLa méthode ne converge pas car df(%f)= 0",x0);
}

/* @corde2 la méthode de corde2
*/
void corde2(float x0, int n)
{

    double prec;
    n = fabs(n);
    prec = 1/pow(10,n);
    float y, qk;
    float x1 = x0;
    int i=10000, iter = 1;
    if (fonctDeriver(x0))
    {
        y= fonction(x0);
        qk= fonctDeriver(x0);
        while (fabs(y)>prec && i>0)
        {
            printf("\nItération N°: %d, f(%f) = %f, df(%f) = %f",iter,x0,y,x1,qk);
            x0 -= y/qk;
            y= fonction(x0);
            i--;
            iter++;
        }

        if(i)
        {
            printf("\nLa solution est X= %f",x0);
        }
        else
            printf("\nLa méthode ne converge pas après %d itérations", iter);
    }
    else
        printf("\nLa méthode ne converge pas car df(%f)= 0",x0);
}


/* @lagrange la méthode de lagrange
*/

void balayage(float a, float b, int n)
{
    double erreur = fabs(n);
    erreur = 1/pow(10,n);
    float pas;
    float inf,sup;
    int trouve=0;
    pas=1;
    inf=a;
    sup=b;
    if (fonction(a)==0 && fonction(b)==0)
    {
        printf("\nL'intervalle [%.2f ; %.2f] est solution", a,b);
    }
    if ((fonction(a)*fonction(b)>0)||(fonction(a)*fonction(b)<0))
    {
        while (((pas+inf)<b) && ((sup-inf)>erreur) && (trouve==0))
        {
            if ((fonction(inf)*fonction(inf+pas))<0)
            {
                sup= inf+pas;
                pas/=10;
            }
            else if ((fonction(inf)*fonction(inf+pas))>0)
            {

                inf+=pas;
            }
            else
            {
                if (fonction(inf)== 0)
                    printf("\nLa solution est x= %f",inf);
                if (fonction(inf+pas)==0)
                    printf("\nla solution est x= %f",inf+pas);
                trouve= 1;
            }
        }
    }
    else
    {
        if(fonction(a)==0)
            printf("\nLa solution est x= %.0f",a);
        if (fonction(b)==0)
            printf("\nLa solution est x= %.0f",b);
        trouve= 1;
    }
    if(!trouve)
        printf("\nPas de solution dans cet intervalle");
}


