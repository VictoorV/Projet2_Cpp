#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <ctime>
#include <cassert>
#include "matrice.h"
#include <random>
using namespace std;

class matricebande : public matrice
{
private:
    vector<int> indice;
public:
    matricebande() {};                         // constructeur par defaut
    matricebande(int,int);                     // constructeur en donnant 2 dimensions
    matricebande(const matricebande&);         // constructeur par recopie
    ~matricebande() {};                        // destructeur par defaut
    
    vector<int> ind() const {return this->indice;}    // retourne le champ privé indice
    vector<int> ind(vector<int>);                     // remplit le champ privé indice

    matricebande laplacien(int);                      // assemblage du laplacien TP3, question (1.1)
    vector<double> operator*(const vector<double>&);
    vector<double> steepest_descent(vector<double>&, vector<double>&, double&);
    vector<double> gradient_pas_fixe(vector<double>&, vector<double>&, double&, double&);

    friend vector<double> smbr(int,double,double);
};

//------------------------------------------------------
// Constructeur
matricebande::matricebande(int n, int d) : matrice(n,d)
{this->indice.resize(d);} // remplit le champ privé indice avec 0
//------------------------------------------------------
// Constructeur par recopie
matricebande::matricebande(const matricebande& A) : matrice(A)
{this->indice=A.ind();}
//-------------------------------------------------------
// Remplit le champ privé indice avec le vecteur v
vector<int> matricebande::ind(vector<int> v)
{return this->indice=v;} 
//-------------------------------------------------------
// Construit le stockage bande de la matrice du laplacien 2D
matricebande matricebande::laplacien(int n){
    
    int N=(n+2)*n;
    double h = 1./(n+1), hh=pow(h,2);
    vector<double> x(n+2), y(n);
    for (int i=0; i<n+2; i++)
        x[i]=i*h;
    for (int i=1; i<n+1; i++)
        y[i]=i*h;

    vector<int> v;
    v.push_back(-n-2);v.push_back(-1); v.push_back(0);
    
    matricebande A(N,3);
    A.ind(v);
        
    int ligne;
 
    // Ici j==0
    for (int i=0; i<n+2; i++){ // CL de Neumann sur les bords verticaux
        ligne=i;
        if (i==0)
            A(ligne,2)= 4./hh;
        else if (i<n+1){
            A(ligne,1)=-1./hh;
            A(ligne,2)= 4./hh;
        }
        else{
            A(ligne,1)=-2./hh;
            A(ligne,2)= 4./hh;
        }
    }
    // CL de Dirichlet sur les bords horizontaux, j==0 déjà traité
    for (int j=1; j<n; j++){
        for (int i=0; i<n+2; i++){ // CL de Neumann sur les bords verticaux
            ligne=(n+2)*j+i;
            if (i==0){
                A(ligne,0)=-1./hh;
                A(ligne,2)= 4./hh;
            }
            else if (i<n+1){
                A(ligne,0)=-1./hh;
                A(ligne,1)=-1./hh;
                A(ligne,2)= 4./hh;
            }
            else{
                A(ligne,0)=-1./hh;
                A(ligne,1)=-2./hh;
                A(ligne,2)= 4./hh;
            }
        }
    }
    //cout << "matrice du laplacien avant symetrisation \n";
    //cout << A << endl;

    // Symetrisation de A
    for (int j=0; j<n; j++){
        ligne=(n+2)*j;
        A[ligne]=A[ligne]*0.5; // modification de toute la ligne
        ligne=(n+2)*j+n+1;
        A[ligne]=A[ligne]*0.5;
    }
    
    return A;
}
//-------------------------------------------------------
// Surcharge l'operateur * 
vector<double> matricebande::operator*(const vector<double>& x)
{
    int n=this->dim1();
    int p=this->ind().size();
    assert(x.size()==n);
    vector<double> y(n,0);
    int j;
    // Partie triangulaire inférieure
    for(int k=0;k<p;k++)
    {
        for(int i=0;i<n;i++){
            j=i+this->ind()[k];      
            if(j>=0)
            {
                y[i]=y[i]+(*this)(i,k)*x[j];
            }
        }
    }
    // ¨Partie triangulaire supérieure
    for(int k=0;k<p-1;k++){
        for(int i=0;i<n;i++){
            j=i+this->ind()[k];
            if(j>=0)
            {
                y[j]=y[j]+(*this)(i,k)*x[i];
            }
        }
    }
    return y;
}
//-------------------------------------------------------
// Algorithme steepest descent
vector<double> matricebande::steepest_descent(vector<double>& b,vector<double>& x,double& tol)
{
    assert((b.size()==this->dim1())&&(x.size()==this->dim1()));
    int stop=10*x.size(),k=0;
    double alpha;
    vector<double> y;
    vector<double> r=b-(*this)*x;
    vector<double> r0(r);
    vector<double> res; // res stocke les ||rk||/||b|| 
    res.push_back(L2(r)/L2(b));
    while((L2(r)>L2(r0)*tol)&&(k<stop))
    {
        k=k+1;
        y=(*this)*r;
        alpha=pow(L2(r),2)/(y*r); // produit scalaire au denominateur
        x=x+r*alpha;
        r=r-y*alpha;
        res.push_back(L2(r)/L2(b));
        if(L2(r)>(L2(b)/tol))
        {
            cout<<"OVERFLOW"<<endl;
            cout<<"Nombre d'iterations : "<<k<<endl;
            break;
        }
    }
    cout<<"Nombre d'iterations : "<<k<<endl;
    return res;
}
//-------------------------------------------------------
// Remplissage du vecteur b
vector<double> smbr(int n,double Ti0,double TiN1)
{
    double h=1./(n+1);
    vector <double> b(n*(n+2),0); // Taille d'apres (1.1)
    // Les premiers et derniers N+1 elements sont non nuls 
    for(int i=0;i<=n+1;i++)
    {
        b[i]=Ti0;
    }

    for(int i=n*(n+1)-2;i<n*(n+2);i++)
    {
        b[i]=TiN1;
    }
    // Symetrisation 
    b[n*(n+2)-1]=b[n*(n+2)-1]/2;  
    b[n*(n+1)-2]=b[n*(n+1)-2]/2;
    b[n+1]=b[n-1]/2;
    b[0]=b[0]/2;
    return b;
}
//-------------------------------------------------------
// Algorithme gradient pas fixe
vector<double> matricebande::gradient_pas_fixe(vector<double>& b,vector<double>& x,double& alpha,double& tol)
{
    assert((b.size()==this->dim1())&&(x.size()==this->dim1()));
    int stop=10*x.size(),k=0;
    vector<double> y;
    vector<double> r=b-(*this)*x;
    vector<double> r0(r); 
    vector<double> res; // res stocke les ||rk||/||b|| 
    while((L2(r)>L2(r0)*tol)&&(k<stop))
    {
        k=k+1;
        y=(*this)*r;
        x=x+r*alpha;
        r=r-y*alpha;
        res.push_back(L2(r)/L2(b));
        if(L2(r)>(L2(b)/tol))
        {
            cout<<"OVERFLOW"<<endl;
            cout<<"Nombre d'iterations : "<<k<<endl;
            break;
        }
    }
    cout<<"Nombre d'iterations : "<<k<<endl;
    return res;
}

double diag(matricebande A) // Pour trouver le plus grand element sur la diagonale de A
{
    double s=-1;
    int n=A.dim1();
    int p=A.dim2();
    for(int i=0;i<n;i++)
    {
        if(s<A(i,p-1))
        {
            s=A(i,p-1);
        }
    }
    return 1./s;
}

int main()
{
    cout<<"Question 1.5 :"<<endl;
    matricebande L;
    L=L.laplacien(3);
    
    cout << "DIM L = " << L.dim1() <<" "<< L.dim2() << endl;
    cout << "Laplacien 2D " <<endl<< L << endl;

    cout<<"Question 1.8 :"<<endl;
    cout<<"N=15"<<endl;
    vector<double> X(15,1);
    cout<<"A*x = "<<L*X<<endl;

    cout<<"Question 2.1 :"<<endl;
    double tol=pow(10,-4);
    double alpha1=diag(L)/5;
    double alpha2=2*diag(L)/5;
    double alpha3=3*diag(L)/5;
    double alpha4=4*diag(L)/5;  
    double alpha5=5*diag(L)/5;
    double alpha6=6*diag(L)/5;

    vector<double> b=smbr(3,100,50);
    cout<<"b="<<b<<endl;
    vector<double> x0(15,0);
    vector<double> x1(15,0);
    vector<double> x2(15,0);
    vector<double> x3(15,0);
    vector<double> x4(15,0);
    vector<double> x5(15,0);
    vector<double> x6(15,0);

    vector<double> res0=L.steepest_descent(b,x0,tol); // Steepest descent 
    cout<<"res0="<<endl;
    cout<<res0<<endl;
    cout<<"x0="<<endl;
    cout<<x0<<endl;

    cout<<"Question 2.2 :"<<endl;
    vector<double> res1=L.gradient_pas_fixe(b,x1,alpha1,tol); // j=1
    cout<<"res1="<<endl;
    cout<<res1<<endl;
    cout<<"x1="<<endl;
    cout<<x1<<endl;
    vector<double> res2=L.gradient_pas_fixe(b,x2,alpha2,tol); // j=2
    cout<<"res2="<<endl;
    cout<<res2<<endl;
    cout<<"x2="<<endl;
    cout<<x2<<endl;
    vector<double> res3=L.gradient_pas_fixe(b,x3,alpha3,tol); // j=3
    cout<<"res3="<<endl;
    cout<<res3<<endl;
    cout<<"x3="<<endl;
    cout<<x3<<endl;
    vector<double> res4=L.gradient_pas_fixe(b,x4,alpha4,tol); // j=4
    cout<<"res4="<<endl;
    cout<<res4<<endl;
    cout<<"x4="<<endl;
    cout<<x4<<endl;
    vector<double> res5=L.gradient_pas_fixe(b,x5,alpha5,tol); // j=5
    cout<<"res5="<<endl;
    cout<<res5<<endl;
    cout<<"x5="<<endl;
    cout<<x5<<endl;
    vector<double> res6=L.gradient_pas_fixe(b,x6,alpha6,tol); // j=6, pas de convergence 
    cout<<"res6="<<endl;
    cout<<res6<<endl;
    cout<<"x6="<<endl;
    cout<<x6<<endl;
	return 0;
}
