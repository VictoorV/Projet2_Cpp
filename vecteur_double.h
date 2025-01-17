#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

std::ostream& operator<<(std::ostream& out,const std::vector<double>& u) //vecteur double
{
    for(int i=0;i<u.size();i++){out<<u[i]<<" ";}
    std::cout<<" "<<std::endl;
    return out;
}

std::ostream& operator<<(std::ostream& out,const std::vector<bool>& u) //vecteur bool
{
    for(int i=0;i<u.size();i++){out<<u[i]<<" ";}
    std::cout<<" "<<std::endl;
    return out;
}

std::ostream& operator<<(std::ostream& out,const std::vector<int>& u) //vecteur int
{
    for(int i=0;i<u.size();i++){out<<u[i]<<" ";}
    std::cout<<" "<<std::endl;
    return out;
}

std::ostream& operator<<(std::ostream& out,const std::vector<std::vector<double>>& u) //vecteur mutlidimensionnel
{
    for(int i=0;i<u.size();i++){out<<u[i];}
    //std::cout<<" "<<std::endl;
    return out;
}

std::istream& operator>>(std::istream& in,std::vector<double>& u)
{
    std::cout<<"Rentre le vecteur:"<<std::endl;
    for(int i=0;i<u.size();i++){in>>u[i];}
    return in;
}

std::vector<double> operator*(const double& a,const std::vector<double>& u) //multiplication a gauche
{
    std::vector<double> w(u);
    for(int i=0;i<w.size();i++){w[i]=a*u[i];}
    return w;
}

std::vector<double> operator*(const std::vector<double>& u,const double& a) //multiplication a droite
{
    std::vector<double> w(u);
    for(int i=0;i<w.size();i++){w[i]=u[i]*a;}
    return w;
}

std::vector<double> operator+(const std::vector<double>& u,const std::vector<double>& v)
{
    std::vector<double> w(u);
    for(int i=0;i<w.size();i++){w[i]=u[i]+v[i];}
    return w;
}

std::vector<double> operator-(const std::vector<double>& u,const std::vector<double>& v)
{
    std::vector<double> w(u);
    for(int i=0;i<w.size();i++){w[i]=u[i]-v[i];}
    return w;
}

double operator*(const std::vector<double>& u,const std::vector<double>& v)
{
    double S=0;
    for(int i=0;i<u.size();i++){S=S+u[i]*v[i];}
    return S;
}

double Ninf(std::vector<double> u)
{
    double a=-1;
    for(int i=0;i<u.size();i++)
    {
        if(std::abs(u[i])>a)
        {
            a=std::abs(u[i]);
        }
    }
    return a;
}

double L1(std::vector<double> u)
{
    double a=0;
    for(int i=0;i<u.size();i++){a=a+std::abs(u[i]);}
    return a;
}

double L1(std::vector<bool> u) //bool
{
    double a=0;
    for(int i=0;i<u.size();i++){a=a+std::abs(u[i]);}
    return a;
}

double N2(std::vector<double> u)
{
    double a=0;
    for(int i=0;i<u.size();i++){a=a+std::abs(u[i])*std::abs(u[i]);}
    return std::sqrt(a);
}

double L2(std::vector<double> u)
{
    return N2(u)/std::sqrt(u.size());
}

void lecture(std::string nom)
{
    std::ifstream fichier(nom, std::ios::in);  // on ouvre en lecture
 
    if(fichier)  // si l'ouverture a fonctionné
    {
        std::string ligne;
        while(getline(fichier, ligne))  
        {
            std::cout<<ligne<<std::endl;  
        }
        fichier.close();
    }
    else
        std::cerr<<"Impossible d'ouvrir le fichier !"<<std::endl;
}

void ecrire_vecteur(std::string nom, const std::vector<double>& u)
{
    std::ofstream fichier(nom, std::ios::out | std::ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert
    if(fichier)
        {
            for(int i=0;i<u.size();i++){fichier<<u[i]<<" ";}
            fichier.close();
        }
        else
            std::cerr << "Impossible d'ouvrir le fichier !" << std::endl;
}

std::vector<double> lire_vect(std::string nom_fichier) //lecture de vecteur dans un fichier adaptee à l'exercice
{
    std::ifstream file;
    std::vector<double> vect;
    double val;
    file.open(nom_fichier,std::ios::in); // ouverture du fichier
    if(file)
    {
        while (!file.eof())
        {
            file >> val;
            vect.push_back(val);
        }
    }
    else
        std::cerr<<"Impossible d'ouvrir le fichier !"<<std::endl;
    file.close();
    return vect;
}
