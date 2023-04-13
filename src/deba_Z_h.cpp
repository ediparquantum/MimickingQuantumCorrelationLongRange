#include<iostream>
#include<cmath>
#include<armadillo>
#include<complex>
#include<fstream>
#define QICLIB_DONT_USE_NLOPT
#include "QIClib-master/include/QIClib"
#include<cmath>

using namespace std;
using namespace arma;
using namespace qic;

string makeString(int i) 
{ 
	std::stringstream sste;
	sste << i;
	string s = sste.str();
	return s;    
} 

string makeString(float i) 
{ 
	std::stringstream sste;
	sste << i;
	string s = sste.str();
	return s;    
} 
string makeString(double i) 
{ 
	std::stringstream sste;
	sste << i;
	string s = sste.str();
	return s;    
} 

int main(int argc, char* argv[])
{
	int n=atof(argv[1]);
	double J = 1.0;
	double gamma = atof(argv[2]);
	ofstream outfileMain;
	outfileMain.open(string("xZ_g_norm_")+makeString(gamma)+string("_n")+makeString(n)+string(".dat"), ios_base::app);		
	double alpha_array[6] = {0.2,0.8, 1.2, 1.5, 2.2, 2.5 };
	double h_array[3] = {0.5,1,1.5}; 
		// for(double alpha = 0; alpha<=3; alpha=alpha+0.1){
	for(int i = 0; i<6; i++){
     	double alpha = alpha_array[i];
		// for(double h = 0; h < 3; h=h+0.1){
		for(int j = 0; j<3; j++){
			double h = h_array[j];
			int Z = n-1;
		    cout << "Z_n...."  << Z <<  "  " << alpha << "  " << h << endl;
			ofstream outfile1;
			outfile1.open(string("Er_r_")+makeString(n)+string("/Er_r_a_")+makeString(alpha)+string("_Zn_")+makeString(Z)+string("_h_")+makeString(h)+string("_g_")+makeString(gamma)+string("_n_")+makeString(n)+string(".dat"), ios_base::app);		
			cx_mat A(n,n,fill::zeros);
			cx_mat B(n,n,fill::zeros);
			double rJ = 0;
			for(int i = 1; i <=Z; i++)
			{
				rJ = rJ+J/pow(i,alpha);
			}
			for(int i = 0; i<= n-1; i++)
			{

				for(int k = 0; k<= n-1; k++)
				{

					if(i==k)
					{
						A(i,k) = -h;
					}
					for(int ite = 1; ite <= Z; ite++){

						if(k == i+ite or k+ite==i)
						{

							// A(i,k) = (J)*(1/pow(abs(k-i),alpha));
							A(i,k) = (J)*(1/pow(abs(k-i),alpha))*(1/rJ);

							if(i>k)
							{
								B(i,k) = -(gamma*J)*(1/pow(abs(k-i),alpha))*(1/rJ);
								// B(i,k) = -(gamma*J)*(1/pow(abs(k-i),alpha));
							}

							if(i<k)
							{
								B(i,k) = (gamma*J)*(1/pow(abs(k-i),alpha))*(1/rJ);	
								// B(i,k) = (gamma*J)*(1/pow(abs(k-i),alpha));
							}

						}

					}
				}	
			}
			cx_mat psi = (A+B)*(A-B);
			cx_mat phi = (A-B)*(A+B);
			vec eigval;
			cx_mat eigvec;
			eig_sym(eigval,eigvec, phi);
			vec eigvalpsi;
			cx_mat eigvecpsi;
			eig_sym(eigvalpsi,eigvecpsi,psi);
			    // sqrt(abs(eigvalpsi)).print("JW Energy Diff");
				//(eigvecpsi.col(0)).print();
			cx_mat eigvecphi(n,n,fill::zeros);
			if(sqrt(abs(eigvalpsi(0)))>1e-6){
				eigvecphi = (A-B)*eigvecpsi*diagmat(1/sqrt(eigvalpsi));
			}
			else{
				for(int i = 0; i<n;i++){
					if(sqrt(abs(eigvalpsi(i)))<1e-5){
						eigvecphi.col(i) = eigvec.col(i);
					}
					else{
						eigvecphi.col(i)=(A-B)*eigvecpsi.col(i)*(1/sqrt(abs(eigvalpsi(i))));

					}
				}
			}
			cx_mat G = -(eigvecpsi*eigvecphi.t());
			complex <double> I(0,1);
			mat sigma_z = {{1,0},{0,-1}};
			mat sigma_x = {{0,1},{1,0}};
			cx_mat sigma_y = {{0,-I},{I,0}};
			mat op_z = kron(sigma_z,sigma_z);
			mat op_x = kron(sigma_x,sigma_x);
			cx_mat op_y = kron(sigma_y,sigma_y);
			int count_x=0;
			int l = 0; 
			for(int m = 1; m < n; m++){
				double czz = real(0.25*(G(l,l)*G(m,m)-G(l,m)*G(m,l)));
				cx_mat Gx(abs(l-m),abs(l-m), fill::zeros);
				for(int i = 0; i < abs(l-m); i++){
					for(int j = 0; j < abs(l-m); j++){
						Gx(i,j)=G(i+l,j+l+1);
					}
				}
				double cxx = real(0.25*det(Gx));
				cx_mat Gy(abs(l-m),abs(l-m), fill::zeros);
				for(int i = 0; i < abs(l-m); i++){
					for(int j = 0; j < abs(l-m); j++){
						Gy(i,j)=G(i+l+1,j+l);
					}
				}
				double cyy = real(0.25*det(Gy));
				double mz_l = -0.5*(real(G(l,l)));
				double mz_m = -0.5*(real(G(m,m)));
				mat rho = 0.25*(eye(4,4)+2*mz_m*(kron(eye(2,2),sigma_z))+2*mz_l*(kron(sigma_z,eye(2,2)))+4*cxx*(kron(sigma_x,sigma_x))+4*cyy*(real(kron(sigma_y,sigma_y)))+4*czz*(kron(sigma_z,sigma_z)));
				double entanglement = log_neg(rho,{1});
				outfile1 << l << "  " << m << "  " << "  " << m-l << "  " << entanglement << endl;
				if(entanglement!=0){
					count_x++;
				}
			}
			for(int Z = 1; Z <=n-1; Z++){
				// cout << "Z..." << Z << endl;
				ofstream outfile2;
				outfile2.open(string("Er_r_")+makeString(n)+string("/Er_r_a_")+makeString(alpha)+string("_Z_")+makeString(Z)+string("_h_")+makeString(h)+string("_g_")+makeString(gamma)+string("_n_")+makeString(n)+string(".dat"), ios_base::app);		
				cx_mat A(n,n,fill::zeros);
				cx_mat B(n,n,fill::zeros);
				double rJ = 0;
				for(int i = 1; i <=Z; i++)
				{
					rJ = rJ+J/pow(i,alpha);
				}
				// cout << "normH.," << rJ << endl;
				for(int i = 0; i<= n-1; i++)
				{

					for(int k = 0; k<= n-1; k++)
					{

						if(i==k)
						{
							A(i,k) = -h;
						}
						for(int ite = 1; ite <= Z; ite++){

							if(k == i+ite or k+ite==i)
							{

								A(i,k) = (J)*(1/pow(abs(k-i),alpha))*(1/rJ);

								if(i>k)
								{
									B(i,k) = -(gamma*J)*(1/pow(abs(k-i),alpha))*(1/rJ);
								}

								if(i<k)
								{
									B(i,k) = (gamma*J)*(1/pow(abs(k-i),alpha))*(1/rJ);
								}

							}

						}
					}
				}
				cx_mat psi = (A+B)*(A-B);
				cx_mat phi = (A-B)*(A+B);
				vec eigval;
				cx_mat eigvec;
				eig_sym(eigval,eigvec, phi);
				// real(psi).print("A+B*A-B");
				vec eigvalpsi;
				cx_mat eigvecpsi;
				eig_sym(eigvalpsi,eigvecpsi,psi);
				cx_mat eigvecphi(n,n,fill::zeros);
				if(abs(eigvalpsi(0))>1e-6){
					eigvecphi = (A-B)*eigvecpsi*diagmat(1/sqrt(eigvalpsi));
				}
				else{
					for(int i = 0; i<n;i++){
						if(sqrt(abs(eigvalpsi(i)))<1e-5){
							eigvecphi.col(i) = eigvec.col(i);
						}
						else{
							eigvecphi.col(i)=(A-B)*eigvecpsi.col(i)*(1/sqrt(abs(eigvalpsi(i))));
							
						}
					}
				}

				cx_mat G = -(eigvecpsi*eigvecphi.t());
				complex <double> I(0,1);
				mat sigma_z = {{1,0},{0,-1}};
				mat sigma_x = {{0,1},{1,0}};
				cx_mat sigma_y = {{0,-I},{I,0}};
				mat op_z = kron(sigma_z,sigma_z);
				mat op_x = kron(sigma_x,sigma_x);
				cx_mat op_y = kron(sigma_y,sigma_y);
				int count=0;
				l = 0; 
				for(int m = 1; m < n; m++){
					// cout << l << "..." << m;
					double czz = real(0.25*(G(l,l)*G(m,m)-G(l,m)*G(m,l)));
					cx_mat Gx(abs(l-m),abs(l-m), fill::zeros);
					for(int i = 0; i < abs(l-m); i++){
						for(int j = 0; j < abs(l-m); j++){
							Gx(i,j)=G(i+l,j+l+1);
						}
					}
					double cxx = real(0.25*det(Gx));
					cx_mat Gy(abs(l-m),abs(l-m), fill::zeros);
					for(int i = 0; i < abs(l-m); i++){
						for(int j = 0; j < abs(l-m); j++){
							Gy(i,j)=G(i+l+1,j+l);
						}
					}
					double cyy = real(0.25*det(Gy));
					double mz_l = -0.5*(real(G(l,l)));
					double mz_m = -0.5*(real(G(m,m)));
					mat rho = 0.25*(eye(4,4)+2*mz_m*(kron(eye(2,2),sigma_z))+2*mz_l*(kron(sigma_z,eye(2,2)))+4*cxx*(kron(sigma_x,sigma_x))+4*cyy*(real(kron(sigma_y,sigma_y)))+4*czz*(kron(sigma_z,sigma_z)));
					double entanglement= log_neg(rho,{1});
					outfile2 << l << "  " << m << "  " << "  " << m-l << "  " << entanglement << endl;
					if(entanglement!=0){
						count++;
					}
				}
				if(count >= count_x){
					outfileMain << alpha << "\t" << h << "\t" << Z << " \t" << count_x << "\t" << count << endl;
					break;
				}
			}

		}
		outfileMain << endl;
		
	}
	return 0;
}
