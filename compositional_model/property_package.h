#pragma once
#include <fstream>
#include <vector>
#include "adetl/scalars/ADscalar.hpp"
#include <iostream>
#include "dlib/matrix.h"

class Flash_Calculator
{
	typedef dlib::matrix<double, 0, 1> col_vec;
	typedef adetl::ADscalar<> scalar;
public:
	Flash_Calculator(int np,double R_factor,double temp=160+459.67):num_Comp(np),R(R_factor){
		OmegaA = 0.45724;//constant for Peng-Robinson EOS
		OmegaB = 0.07780;
		Z.set_size(np);
		MW.set_size(np);
		Pc.set_size(np);
		Tc.set_size(np);
		W.set_size(np);
		aTc.set_size(np);
		bTc.set_size(np);
		m.set_size(np);
		expwils.set_size(np); 
		load_Comp_data();
		aTc = dlib::pointwise_multiply(OmegaA*pow(R, 2)*pow(Tc, 2), 1 / Pc);
		bTc = dlib::pointwise_multiply(OmegaB*R*Tc, 1 / Pc);
		m=0.37464 + 1.54226*W- 0.26992*pow(W, 2);
		set_T(temp);
		P_adj = 10;
		n_miter = 50;
		x0 = 0.5;// initial guess for V/F ratio;
		toler = 1e-12; //error for Newton method;
	}


	void set_T(double temperature){
		T = temperature;
	}
	void flash_cal(const double &P){
		double pi=WilsonAdj(P);
		ki = dlib::pointwise_multiply(Pc / pi, exp(dlib::pointwise_multiply(5.37*(1+W),(1-Tc/T))));
		bool flag = true;
		double ncom_f = 1 / double(num_Comp);
		fl = dlib::ones_matrix<double>(1,num_Comp);
		fv = dlib::ones_matrix<double>(1, num_Comp);
		while (ncom_f*sum(pow(log(dlib::pointwise_multiply(fv,1/fl)),2))>toler || flag)
		{
			flag = false;

		}

	}
	
private:
	double WilsonAdj(const double &P){	
		double pi;
		expwils = exp(5.37*dlib::pointwise_multiply((1 + W),(1 - Tc / T)));//for isotherm reservoir
		col_vec tep = 1 / dlib::pointwise_multiply(Pc, expwils);
		PminW = 1 / sum(dlib::pointwise_multiply(Z,1 / dlib::pointwise_multiply(Pc, expwils)));
		PmaxW = sum(dlib::pointwise_multiply(Z, dlib::pointwise_multiply(expwils,Pc)));
		if (P>PminW && P<PmaxW){
			pi = P;
		}
		else if (P>=PmaxW){
			pi = PmaxW - P_adj;
		}
		else if (P<=PminW){
			pi = PminW + P_adj;
		}
		return pi;
	}
	double WilsonK(){
		if (dlib::min(ki)<1 && dlib::max(ki)>1){
			std::cout << "No problem with Ki" << std::endl;
		}
		else{
			std::cout << "Problem with Ki" << std::endl;
		}
		col_vec asymp = -1 / (ki-1);
		if (sum(asymp > 1)>0 || sum(asymp>0))
		{
			std::cout << "Take care of possible error" << std::endl;
		}
	}



	void load_Comp_data()
	{	double z,mw,w,tc,pc,tmp;
		std::ifstream input("input/Comp_inform.dat");
		std::ifstream input1("input/Inter_inform.dat");
		inter_cor.set_size(num_Comp,num_Comp);
		for (int i = 0; i < num_Comp;i++)
		{
			input >> z;
			input >> pc;
			input >> tc;
			input >> mw;
			input >> w;
			Z(i)=z;
			MW(i) = mw;
			Tc(i) = tc;
			Pc(i) = pc;
			W(i) = w;
			for (int j = 0; j < num_Comp; j++)
			{
				input1 >> tmp;
				inter_cor(i,j)=tmp;
			}
		}
		input.close();
		input1.close();
	}

private:
	col_vec Z, MW, Tc, Pc, W;
	col_vec aTc, bTc, m;
	col_vec expwils,ki;
	dlib::matrix<double> inter_cor;
	const int num_Comp;
	double OmegaA, OmegaB,R,T;
	double PminW, PmaxW,P_adj;
	double x0, toler;
	int n_miter;
	col_vec fl, fv;
};