#ifndef _NUTS_
#define _NUTS_

template<class T,class U,class D>
class NUTS{

	private:
		double epsilon;
	    double epsilon_bar;
		double t_naught;
		double p;
		double H_bar;
		double gamma;
		double kappa;
		double log_z;
		double mu_beta;
		double n_alpha;
		double n_alpha_prime;
		double adapt_delta; 
		int n;
		int s;
		int j;
		int vj;
		int treedepth;
		const int chain = 1;
		bool accepted = false;
		bool UTI_one;
		bool UTI_two;
		T model;
		U par;
		U parLeft;
		U parRight;


	public:
		NUTS(
			const D &data,
			std::mt19937 &rng
			 );


		void preTreePrep(
				const int &iter_ix,
				const int &warm_up,
				const int &iter_max,
				std::mt19937 &rng);

		void traverseTree(const D &data,
                      std::mt19937 &rng);

	   void updateEpsilon(int &iter_ix);

	   template<class Z>
	   void store_samples(int &iter_ix,Z &storage);

};

#include "NUTS.inl"

#endif 
