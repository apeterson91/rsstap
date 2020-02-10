#ifndef _BbNet_
#define _BbNet_

class BbNet 
{

	private:
		Eigen::VectorXd residual;
		Eigen::VectorXd eta;
	public:
	  const BbData &data;
		BbNet(const BbData &data):
			data(data){};
		

		double calculate_ll(BbPar &par);

		double calculate_energy(BbPar &par);

		double gaussian_ll(BbPar &par);

		void calculate_residual(BbPar &par);

		void calculate_eta(BbPar &par);
};

#include "BbNet.inl"

#endif 
