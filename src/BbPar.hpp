#ifndef _BbPar_
#define _BbPar_

class BbPar
{

	private:
		const BbData &data;
	public:
		double alpha;
		Eigen::VectorXd delta;
		Eigen::VectorXd beta;
		double scale;
		double precision;
		double am;
		Eigen::VectorXd dm;
		Eigen::VectorXd bm;
		double sm;
		BbPar(const BbData &data,std::mt19937 &rng) : data(data) {
		  
		  alpha = initialize_scalar(rng);
		  delta = initialize_vec(data.Z.cols(),rng);
		  beta = initialize_vec(data.get_num_beta(),rng);
		  scale = initialize_scalar(rng);

		}

		Eigen::VectorXd get_alpha_vector(){
			return(Eigen::VectorXd::Ones(data.y.size()) * alpha );
		}

		void calculate_precision(){
			
			precision = pow(exp(scale),-2);

		}

};


#endif
