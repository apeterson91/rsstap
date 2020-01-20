#ifndef _BData_
#define _BData_

// Holds references to Bnet Data components
// dists - dists matrix
// y - vector of outcomes
// Z - matrix of subject constant coefficients
// rpc - regression parameter code vector
// prior_means - self explanatory
// prior_scales - self explanatory
// prior_scales - self explanatory
class BData
{
	public:
		const Eigen::MappedSparseMatrix<double> &D;
		const Eigen::VectorXd &y;
		const Eigen::MatrixXd &Z;
		const Eigen::ArrayXi  &rpc;
		const Eigen::ArrayXi &start_stop;
		const Eigen::ArrayXd &prior_means;
		const Eigen::ArrayXd &prior_scales;
		BData(const Eigen::VectorXd &y,
				const Eigen::MatrixXd &Z,
				const Eigen::MappedSparseMatrix<double> &D,
				const Eigen::ArrayXi &reg_code,
				const Eigen::ArrayXi &start_stop,
				const Eigen::ArrayXd &prior_means,
				const Eigen::ArrayXd &prior_scales) :
			y(y),
			Z(Z),
			D(D),
			rpc(reg_code),
			start_stop(start_stop),
			prior_means(prior_means),
			prior_scales(prior_scales){};

		const int get_num_obs() const {
			return(y.size());
		}

		const int get_outcome_dist() const {
			return(rpc(0));
		}

		const bool use_intercept() const {
			return(rpc(1));
		}

		const int get_p() const {
			return(Z.cols());
		}

		const int get_q() const {
			return(D.cols());
		}

};

#endif 
