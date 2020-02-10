
double BbNet::calculate_ll(BbPar &par){

	double out = 0;
	switch(data.get_outcome_dist()){

		case 0:
			return(gaussian_ll(par));

	}
	return(out);
}



double BbNet::calculate_energy(BbPar &par){

	double out = calculate_ll(par);
	return(out);

}


double BbNet::gaussian_ll(BbPar &par){

	double out = 0;
	calculate_residual(par);

	out +=  - data.get_num_obs() * par.scale ; 

    out += - .5 * par.precision * (residual).dot(residual);

	return(out);
}

void BbNet::calculate_residual(BbPar &par){


	residual = data.y - eta;

}

void BbNet::calculate_eta(BbPar &par){

	eta = par.get_alpha_vector() + data.Z * par.delta + data.D * data.T * par.beta;
}
