template<class T,class U, class D>
NUTS<T,U,D>::NUTS(
			const D &data,
			std::mt19937 &rng):

			model(data),
			par(data,rng),
			parLeft(data,rng),
			parRight(data,rng)
{

			epsilon = 1.0;
			p = 0;
			epsilon_bar = 1.0;
			H_bar = 0.0;
			gamma = 0.05;
			t_naught = 10;
			kappa = 0.75;
			adapt_delta = .65;
			epsilon = model.FindReasonableEpsilon(par,rng);
			mu_beta = log(10 * epsilon);
}


template<class T,class U, class D>
void NUTS<T,U,D>::preTreePrep(std::mt19937 &rng){

   par.initialize_momenta(rng);
   parLeft.copy(par);
   parRight.copy(par);
   log_z = model.sample_u(par,rng);
   n = 1;
   s = 1;
   j = 0;
   accepted = false;

}

template<class T,class U, class D>
void NUTS<T,U,D>::traverseTree(const D &data,
						std::mt19937 &rng){

	std::uniform_real_distribution<double> runif(0.0,1.0);

	Tree<U,D> tree(data,rng);
	while(s == 1){
		vj = runif(rng) <= .5 ? 1 : -1;
	   if(vj == -1){
		   tree.BuildTree(model,parLeft,par,log_z,vj,j,epsilon,rng);
		   parLeft.copy(tree.get_left());
	   }
	   else{
		    tree.BuildTree(model,parRight,par,log_z,vj,j,epsilon,rng);
		    parRight.copy(tree.get_right());
	   }
	   if(tree.get_s_prime() == 1){
		   p = std::min(1.0, tree.get_n_prime() / n);
		   if(runif(rng)<=0.5)
			   accepted = true;
	   }

	   n = n + tree.get_n_prime();
	   UTI_one = parLeft + parRight;
	   UTI_two = parLeft % parRight;
	   s = (UTI_one && UTI_two) ? tree.get_s_prime() : 0;
	   j++;
	}

	n_alpha = tree.get_n_alpha();
	n_alpha_prime = tree.get_alpha_prime();
	if(accepted)
		par.copy(tree.get_main());

}

template<class T,class U, class D>
void NUTS<T,U,D>::updateEpsilon(int &iter_ix){

	H_bar = (1.0 - 1.0 / (iter_ix + t_naught)) * H_bar +
		(1.0 /(iter_ix + t_naught)) * (adapt_delta - n_alpha_prime / n_alpha);

	epsilon = exp(mu_beta - (sqrt(iter_ix) / gamma) * H_bar);

	epsilon_bar = exp(pow(iter_ix,-kappa) * log(epsilon) + 
		   (1.0 - pow(iter_ix,-kappa)) * log(epsilon_bar));
}


template<class T,class U, class D>
template<class Z>
void NUTS<T,U,D>::store_samples(int &iter_ix, Z &storage){

   storage.treedepth(iter_ix - 1) = j;
   storage.epsilons(iter_ix - 1) = epsilon;
   storage.acceptance(iter_ix - 1) = accepted;


   if(accepted)
	   storage.get_sample(iter_ix,par,model);

}
