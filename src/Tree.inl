template<class T,class D>
Tree<T,D>::Tree(
		const D &data,
		std::mt19937 &rng):

	data(data),
	left(data,rng),
	right(data,rng),
	main(data,rng){

	main.initialize_momenta(rng);
	left.initialize_momenta(rng);
	right.initialize_momenta(rng);

}

template<class T,class D>
template<class S>
void Tree<T,D>::BuildTree(S &model,
			   T &proposed,
			   T &init,
			   double &u,
			   int v,
			   int j,
			   double &epsilon,
			   std::mt19937 &rng){

	if(j == 0){
		energy_init = model.calculate_energy(init);
		Leapfrog(model,proposed,epsilon);
		energy_prop = model.calculate_energy(proposed);
		n_prime = u <= energy_prop ? 1 : 0;
		s_prime = u < (1000 + energy_prop) ? 1 :0;
		left.copy(main);
		right.copy(main);
		alpha_prime = std::min(1.0,exp(energy_prop - energy_init));
		n_alpha = 1.0;
	}else{
		Tree<T,D> subtree(data,rng);
		subtree.BuildTree(model,proposed,init,u,v,j-1,epsilon,rng);
		//TODO refactor to line below with separate function
		s_prime = subtree.get_s_prime();
		n_prime = subtree.get_n_prime();
		n_alpha = subtree.get_n_alpha();
		main.copy(subtree.get_main());
		left.copy(subtree.get_left());
		right.copy(subtree.get_right());
		alpha_prime = subtree.get_alpha_prime();
		// --------------
		if(subtree.get_s_prime() == 1){
			Tree<T,D> subsubtree(data,rng);
			if( v == -1 ){
				subsubtree.BuildTree(model,left,main,u,v,j-1,epsilon,rng);
				left.copy(subsubtree.get_left());
			}else{
				subsubtree.BuildTree(model,right,main,u,v,j-1,epsilon,rng);
				right.copy(subsubtree.get_right());
			}
			// TODO refactor to line below with separate post processing function
			p = subsubtree.get_n_prime() / (subsubtree.get_n_prime() + subtree.get_n_prime());
			std::uniform_real_distribution<double> die(0.0,1.0);
			if(die(rng) <= p)
				main.copy(subsubtree.get_main());
			alpha_prime = subsubtree.get_alpha_prime() + subtree.get_alpha_prime();
			n_alpha = subtree.get_n_alpha() + subsubtree.get_n_alpha();
			UTI_one = left + right; // overloaded operators that calculate U-Turn
			UTI_two = left % right;
			s_prime = (UTI_one && UTI_two ) ? subsubtree.get_s_prime() : 0 ;
			n_prime = subtree.get_n_prime() + subsubtree.get_n_prime();
			// ------------
		}
	}
}

template<class T, class D>
template<class S>
void Tree<T,D>::Leapfrog(S &model,
			  T &par,
			  double epsilon){

	model.calculate_gradient(par);

	main.momenta_leapfrog_other(par,epsilon,model.sg);

	main.momenta_leapfrog_position(par,epsilon);

	model.calculate_gradient(main);

	main.momenta_leapfrog_self(epsilon,model.sg);

}

template<class T, class D>
int Tree<T,D>::get_s_prime() const{
	return(s_prime);
}

template<class T, class D>
double Tree<T,D>::get_n_prime() const{
	return(n_prime);
}

template<class T, class D>
double Tree<T,D>::get_alpha_prime() const{
	return(alpha_prime);
}

template<class T, class D>
double Tree<T,D>::get_n_alpha() const{
	return(n_alpha);
}


template<class T, class D>
const T Tree<T,D>::get_left() const{
	return(left);
}


template<class T, class D>
const T Tree<T,D>::get_right() const{
	return(right);
}


template<class T, class D>
const T Tree<T,D>::get_main() const{
	return(main);
}
