#ifndef _TREE_
#define _TREE_

template<class T,class D>
class Tree
{
	private:
		bool UTI_one;
		bool UTI_two;
		double p;
		const D &data;
		double energy_init;
		double energy_prop;
		T left;
		T right;
		T main;
		double n_prime;
		int s_prime;
		double alpha_prime;
		double n_alpha;
    
	public:
		Tree(const D &data,
			 std::mt19937 &rng); 
    
    template<class S>
    void BuildTree(S &model,
                   T &proposed,
                   T &init,
                   double &u,
                   int v,
                   int j,
                   double &epsilon,
                   std::mt19937 &rng);
    
    template<class S>
    void Leapfrog(S &model,
                  T &par,
                  double epsilon);

    int get_s_prime() const;
    
    double get_n_prime() const; 
    
    double get_alpha_prime() const;
    
    double get_n_alpha() const;
    
    const T get_left() const; 
    
    const T get_right() const; 
    
    const T get_main() const; 
    
};

#include "Tree.inl"

#endif
