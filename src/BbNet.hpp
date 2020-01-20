#ifndef _BbNet_
#define _BbNet_

class BbNet
{

	private:
		const BData &data;
	public:
		BbNet(const BData &data):
			data(data){};

		template<class S>
		double calculate_energy(S &par);
};

#include "BbNet.inl"

#endif 