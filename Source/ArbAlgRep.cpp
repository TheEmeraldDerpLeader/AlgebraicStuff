#include "ArbAlgRep.hpp"

std::vector<float> MultTable::Mult(std::vector<float> i, std::vector<float> j)
{
	if (i.size() != count+1 || j.size() != count+1)
		return std::vector<float>();
	std::vector<float> hold; hold.resize(count+1);
	for (int k = 0; k < count+1; k++)
		hold[k] = (i[0]*j[k]) + (i[k]*j[0]);

	for (int i1 = 1; i1 < count+1; i1++)
	{
		float ci = i[i1];
		for (int j1 = 1; j1 < count+1; j1++)
		{
			float cj = j[j1];
			for (int k = 0; k < count+1; k++)
			{
				hold[k] += At(i1-1, j1-1, k)*ci*cj;
			}
		}
	}

	return hold;
}

std::vector<SymExp> MultTable::Mult(std::vector<SymExp>& i, std::vector<SymExp>& j)
{
	if (i.size() != count+1 || j.size() != count+1)
		return std::vector<SymExp>();
	std::vector<SymExp> hold; hold.resize(count+1);

	hold[0] = i[0]*j[0];
	for (int k = 1; k < count+1; k++)
		hold[k] = (i[0]*j[k]) + (i[k]*j[0]);

	for (int i1 = 1; i1 < count+1; i1++)
	{
		SymExp& ci = i[i1];
		for (int j1 = 1; j1 < count+1; j1++)
		{
			SymExp& cj = j[j1];
			for (int k = 0; k < count+1; k++)
			{
				hold[k] += (ci*cj)*Product(At(i1-1, j1-1, k-1));
			}
		}
	}

	return hold;
}