#pragma once

#include <SymExp.hpp>

#include <vector>

class MultTable
{
public:
	
	int count = 0;
	std::vector<float> table;
	
	MultTable() {}
	MultTable(int countInt) { count = countInt; table.resize(count*count*(count+1)); }

	float& At(int i, int j, int o) 
	{ 
		return table[(count*(count+1)*i)+((count+1)*j)+o+1];
	}

	std::vector<float> Mult(std::vector<float> i, std::vector<float> j);

	std::vector<SymExp> Mult(std::vector<SymExp>& i, std::vector<SymExp>& j);
};

struct AlgebraName
{
public:

};