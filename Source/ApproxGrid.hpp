#pragma once

#include <QuatRep.hpp>

//waste of time
class ApproxGrid
{
public:
	int count = 0;
	float gridWidth = 0.25f;
	int sCount = 0;
	float sGridWidth = 0.25f;

	std::vector<QuatRep> cGrid;
	std::vector<int> cGridIndices;
	std::vector<int> cGridLength;
	std::vector<int> cGridStart;

	std::vector<int> GetCoordsFull(QuatRep e); //Gets coordinates for full grid in alg space, only Generate should use this
	int FullCoordIndex(std::vector<int>& indexes);
	int FullCoordIndex(int x, int y, int z, int w);
	int CondensedIndex(int x, int y, int z);
	std::vector<int> GetCoords(QuatRep e); //Gets coords in compressed grid to use in CGridAt
	QuatRep GetApprox(QuatRep e); //basically CGridAt(GetCoords(e)) doesn't create a temp vector
	QuatRep& CGridAt(std::vector<int>& indexes); //Finds point in cGrid using start and length
	void Generate(int countInt, float gridWidthFloat, int sCountInt, float sGridWidthFloat); //Initialize grid

	//At the end of Generate, assign undefined compressed points to the closest defined point, this fixes the edge case problems
};


class LoopApprox
{
public:
	int gridCount = 96;

	std::vector<float> loopAssocMovMag;
	std::vector<QuatRep> loopAssocMov;
	//potentially another vector that describes how aligned tangent spaces compare, e.g. for torus, tangent space of <1,1> is same as <1,-1> but reflected along y axis\
	so when mapping, <1,1>*1.1 should get mapped to <1,-1> + <1,1>*0.1, not <1,-1>*1.1
	std::vector<float> loopAssocRotMag;
	std::vector<QuatRep> loopAssocRot;

	std::vector<float> loopAssocQuotMag;
	std::vector<QuatRep> loopAssocQuot;

	std::vector<QuatRep> _movBase;

	void Generate();
	QuatRep TryLoop(QuatRep localTrans);
	QuatRep TryLoopQuot(QuatRep localTrans);
};