#include "ApproxGrid.hpp"

std::vector<int> ApproxGrid::GetCoordsFull(QuatRep e)
{
	std::vector<int> hold; hold.resize(4);
	int axisCellCount = count/2;
	e = (e+QuatRep((sGridWidth/2)+((sCount/2)*sGridWidth),(gridWidth/2)+(axisCellCount*gridWidth),(gridWidth/2)+(axisCellCount*gridWidth),(gridWidth/2)+(axisCellCount*gridWidth))); //Move into only positive coords to unify floor rounding
	e.s /= sGridWidth; e.i /= gridWidth; e.j /= gridWidth; e.k /= gridWidth;
	hold[0] = int(e.s); hold[1] = int(e.i); hold[2] = int(e.j); hold[3] = int(e.k); //round towards closest to 0
    return hold;
}

int ApproxGrid::FullCoordIndex(std::vector<int>& indexes)
{
	return (count*count*count*indexes[0])+(count*count*indexes[1])+(count*indexes[2])+indexes[3];
}

int ApproxGrid::FullCoordIndex(int x, int y, int z, int w)
{
	return (count*count*count*x)+(count*count*y)+(count*z)+w;
}

int ApproxGrid::CondensedIndex(int x, int y, int z)
{
	return (count*count*x)+(count*y)+z;
}

std::vector<int> ApproxGrid::GetCoords(QuatRep e)
{
    return std::vector<int>();
}

QuatRep ApproxGrid::GetApprox(QuatRep e)
{
	thread_local std::vector<int> indexes;
	indexes = GetCoordsFull(e);
	if (indexes[1] < 0) indexes[1] = 0; if (indexes[1] >= count) indexes[1] = count-1;
	if (indexes[2] < 0) indexes[2] = 0; if (indexes[2] >= count) indexes[2] = count-1;
	if (indexes[3] < 0) indexes[3] = 0; if (indexes[3] >= count) indexes[3] = count-1;
	int condIndex = CondensedIndex(indexes[1], indexes[2], indexes[3]);

	if (cGridLength[condIndex] == 0)
		return QuatRep(0, 0, 0, 0);

	int sIndex = indexes[0]-cGridStart[condIndex];
	if (sIndex < 0) sIndex = 0; if (sIndex >= cGridLength[condIndex]) sIndex = cGridLength[condIndex];
    return cGrid[cGridIndices[condIndex]+sIndex];
}

QuatRep& ApproxGrid::CGridAt(std::vector<int>& indexes)
{
	int i = (count*count*indexes[0])+(count*indexes[1])+indexes[2];
	int w = indexes[3]-cGridStart[i];
	if (w < 0 || w >= cGridLength[i])
		abort();
	return cGrid[cGridIndices[i]+w];
}

void ApproxGrid::Generate(int countInt, float gridWidthFloat, int sCountInt, float sGridWidthFloat)
{
	std::vector<QuatRep> fullGrid;
	gridWidth = gridWidthFloat;
	count = countInt;
	sCount = sCountInt;
	sGridWidth = sGridWidthFloat;

	if (count % 2 == 0)
		count++;
	
	fullGrid.resize(sCount*count*count*count);
	for (int i = 0; i < fullGrid.size(); i++)
		fullGrid[i] = QuatRep(1,0,0,0); //scaler defined is uninitialized
	std::vector<float> gridDists; gridDists.resize(fullGrid.size()); //also check magnitude to prevent far points that loop
	for (int i = 0; i < fullGrid.size(); i++)
		gridDists[i] = 1000;
	for (float x = -2.0f; x < 2.1f; x += 0.05f)
	{
		for (float y = -2.0f; y < 2.1f; y += 0.05f)
		{
			for (float z = -2.0f; z < 2.1f; z += 0.05f)
			{
				QuatRep expP = QuatGeoExp(x, y, z);
				std::vector<int> index = GetCoordsFull(expP);
				if (index[0] < 0 || index[1] < 0 || index[2] < 0 || index[3] < 0 || index[0] >= sCount || index[1] >= count || index[2] >= count || index[3] >= count)
					continue;
				expP -= QuatRep((index[0]-(sCount/2))*sGridWidth, (index[1]-(count/2))*gridWidth, (index[2]-(count/2))*gridWidth, (index[3]-(count/2))*gridWidth);
				float dist = (expP.s*expP.s)+(expP.i*expP.i)+(expP.j*expP.j)+(expP.k*expP.k);
				if (dist < gridDists[FullCoordIndex(index)])
				{
					gridDists[FullCoordIndex(index)] = dist;
					fullGrid[FullCoordIndex(index)] = QuatRep(0, x, y, z);
				}
			}
		}
	}

	cGridIndices.resize(count*count*count);
	cGridLength.resize(cGridIndices.size());
	cGridStart.resize(cGridIndices.size());
	cGrid.reserve(cGridIndices.size());
	//condense along s axis
	for (int x = 0; x < count; x++)
	{
		for (int y = 0; y < count; y++)
		{
			for (int z = 0; z < count; z++)
			{
				int l = 0;
				int h = sCount-1;
				while (l < sCount)
				{
					if (gridDists[FullCoordIndex(l, x, y, z)] != 1000)
						break;
					l++;
				}
				while (h >= 0)
				{
					if (gridDists[FullCoordIndex(h, x, y, z)] != 1000)
						break;
					h--;
				}

				cGridIndices[CondensedIndex(x, y, z)] = cGrid.size();
				if ((h-l)+1 <= 0)
				{
					cGridLength[CondensedIndex(x, y, z)] = 0;
					cGridStart[CondensedIndex(x, y, z)] = 0;
				}
				else
				{
					cGridLength[CondensedIndex(x, y, z)] = (h-l)+1;
					cGridStart[CondensedIndex(x, y, z)] = l;
					for (int i = l; i <= h; i++)
						cGrid.push_back(fullGrid[FullCoordIndex(i, x, y, z)]);
				}
			}
		}
	}
	cGrid.shrink_to_fit();

	//assign closest value to empty cGrids
	for (int x = 0; x < count; x++)
	{
		for (int y = 0; y < count; y++)
		{
			for (int z = 0; z < count; z++)
			{
				int length = cGridLength[CondensedIndex(x,y,z)];
				int index = cGridIndices[CondensedIndex(x,y,z)];
				int lastNon = -1;
				
				//Assign initial cGrid to first valid approx
				int i = 0;
				lastNon = 0;
				/* unnecessary because of condensementfor (; i < length; i++)
				{
					if (cGrid[index+i].s != 1)
					{
						lastNon = i;
						for (int j = i-1; j >= 0; j--)
							cGrid[index+j] = cGrid[index+i];
						i++;
						break;
					}
				}*/
				//Assign invalid cGrid to closest valid
				for (; i < length; i++)
				{
					if (cGrid[index+i].s != 1)
					{
						int j = lastNon+1;
						for (; j < (lastNon+i)/2; j++)
							cGrid[index+j] = cGrid[index+lastNon];
						for (; j < i; j++)
							cGrid[index+j] = cGrid[index+i];
						lastNon = i;
					}
				}
				/* unnecessary because of condensement //Assign trailing cGrid to last valid
				if (lastNon != -1)
				{
					for (int j = lastNon+1; j < length; j++)
						cGrid[index+j] = cGrid[index+lastNon];
				}*/

			}
		}
	}
}

void LoopApprox::Generate()
{
	float startDist = 0.2f;

	{
		std::vector<QuatRep> points;
		std::vector<QuatRep> pointsBase;
		std::vector<QuatRep> pointsExp;
		std::vector<QuatRep> pointsBaseExp;
		for (int x = 0; x <= 24; x++) //maybe separate into multiple loops when trying to check approx against array
		{
			QuatRep p(0, -1.0f+(x/12.0f), -1, 0);
			pointsBase.push_back(p*startDist/glm::sqrt(p.SqrMag()));
			p = QuatRep(0, -1.0f+(x/12.0f), 1, 0);
			pointsBase.push_back(p*startDist/glm::sqrt(p.SqrMag()));
		}
		for (int x = 1; x < 24; x++) //maybe separate into multiple loops when trying to check approx against array
		{
			QuatRep p(0, -1, -1.0f+(x/12.0f), 0);
			pointsBase.push_back(p*startDist/glm::sqrt(p.SqrMag()));
			p = QuatRep(0, 1, -1.0f+(x/12.0f), 0);
			pointsBase.push_back(p*startDist/glm::sqrt(p.SqrMag()));
		}
		//                   (furthest dif in 1 axis) * (dim), actual dist is related to sqrt of dim so this works as a maximum
		float ignoreDistBase = (pointsBase[0]-pointsBase[2]).SqrMag()*2; //if sqr mag of points is less than this, then these points are adjacent
		float ignoreDist = 0;

		points.resize(pointsBase.size());
		pointsBaseExp.resize(pointsBase.size());
		pointsExp.resize(pointsBase.size());
		for (int i = 0; i < pointsBase.size(); i++)
		{
			points[i] = QuatRep(0, 0, 0, 0);

			pointsBaseExp[i] = QuatExp(pointsBase[i]);
			pointsExp[i] = QuatRep(1, 0, 0, 0);
		}

		loopAssocMov.resize(points.size());
		for (int i = 0; i < loopAssocMov.size(); i++) loopAssocMov[i] = QuatRep(1, 0, 0, 0);

		for (int count = 0; count < 40; count++)
		{
			ignoreDist += ignoreDistBase;
			for (int i = 0; i < points.size(); i++)
			{
				points[i] += pointsBase[i];
				pointsExp[i] *= pointsBaseExp[i];
			}

			for (int i = 0; i < points.size(); i++)
			{
				if (loopAssocMov[i].s != 1)
					continue;

				QuatRep iQ = points[i];
				QuatRep eiQ = pointsExp[i];

				float maxAdjacent = 0;
				float minNonAdj = INFINITY;
				int minNonAdjIndex = i;
				for (int j = 0; j < points.size(); j++)
				{
					QuatRep jQ = points[j];
					QuatRep ejQ = pointsExp[j];
					if (j == i)
						continue;
					QuatRep dif = iQ-jQ;
					QuatRep eDif = eiQ-ejQ;
					if (dif.SqrMag() <= ignoreDist)
					{
						if (eDif.SqrMag() > maxAdjacent)
							maxAdjacent = eDif.SqrMag();
					}
					else
					{
						if (eDif.SqrMag() < minNonAdj)
						{
							minNonAdj = eDif.SqrMag();
							minNonAdjIndex = j;
						}
					}
				}

				if (minNonAdj < maxAdjacent)
				{
					loopAssocMov[i] = points[minNonAdjIndex];
					loopAssocMovMag[i] = glm::sqrt(points[i].SqrMag());
				}
			}
		}
	}

	{
		std::vector<QuatRep> points;
		std::vector<QuatRep> pointsBase;
		std::vector<QuatRep> pointsExp;
		std::vector<QuatRep> pointsBaseExp;
		{
			QuatRep p(0, 0, 0, -1);
			pointsBase.push_back(p*startDist/glm::sqrt(p.SqrMag()));
			p = QuatRep(0, 0, 0, 1);
			pointsBase.push_back(p*startDist/glm::sqrt(p.SqrMag()));
		}
		//                   (furthest dif in 1 axis) * (dim), actual dist is related to sqrt of dim so this works as a maximum
		float ignoreDistBase = 0; //if sqr mag of points is less than this, then these points are adjacent
		float ignoreDist = 0;

		points.resize(pointsBase.size());
		pointsBaseExp.resize(pointsBase.size());
		pointsExp.resize(pointsBase.size());
		for (int i = 0; i < pointsBase.size(); i++)
		{
			points[i] = QuatRep(0,0,0,0);

			pointsBaseExp[i] = QuatExp(pointsBase[i]);
			pointsExp[i] = QuatRep(1,0,0,0);
		}

		loopAssocRot.resize(points.size());
		for (int i = 0; i < loopAssocRot.size(); i++) loopAssocRot[i] = QuatRep(1, 0, 0, 0);

		for (int count = 0; count < 40; count++)
		{
			ignoreDist += ignoreDistBase;
			for (int i = 0; i < points.size(); i++)
			{
				points[i] += pointsBase[i];
				pointsExp[i] *= pointsBaseExp[i];
			}

			for (int i = 0; i < points.size(); i++)
			{
				if (loopAssocRot[i].s != 1)
					continue;

				QuatRep iQ = points[i];
				QuatRep eiQ = pointsExp[i];

				float maxAdjacent = (pointsBaseExp[0]-pointsBaseExp[1]).SqrMag()*0.9f;
				float minNonAdj = INFINITY;
				int minNonAdjIndex = i;
				for (int j = 0; j < points.size(); j++)
				{
					QuatRep jQ = points[j];
					QuatRep ejQ = pointsExp[j];
					if (j == i)
						continue;
					QuatRep dif = iQ-jQ;
					QuatRep eDif = eiQ-ejQ;
					if (dif.SqrMag() <= ignoreDist)
					{
						if (eDif.SqrMag() > maxAdjacent)
							maxAdjacent = eDif.SqrMag();
					}
					else
					{
						if (eDif.SqrMag() < minNonAdj)
						{
							minNonAdj = eDif.SqrMag();
							minNonAdjIndex = j;
						}
					}
				}

				if (minNonAdj < maxAdjacent)
				{
					loopAssocRot[i] = points[minNonAdjIndex];
					loopAssocRotMag[i] = glm::sqrt(points[i].SqrMag());
				}
			}
		}
	}
}

QuatRep LoopApprox::TryLoop(QuatRep trans)
{
	//convert trans to cubal coords
	//Get loop mag and assoc using grid indexing depending on face of cube
	//if no assoc, ignore
	//if glm::sqrt(trans.SqrMag()) < loop mag * 1.1, then ignore
	//otherwise, trans = assoc*(2-transMag) (assumes aligned tangent space is direct negation, e.g. sphere or axes of torus : exp(1 + 0.1) = exp(1 - 0.1) )
	return QuatRep();
}
