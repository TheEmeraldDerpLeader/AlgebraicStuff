#include <vector>
#include <FlatSymExp.hpp>
#include <SymExp.hpp>

#include <Helpers.hpp>

std::vector<SymExp> GenRule(int p, int i1, int i2, int i3)
{
	std::vector<SymExp> hold; for (int i = 0; i < p+1; i++) hold.push_back(SymExp());

	for (int i = 0; i < p+1; i++) //component of first product i1*i2
	{
		for (int j = 0; j < p+1; j++) //component of first product times i3
		{
			if (i == 0) //first component is the scalar
			{
				if (i3+1 == j) //scalar * i3 only ever has a component of i3
					hold[j].terms.push_back(Product(1,(i1*p*(p+1))+(i2*(p+1))+i));
				continue;
			}
			// a_(i1,i2,i)*a_(i,i3,j)
			//hold[j].terms.push_back(Product(prodVals[(i1*p*(p+1))+(i2*(p+1))+i]*prodVals[(i*p*(p+1))+(i3*(p+1))+j]));
			hold[j].terms.push_back(Product());
			Product& prod = hold[j].terms.back();
			prod.MultId((i1*p*(p+1))+(i2*(p+1))+i);
			prod.MultId(((i-1)*p*(p+1))+(i3*(p+1))+j);
		}
	}

	//repeat for i1*(i2*i3)
	for (int i = 0; i < p+1; i++) //component of first product i2*i3
	{
		for (int j = 0; j < p+1; j++) //component of i1 times first product
		{
			if (i == 0) //first component is the scalar
			{
				if (i1+1 == j) //i1 * scalar only ever has a component of i1
					hold[j].terms.push_back(Product(-1,(i2*p*(p+1))+(i3*(p+1))+i));
				continue;
			}
			// a_(i1,i2,i)*a_(i,i3,j)
			//hold[j].terms.push_back(Product(prodVals[(i1*p*(p+1))+(i2*(p+1))+i]*prodVals[(i*p*(p+1))+(i3*(p+1))+j]));
			hold[j].terms.push_back(Product(-1));
			Product& prod = hold[j].terms.back();
			prod.MultId((i2*p*(p+1))+(i3*(p+1))+i);
			prod.MultId((i1*p*(p+1))+((i-1)*(p+1))+j);
		}
	}

	for (int i = 0; i < hold.size(); i++)
		hold[i].Simplify();

	return hold;
}
void LeTest();
bool VerifyDet(Vector<float>& trans, int dim)
{
	Vector<float> mat = trans;
	VectorRef<float> cur(0,nullptr);
	for (int i = 0; i < dim; i++)
	{
		cur = VectorRef<float>(dim+1,mat.data()+((dim+1)*i));
		cur[0] = 0; //project onto 1
		for (int j = i-1; j >= 0; j--) //project onto previous vecs
		{
			VectorRef<float> col(dim+1, mat.data()+(j*(dim+1)));
			Vector<float> dif = Mul(col, Dot(col, cur)/Dot(col, col));
			SubEq(cur, dif);
		}
		bool isNonZero = false;
		for (int i = 0; i < dim+1; i++)
		{
			if (cur[i] > 0.01f || cur[i] < -0.01f)
			{
				isNonZero = true;
				break;
			}
		}
		if (isNonZero == false)
			return false;
	}
	return true;
}
int GetKernelDim (Vector<float>& mat, int height)
{
	int width = mat.size()/height;
	int kernDim = width;
	Vector<float> matCopy = mat;
	for (int i = 0; i < width; i++) //project vectors onto previous vectors, if nonzero reduce kernel dim
	{
		VectorRef<float> cur(height, matCopy.data()+(i*height));
		for (int j = i-1; j >= 0; j--)
		{
			VectorRef<float> col(height, matCopy.data()+(j*height));
			float div = Dot(col, col);
			if (div == 0)
				continue;
			Vector<float> dif = Mul(col, Dot(col, cur)/div);
			SubEq(cur, dif);
		}
		bool isNonZero = false;
		for (int i = 0; i < height; i++)
		{
			if (cur[i] > 0.01f || cur[i] < -0.01f)
			{
				isNonZero = true;
				break;
			}
		}
		if (isNonZero == false)
			kernDim--;
	}
	return kernDim;
}
int GetKernelDimTransp (Vector<float>& mat, int height)
{
	int width = mat.size()/height;
	Vector<float> matCopy(mat.size()); for (int i = 0; i < width; i++) for (int j = 0; j < height; j++) matCopy[(j*width)+i] = mat[(i*height)+j];
	height = width;
	width = mat.size()/height;
	int kernDim = width;
	for (int i = 0; i < width; i++) //project vectors onto previous vectors, if nonzero reduce kernel dim
	{
		VectorRef<float> cur(height, matCopy.data()+(i*height));
		for (int j = i-1; j >= 0; j--)
		{
			VectorRef<float> col(height, matCopy.data()+(j*height));
			float div = Dot(col, col);
			if (div == 0)
				continue;
			Vector<float> dif = Mul(col, Dot(col, cur)/div);
			SubEq(cur, dif);
		}
		bool isNonZero = false;
		for (int i = 0; i < height; i++)
		{
			if (cur[i] > 0.01f || cur[i] < -0.01f)
			{
				isNonZero = true;
				break;
			}
		}
		if (isNonZero == false)
			kernDim--;
	}
	return kernDim;
}

std::vector<SymExp> GenIsoRule(int p, VectorRef<float> prodRulesBase, VectorRef<float> prodRulesTarget) //gives system of quadratic expressions which describe the general product rules after applying a linear trans
{

	std::vector<SymExp> hold; hold.resize(p*p*(p+1));
	
	// if i'^2 = i' + j', then we need f(i)^2 = f(i) + f(j) => f(i)^2 - f(i) - f(j) = 0

	//complex case
	for (int i = 0; i < p; i++) //first basis vec being mapped
	{
		for (int j = 0; j < p; j++) //second basis vec being mapped
		{
			for (int k = 0; k < p+1; k++) //components of mult
			{
				SymExp& prodRule = hold[(i*p*(p+1))+(j*(p+1))+k];
				for (int il = 0; il < p+1; il++) //component of first basis vec trans
				{
					for (int jl = 0; jl < p+1; jl++) //component of second basis vec trans
					{
						Product prod;
						prod.MultId((i*(p+1))+il); prod.MultId((j*(p+1))+jl);

						//handle scalar cases
						if (il == 0)
						{
							if (jl == k)
								prodRule.terms.push_back(prod);
							continue;
						}
						else if (jl == 0)
						{
							if (il == k)
								prodRule.terms.push_back(prod);
							continue;
						}

						prod.coeff = prodRulesBase[((il-1)*(p)*(p+1))+((jl-1)*(p+1))+k];
						prodRule.terms.push_back(prod);

					}
				}
				//subtract linear terms
				if (k == 0)
					prodRule.scalar -= prodRulesTarget[(i*p*(p+1))+(j*(p+1))+0];
				for (int l = 0; l < p; l++) //index in sum
				{
					//product of trans should equal target product sum, e.g. f(i)f(j) = f(i) + f(j)
					prodRule.terms.push_back(Product(-prodRulesTarget[(i*p*(p+1))+(j*(p+1))+(l+1)], (l*(p+1))+k));
				}
				prodRule.Simplify(); //get rid of 0 terms
			}
		}
	}
	
	return hold;
}

SymExp ApplyDivSub(SymExp& exp, int idToReplace, SymExp& num, SymExp& denom)
{
	SymExp out;

	//get maximum power of id in all terms
	int idMaxPow = 0;
	for (int i = 0; i < exp.terms.size(); i++)
		for (int j = 0; j < exp.terms[i].ids.size(); j++)
			if (exp.terms[i].ids[j] == idToReplace && exp.terms[i].pows[j] > idMaxPow)
				idMaxPow = exp.terms[i].pows[j];

	//substitute id for num, multiply by denom for missing y (e.g. x + y - 1 => x + 1/y' - 1 => xy' + 1 - y'
	SymExp hold(exp.scalar);
	for (int j = 0; j < idMaxPow; j++)
		hold *= denom;
	out += hold;

	for (int i = 0; i < exp.terms.size(); i++)
	{
		int pow = 0;
		hold = exp.terms[i];
		for (int j = 0; j < exp.terms[i].ids.size(); j++)
			if (exp.terms[i].ids[j] == idToReplace)
			{
				pow = exp.terms[i].pows[j];
				hold.terms[0].ids.erase(hold.terms[0].ids.begin()+j);
				hold.terms[0].pows.erase(hold.terms[0].pows.begin()+j);
			}
		
		int j = 0;
		for (; j < pow; j++)
			hold *= num;
		for (; j < idMaxPow; j++)
			hold *= denom;
		
		out += hold;

	}

	return out;
}

// y = num/denom, where the old id of y is reused for y', with y' = 1/exp
void SolveForVar(SymExp& exp, int idToSolve, SymExp& numOut, SymExp& denomOut)
{
	numOut = SymExp(1-exp.scalar);
	denomOut = SymExp();

	for (int i = 0; i < exp.terms.size(); i++)
	{
		int powOfId = 0;
		Product hold;
		for (int j = 0; j < exp.terms[i].ids.size(); j++)
			if (exp.terms[i].ids[j] == idToSolve)
			{
				powOfId = exp.terms[i].pows[j];
				hold = exp.terms[i];
				hold.ids.erase(hold.ids.begin()+j);
				hold.pows.erase(hold.pows.begin()+j);
			}

		if (powOfId == 0)
			numOut -= exp.terms[i].MultId(idToSolve);
		if (powOfId == 1)
			denomOut += hold;
		if (powOfId != 0 && powOfId != 1)
		{
			numOut = SymExp();
			denomOut = SymExp();
			return;
		}
		
	}
	denomOut *= Product(1, 1, 1);
}

int CPPBindingTest()
{
	int p = 2;
	std::vector<SymExp> rules = GenRule(p, 0, 0, 0);
	SymParser multParser;
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < p+1; k++)
			{
				multParser.Add((i*p*(p+1))+(j*(p+1))+k, std::string("a_(")+std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)+")");
			}
		}
	}

	std::cout << "Associative rules for iii:\n" << rules[0].ToString(multParser) << '\n' << rules[1].ToString(multParser) << '\n';

	std::cout << "\n\n" << "All associative rules for p = " << p << ":\n";
	rules.clear();
	for (int i1 = 0; i1 < p; i1++) 
	{
		for (int i2 = 0; i2 < p; i2++) 
		{
			for (int i3 = 0; i3 < p; i3++)
			{
				std::vector<SymExp> rulesTemp = GenRule(p, i1, i2, i3);
				for (int i = 0; i < p+1; i++)
				{
					rules.push_back(rulesTemp[i]);
					std::cout << rules.back().ToString(multParser) << '\n';
				}
			}
		}
	}

	std::cout << '\n';

	Vector<float> ruleTest(p*p*(p+1));
	Vector<int> ruleIds(p*p*(p+1));
	Vector<float> out(rules.size());
	bool foundSol = false;
	int count = 0;

	while (foundSol == false && count < 800) //usually finds sol first try
	{
		for (int i = 0; i < ruleTest.size(); i++)
		{
			ruleTest[i] = randF(-6, 6);
			//ruleTest[i] = 0;
			ruleIds[i] = i;
		}
		//ruleTest[1] = 1;
		//ruleTest[5] = 1;
		//ruleTest[7] = -2.1;
		//ruleTest[11] = -2;

		ruleTest = NMnTom(rules, ruleIds, ruleTest);
		out = SclEvalVec(rules, ruleIds, ruleTest);
		
		foundSol = true;
		for (int i = 0; i < out.size(); i++)
		{
			if (out[i] > 0.01 || out[i] < -0.01)
			{
				foundSol = false;
				break;
			}
		}
		count++;
	}
	//create equations to test if two polys are iso

	Vector<float> baseRules(2); baseRules[0] = 0; baseRules[1] = 0;
	Vector<float> targetRules(2); targetRules[0] = -1; targetRules[1] = 2;
	Vector<float> trans(2); trans[0] = 0; trans[1] = 1;
	Vector<int> transIds(2); transIds[0] = 0; transIds[1] = 1;
	std::vector<SymExp> isoRules = GenIsoRule(1, baseRules, targetRules);
	std::vector<SymExp> isoRulesR = GenIsoRule(1, targetRules, baseRules);
	for (int i = 0; i < isoRules.size(); i++)
		std::cout << isoRules[i].ToString() << '\n';

	//To Do: get a random rule and check if it is iso to any previous rules. If not, add it to the list. If so, keep the one with the abs mag closest to 1.

	//works
	bool isIso = true;
	trans = NMnTom(isoRules, transIds, trans);
	for (int i = 0; i < isoRules.size(); i++)
	{
		if (std::abs(isoRules[i].SclEval(transIds, trans)) > 0.0001f)
		{
			isIso = false;
			break;
		}
	}

	trans[0] = 0; trans[1] = 1;
	trans = NMnTom(isoRulesR, transIds, trans);
	for (int i = 0; i < isoRulesR.size(); i++)
	{
		if (std::abs(isoRulesR[i].SclEval(transIds, trans)) > 0.0001f)
		{
			isIso = false;
			break;
		}
	}
	if (isIso == true)
		std::cout << "Rules are isomorphic\n";
	else
		std::cout << "Rules aren't isomorphic\n";
	LeTest();
	return 0;
}

Vector<float> Basis1FromArbUnitial(Vector<float>& arbUnitialRule, int dim)
{
	if (dim*dim*dim != arbUnitialRule.size())
	{
		std::cout << "Incorrect size for arbUnitialRule in Basis1FromArbUnitial()\n";
		return Vector<float>();
	}

	std::vector<SymExp> identSys; identSys.resize(dim*dim);
	for (int i = 0; i < dim; i++) //initial component for i*(...)
	{
		for (int j = 0; j < dim; j++) //evaluated component, e.g. p<...>_j
		{
			if (i == j)
				identSys[(i*dim)+j].scalar = -1;
			for (int k = 0; k < dim; k++) //component of basisI to consider
			{
				identSys[(i*dim)+j].terms.push_back(Product(arbUnitialRule[(i*dim*dim)+(k*dim)+j],k));
			}
		}
	}
	// i*(ai + bj + ...) = i
	Vector<float> basis1(dim);
	Vector<int> ids; for (int i = 0; i < dim; i++) ids.push_back(i);
	basis1 = NMnTom(identSys, ids, basis1);

	Vector2D<float> basisO((dim-1),dim);
	for (int i = 0; i < dim-1; i++)
	{
		VectorRef<float> b = basisO.GetCol(i);
		b[i+1] = 1;
		for (int j = i-1; j >= -1; j--)
		{
			VectorRef<float> ba(0, nullptr);
			if (j == -1)
				ba = basis1;
			else
				ba = basisO.GetCol(j);
			float coeff = Dot(b, ba)/Dot(ba, ba);
			for (int i = 0; i < dim; i++)
				b[i] -= ba[i]*coeff;
		}
	}

	Vector<float> rule((dim-1)*(dim-1)*dim);
	//i*j
	for (int i = 0; i < dim-1; i++) 
	{
		for (int j = 0; j < dim-1; j++)
		{
			Vector<float> prod(dim);
			VectorRef<float> v1 = basisO.GetCol(i);
			VectorRef<float> v2 = basisO.GetCol(j);
			for (int k = 0; k < dim; k++)
				for (int i = 0; i < dim; i++)
					for (int j = 0; j < dim; j++)
						prod[k] += arbUnitialRule[(i*dim*dim)+(j*dim)+k]*v1[i]*v2[j];

			for (int i2 = -1; i2 < dim-1; i2++)
			{
				VectorRef<float> v(0, nullptr);
				if (i2 == -1)
					v = basis1;
				else
					v = basisO.GetCol(i2);
				float coeff = Dot(prod,v)/Dot(v,v);
				for (int j = 0; j < dim; j++)
					prod[j] -= coeff*v[j];
				rule[(i*dim*(dim-1))+(j*dim)+(i2+1)] = coeff;
			}
		}
	}
	return rule;
}

//true if found valid and nontrivial trans
bool SearchIso(int dim, Vector<float>& trans, Vector<float>& transTest, FlatSymExp& isoRules, FlatSymExp& isoRulesGrad)
{
	thread_local Vector<int> transIds;
	transIds.resize(trans.size()); for (int i = 0; i < transIds.size(); i++) transIds[i] = i;
	trans = isoRules.NewtonsMethodSolve(isoRulesGrad, trans); //look for a possible trans between alges
	transTest = isoRules.SclEval(trans); //check if trans is valid
	bool wasIso = true;
	for (int i = 0; i < transTest.size(); i++)
	{
		if (transTest[i] < -0.00001f || transTest[i] > 0.00001f)
		{
			wasIso = false;
			break;
		}
	}
	//if (wasIso == false) //algo should converge if possible, but didn't, so these must be different alges
	//	break; //move on to next space
	if (wasIso == true && VerifyDet(trans, dim) == true) //algo converged and valid test, spaces must be iso
		return true;
	else
		return false;
}

std::vector<float> ConvergenceAtPoint(Vector2D<SymExp> grad, std::vector<int> ids, std::vector<float> p)
{

	int idC = grad.width;
	int expC = grad.height;
	Vector2D<float> vBasis(idC, idC);
	Vector2D<float> wBasis(expC, expC);
	Vector<float> wDivForV(expC);

	Vector2D<float> gradEval = SclEvalVec2D(grad, ids, p);

	//transform gradEval into new basis with ortho wi
	int vCount = 0;
	for (int i = 0; i < expC && vCount < idC; i++) //i is the grad row used, vCount is the number of vis currently found
	{
		//calculate starting vi and wi
		VectorRef<float> initialV = vBasis.GetCol(vCount);
		for (int j = 0; j < idC; j++)
			initialV[j] = gradEval[(j*expC)+i];
		VectorRef<float> initialW = wBasis.GetCol(vCount);
		for (int i = 0; i < expC; i++)
			initialW[i] = 0;

		//calculate initialW from vs and ws
		for (int i = 0; i < idC; i++)
		{
			VectorRef<float> wTemp = gradEval.GetCol(i);
			for (int j = 0; j < expC; j++)
				initialW[j] += wTemp[j]*initialV[i];
		}

		//antiproject wi and vi
		for (int i = vCount-1; i >= 0; i--)
		{
			VectorRef<float> antiV = vBasis.GetCol(i);
			VectorRef<float> antiW = wBasis.GetCol(i);
			float coeff = 0;
			for (int i = 0; i < expC; i++)
				coeff += initialW[i]*antiW[i];
			coeff /= wDivForV[i];
			//project wi
			for (int i = 0; i < expC; i++)
				initialW[i] -= antiW[i]*coeff;
			//project vi
			for (int i = 0; i < idC; i++)
				initialV[i] -= antiV[i]*coeff;
		}

		wDivForV[vCount] = 0;
		for (int i = 0; i < expC; i++)
			wDivForV[vCount] += initialW[i]*initialW[i];

		//if wDivForV = 0, then vi and wi would be 0, so this vector should be skipped 
		if (wDivForV[vCount] > 0.0000001 || wDivForV[vCount] < -0.0000001) //ferr is relevant
			vCount++;
	}

	std::vector<float> out; out.resize(vCount);

	//magnitude of vBasis direction determines convergence in that direction
	for (int i = 0; i < vCount; i++)
	{
		float mag = 0;
		for (int j = 0; j < idC; j++)
			mag += vBasis[j+(i*idC)]*vBasis[j+(i*idC)];

		out[i] = mag;
	}

	return out;
}

void LeTest()
{
	Vector<float> solver(27); 
	//solver[0] = 1; solver[13] = 1; solver[26] = 1;
	solver[0] = 1; solver[13] = 1; solver[17] = 1; solver[23] = 1; solver[25] = -1;
	Vector<float> newRules = Basis1FromArbUnitial(solver, 3);


	int dim = 2;
	Vector<float> rules;
	Vector<float> rule(dim*dim*(dim+1));
	Vector<int> ruleIds(rule.size()); for (int i = 0; i < rule.size(); i++) ruleIds[i] = i;
	Vector<int> ruleFails;
	Vector<int> ruleIsos;

	//setup associativity search
	Vector<SymExp> assocEq; assocEq.reserve(dim*dim*dim*(dim+1));
	for (int i1 = 0; i1 < dim; i1++) 
	{
		for (int i2 = 0; i2 < dim; i2++) 
		{
			for (int i3 = 0; i3 < dim; i3++)
			{
				std::vector<SymExp> rulesTemp = GenRule(dim, i1, i2, i3);
				for (int i = 0; i < dim+1; i++)
				{
					assocEq.push_back(rulesTemp[i]);
				}
			}
		}
	}
	Vector<float> assocOut(assocEq.size());

	Vector<float> trans(dim*(dim+1));
	Vector<int> transIds(trans.size()); for (int i = 0; i < transIds.size(); i++) transIds[i] = i;
	Vector<float> transTest(dim*(dim+1));

	//for (int i = 0; i < 12; i++)
	//	rules.push_back(newRules[i]);
	//rules.push_back(1); rules.push_back(0); //hypercomplex
	//rules.push_back(-1); rules.push_back(0);
	//rules.push_back(-1); rules.push_back(2);

	//for (int i = 0; i < 12; i++) rules.push_back(0); rules[2] = 1; rules[3] = 1; rules[6] = 1; rules[10] = 1;
	//for (int i = 0; i < 12; i++) rules.push_back(0); rules[2] = 0; rules[4] = 1; rules[7] = -1; rules[9] = 1;
	
	SymExpTable unitRestrict; //only consider i^2 = i + ..., j^2 = j + ...
	//every rule is iso to rules in this form, a path within rules corresponds to another path in rules of this form as long as the path never contains i^2 = 0*i + ...
	//for (int i = 0; i < dim; i++)
	//	unitRestrict.Add((i*dim*(dim+1))+(i*(dim+1))+(i+1), 1);
	 
	//for each (ij) = -(ji), there is redundancy in the transformation subsurface
	//i^2 = - i^2 = 0 (reduction of 1)
	//for (int i = 0; i < dim+1; i++)
	//	unitRestrict.Add(i+0, 0);
	//ij = -ji (reduction of 2)
	for (int i = 0; i < dim+1; i++)
		unitRestrict.Add(i+3, SymExp(Product(-1,i+6)));
	//j^2 = 0
	//for (int i = 0; i < dim+1; i++)
	//	unitRestrict.Add(i+9, 0);
	//i reducible
	//unitRestrict.Add(10, 0);

	//adding eq constraints rather than substituting variables
	//i^2 = - i^2 = 0 (reduction of 1)
	//for (int i = 0; i < dim+1; i++)
	//	assocEq.push_back(SymExp(Product(1,i+0)));
	//ij = -ji (reduction of 2)
	//for (int i = 0; i < dim+1; i++)
	//{
	//	assocEq.push_back(SymExp(Product(1,i+3)));
	//	assocEq.back().terms.push_back(Product(1,i+6));
	//}
	//j^2 = 0
	//for (int i = 0; i < dim+1; i++)
	//	assocEq.push_back(SymExp(Product(1,i+9)));

	for (int i = 0; i < assocEq.size(); i++)
		assocEq[i] = assocEq[i].Eval(unitRestrict);

	for (int count = 0; count < 16000; count++)
	{
		//reset rule fail counts
		for (int i = 0; i < ruleFails.size(); i++) ruleFails[i] = 0;

		//generate a random rule
		for (int i = 0; i < rule.size(); i++)
			rule[i] = randF(4.0f, 8.0f)*((randI(0,1)*2)-1);

		if (count == 0) //i^2 = i, j^2 = j, k^2 = k
		{
			for (int i = 0; i < 12; i++) rule[i] = 0;
			rule[0] = 2; rule[1] = 3; rule[5] = -3; rule[8] = -3; rule[9] = 1.5f; rule[10] = -2.25;
			for (int i = 0; i < 12; i++) rule[i] /= 9.0f;
			
		}
		if (count == 1) //https://arxiv.org/pdf/1903.01623 unitial algebra from iso over R part in summary
		{
			for (int i = 0; i < 12; i++) rule[i] = 0;
			rule[0] = 1; rule[5] = 2; rule[8] = 2; rule[9] = -2; rule[10] = -4;
			for (int i = 0; i < 12; i++) rule[i] /= 4.0f;
		}
		if (count == 2)
		{
			Vector<float> solver(27); 
			solver[0] = 1; solver[4] = 1; solver[8] = 1; solver[10] = 1; solver[20] = 1;
			rule = Basis1FromArbUnitial(solver, 3);
			//continue;
		}
		if (count == 3)
		{
			Vector<float> solver(27); 
			solver[0] = 1; solver[4] = 1; solver[8] = 1; solver[10] = 1; solver[16] = 1;
			solver[20] = 1; solver[22] = -1; solver[24] = 1;
			rule = Basis1FromArbUnitial(solver, 3);
		}
		if (count == 4)
		{
			Vector<float> solver(27); 
			solver[0] = 1; solver[13] = 1; solver[17] = 1; solver[23] = 1;
			rule = Basis1FromArbUnitial(solver, 3);
			//continue;
		}
		if (count == 5)
		{
			Vector<float> solver(27); 
			solver[0] = 1; solver[4] = 1; solver[8] = 1; solver[10] = 1;
			solver[14] = 1; solver[20] = 1;
			rule = Basis1FromArbUnitial(solver, 3);
			//continue;
		}
		if (count == 6)
		{
			for (int i = 0; i < 12; i++) rule[i] = 0;
			//rule[0] = -1; rule[1] = 2; rule[3] = -1; rule[4] = 1; rule[5] = 1; rule[6] = 1; rule[7] = -1; rule[8] = 1; rule[9] = 1;
			
			rule[10] = 1; rule[11] = 0.01f; //this hits both 4 and 5
		}

		//rule[0] = 0; rule[1] = 0; //zero div algs may give NaN outputs in NM, but can still determine if spaces are different or same maybe?
		//for (int i = 0; i < 12; i++) rule[i] = rules[i]; //rule[0] = -1; rule[1] = 2; rule[2] = 1; rule[3] = 1; rule[5] = 1; rule[6] = 1; rule[8] = 1; rule[9] = -1; rule[10] = 1;
		//for (int i = 0; i < 12; i++) rule[i] = 0; rule[2] = 1; rule[3] = -1; rule[6] = -1; rule[10] = -1;
		//for (int i = 0; i < 12; i++) rule[i] = 0.3f; rule[2] = 2.5f; rule[4] = 1; rule[7] = -1; rule[9] = 1;

		if (count <= 6)
			continue;

		if (count > 5)
		{
			for (int i = 0; i < unitRestrict.lookup.size(); i++)
				if (unitRestrict.exps[i].terms.size() == 0)
					rule[unitRestrict.lookup[i]] = unitRestrict.exps[i].scalar;
				else
					rule[unitRestrict.lookup[i]] = -rule[unitRestrict.exps[i].terms[0].ids[0]];

			//rule = NMnTom(assocEq, ruleIds, rule);
			rule = NMnTomFailSlowConvergence(assocEq, ruleIds, rule, 0, 0.1f, 4, 0.000001);
			for (int i = 0; i < unitRestrict.lookup.size(); i++) //NM won't update substituted values
				if (unitRestrict.exps[i].terms.size() != 0)
					rule[unitRestrict.lookup[i]] = -rule[unitRestrict.exps[i].terms[0].ids[0]];
		}
		assocOut = SclEvalVec(assocEq, ruleIds, rule);
		bool conv = true;
		for (int i = 0; i < assocOut.size(); i++)
		{
			if (count <= 5)
				break;
			if (assocOut[i] < -0.00001f || assocOut[i] > 0.00001f)
			{
				conv = false;
				break;
			}
		}

		//ignore rules that are close to slow conv areas
		std::vector<float> gradDebug = ConvergenceAtPoint(Gradient(assocEq), ruleIds, rule);
		for (int i = 0; i < gradDebug.size(); i++)
			for (int j = i-1; j >= 0; j--)
				if (gradDebug[j+1] > gradDebug[j])
				{
					float hold = gradDebug[j+1]; gradDebug[j+1] = gradDebug[j]; gradDebug[j] = hold;
				}
				else break;
		//if (count <= 6)
		//	conv = false;
		if (gradDebug.back()/gradDebug[0] < 0.1f && count > 5)
			conv = false;

		if (conv == false)
			continue;

		thread_local Vector<FlatSymExp> isoRulesVec;
		thread_local Vector<FlatSymExp> isoRulesVecGrad;
		thread_local Vector<FlatSymExp> isoRulesRVec;
		thread_local Vector<FlatSymExp> isoRulesRVecGrad;
		isoRulesVec.resize(rules.size()/rule.size());
		isoRulesVecGrad.resize(rules.size()/rule.size());
		isoRulesRVec.resize(rules.size()/rule.size());
		isoRulesRVecGrad.resize(rules.size()/rule.size());
		for (int i = 0; i < rules.size()/rule.size(); i++)
		{
			VectorRef<float> targetRule(rule.size(), rules.data()+(rule.size()*i));
			std::vector<SymExp> isoRule = GenIsoRule(dim, rule, targetRule);
			isoRulesVec[i] = FlatSymExp(isoRule);
			isoRulesVecGrad[i] = isoRulesVec[i].Gradient();
			isoRule = GenIsoRule(dim, targetRule, rule);
			isoRulesRVec[i] = FlatSymExp(isoRule);
			isoRulesRVecGrad[i] = isoRulesRVec[i].Gradient();
		}


		//meant to prevent trivial transformations
		//worked, but also broke convergence for valid transformations, so this was completely pointless
		//Maybe could combine isoRule and inverse, just use NMmTonFromFunction to compute inverse coords and inverse gradient
		//for (int i = 0; i < isoRules.size(); i++)
		//	isoRules[i] = ApplyDivSub(isoRules[i],1,num,denom);
		//for (int i = 0; i < isoRules.size(); i++)
		//	isoRulesR[i] = ApplyDivSub(isoRulesR[i],1,num,denom);

		//reset trans
		std::vector<float> initTrans; initTrans.resize(trans.size());
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim+1; j++)
			{
				initTrans[(i*(dim+1))+j] = 0;
				if (i == j-1)
					initTrans[(i*(dim+1))+j] = 1;
			}
		}
		
		bool wasIso = false;
		for (int f = 0; f < 200000; f++) //200000
		{
			if (count <= 5)
				break;
			
			if (f != 0) //randomized trans
			{
				for (int i = 0; i < dim; i++)
					for (int j = 0; j < dim+1; j++)
						initTrans[(i*(dim+1))+j] = randF(-4, 4);
			}

			int i = 0;
			for (;i < rules.size()/rule.size(); i++)
			{
				if (i >= 10)
					break;
				trans = initTrans;
				if (SearchIso(dim, trans, transTest, isoRulesVec[i], isoRulesVecGrad[i]))
				{
					wasIso = true;
					break;
				}
				trans = initTrans;
				if (SearchIso(dim, trans, transTest, isoRulesRVec[i], isoRulesRVecGrad[i]))
				{
					wasIso = true;
					break;
				}
				//both tests were invalid, probably not iso
				ruleFails[i]++;
			}
			if (wasIso == false)
				continue;

			thread_local int failmax = 0;
			if (ruleFails[0] > failmax)
				failmax = ruleFails[0];
			
			//algo converged more than it didn't, spaces are iso
			wasIso = true;
			ruleIsos[i]++;
			break;
			
		}
		if (wasIso == false)
		{
			for (int i = 0; i < rule.size(); i++)
				rules.push_back(rule[i]);
			ruleFails.push_back(0);
			ruleIsos.push_back(0);
		}
	}

 	std::cout << rules.size()/rule.size() << '\n';
}

/*solutions for p = 2

//paper has i^2 = i, ij = j, j^2 = -i, ji = j
*/