#include "QuatRep.hpp"

#include <SymExp.hpp>
#include <FlatSymExp.hpp>
#include <ArbAlgRep.hpp>

QuatRep QuatRep::operator*(QuatRep quat)
{
	QuatRep hold;
	//hold.s = (s*quat.s)+(i*quat.i)+(j*quat.j)-(k*quat.k);
	//hold.i = (s*quat.i)+(i*quat.s)+(j*quat.k)-(k*quat.j);
	//hold.j = (s*quat.j)-(i*quat.k)+(j*quat.s)+(k*quat.i);
	//hold.k = (s*quat.k)-(i*quat.j)+(j*quat.i)+(k*quat.s);
	//hold.s = (s*quat.s)-(i*quat.i)-(j*quat.j)-(k*quat.k);
	//hold.i = (s*quat.i)+(i*quat.s)+(j*quat.k)-(k*quat.j);
	//hold.j = (s*quat.j)-(i*quat.k)+(j*quat.s)+(k*quat.i);
	//hold.k = (s*quat.k)+(i*quat.j)-(j*quat.i)+(k*quat.s);
	//hold.s = (s*quat.s)-(k*quat.k);
	//hold.i = (s*quat.i)+(i*quat.s)-(j*quat.k)+(k*quat.j);
	//hold.j = (s*quat.j)+(i*quat.k)+(j*quat.s)-(k*quat.i);
	//hold.k = (s*quat.k)+(k*quat.s);
	//hold.s = (s*quat.s)+(i*quat.j)+(j*quat.i);
	//hold.i = (s*quat.i)+(i*quat.s);
	//hold.j = (s*quat.j)+(j*quat.s);
	//hold.k = (s*quat.k)+(k*quat.s);

	//rule[0] = 2; rule[1] = 3; rule[5] = -3; rule[8] = -3; rule[9] = 1.5f; rule[10] = -2.25;
	//hold.s = (s*quat.s)+(2.0f*i*quat.i)+(1.5f*j*quat.j);
	//hold.i = (s*quat.i)+(i*quat.s)+(3.0f*i*quat.i)-(3.0f*i*quat.j)-(3.0f*j*quat.i)-(2.25f*j*quat.j);
	//hold.j = (s*quat.j)+(j*quat.s);
	//hold.k = (s*quat.k)+(k*quat.s);

	//rule[0] = 1; rule[5] = 2; rule[8] = 2; rule[9] = -2; rule[10] = -4;
	//hold.s = (s*quat.s)+(i*quat.i)-(4.0f*j*quat.j);
	//hold.i = (s*quat.i)+(i*quat.s);
	//hold.j = (s*quat.j)+(j*quat.s)+(2.0f*i*quat.j)+(2.0f*j*quat.i);
	//hold.k = (s*quat.k)+(k*quat.s);

	//Cool space
	//hold.s = (s*quat.s)+(j*quat.j);
	//hold.i = (s*quat.i)+(i*quat.s)+(i*quat.j)-(j*quat.i);
	//hold.j = (s*quat.j)+(j*quat.s);
	//hold.k = (s*quat.k)+(k*quat.s);
	 
	//complex numbers but i' = 1 + i
	//hold.s = (s*quat.s)-(2*i*quat.i)-(j*quat.k)+(k*quat.j)-(j*quat.j)-(k*quat.k);
	//hold.i = (s*quat.i)+(i*quat.s)+(2*i*quat.i)+(j*quat.k)-(k*quat.j);
	//hold.j = (s*quat.j)+(j*quat.s)+(i*quat.j)+(j*quat.i)-(i*quat.k)+(k*quat.i);
	//hold.k = (s*quat.k)+(k*quat.s)+(i*quat.j)-(j*quat.i)+(i*quat.k)+(k*quat.i);

	//hold.s = (s*quat.s);
	//hold.i = (s*quat.i)+(i*quat.s);
	//hold.j = (s*quat.j)+(j*quat.s)+(i*quat.i);
	//hold.k = (s*quat.k)+(k*quat.s);

	//not one of the algs from the paper, this is i2 = j, j2 = i, ij = ji = 1
	hold.s = (s*quat.s)+(i*quat.j)+(j*quat.i);
	hold.i = (s*quat.i)+(i*quat.s)+(j*quat.j);
	hold.j = (s*quat.j)+(j*quat.s)+(i*quat.i);
	hold.k = (s*quat.k)+(k*quat.s);

	//hold.s = (s*quat.s)-(i*quat.i);
	//hold.i = (s*quat.i)+(i*quat.s)+(2*i*quat.i);
	//hold.j = (s*quat.j)+(j*quat.s);
	//hold.k = (s*quat.k)+(k*quat.s);


	//circle times line
	//hold.s = (s*quat.s)-(i*quat.i);
	//hold.i = (s*quat.i)+(i*quat.s);
	//hold.j = (s*quat.j)+(j*quat.s)+(i*quat.j)+(j*quat.i);
	//hold.k = (s*quat.k)+(k*quat.s);

	return hold;
}

QuatRep QuatRep::operator/(QuatRep quat)
{
	return (*this)*quat.Inverse();
}

QuatRep QuatRep::operator+(QuatRep quat)
{
	QuatRep hold;
	hold.s = s+quat.s;
	hold.i = i+quat.i;
	hold.j = j+quat.j;
	hold.k = k+quat.k;
	return hold;
}

QuatRep QuatRep::operator-(QuatRep quat)
{
	QuatRep hold;
	hold.s = s-quat.s;
	hold.i = i-quat.i;
	hold.j = j-quat.j;
	hold.k = k-quat.k;
	return hold;
}

QuatRep& QuatRep::operator*=(QuatRep quat)
{
	*this = (*this)*quat;
	return *this;
}

QuatRep& QuatRep::operator/=(QuatRep quat)
{
	*this = (*this)/quat;
	return *this;
}

QuatRep& QuatRep::operator+=(QuatRep quat)
{
	*this = (*this)+quat;
	return *this;
}

QuatRep& QuatRep::operator-=(QuatRep quat)
{
	*this = (*this)-quat;
	return *this;
}

QuatRep QuatRep::operator*(float scale)
{
	return QuatRep(s*scale, i*scale, j*scale, k*scale);
}

QuatRep QuatRep::operator/(float scale)
{
	return QuatRep(s/scale, i/scale, j/scale, k/scale);
}

QuatRep& QuatRep::operator*=(float scale)
{
	*this = (*this)*scale;
	return *this;
}

QuatRep& QuatRep::operator/=(float scale)
{
	*this = (*this)/scale;
	return *this;
}

thread_local std::vector<float> inverseMat; //x coordinate is sub var multiplier, y coordinate is var, z coordinate is eq
thread_local FlatSymExp sqrtSys;
thread_local std::vector<SymExp> invSys;
//thread_local FlatSymExp invExpSys;
thread_local float invExpDist; 
//thread_local FlatSymExp invGeoExpSys;
thread_local FlatSymExp invGeoExpRRSys;

thread_local FlatSymExp sqrtSysGrad;
//thread_local FlatSymExp invExpSysGrad;
thread_local FlatSymExp invGeoExpSysGrad;
thread_local FlatSymExp invGeoExpRRSysGrad;

thread_local FlatSymExp testInvGeoExpSys;
thread_local FlatSymExp testInvGeoExpSysGrad;

thread_local QuatRep d1;
thread_local QuatRep d2; 
thread_local QuatRep r1; //have r1 = exp(del*d1)*exp(del*d2) - exp(del*(d1+d2)) | lim => 0 == d1*d2 - 0.5*(d1d2 + d2d1) = 0.5(d1d2 - d2d1) 
//static QuatRep r2; exp(1) can never result from d1*d2, as it would imply d1 and d2 are comm (d1 = s*d2^-1), so we never have more than 1 rotation dir

thread_local QuatRep QuatRepDirs::d1;
thread_local QuatRep QuatRepDirs::d2;
thread_local QuatRep QuatRepDirs::r1;

struct __SetupQuat_
{
public: 
	__SetupQuat_()
	{
		thread_local std::vector<SymExp> sqrtSysPre;
		//thread_local std::vector<SymExp> invExpSysPre;
		//thread_local std::vector<SymExp> invGeoExpSysPre;
		thread_local std::vector<SymExp> invGeoExpRRSysPre;
		thread_local std::vector<SymExp> testInvGeoExpSysPre;

		MultTable mt(3);
		//mt.At(0, 0, -1) = 1; mt.At(0, 1, 2) = -1; mt.At(0, 2, 1) = -1; 
		//mt.At(1, 1, -1) = 1; mt.At(1, 0, 2) = 1; mt.At(1, 2, 0) = 1;
		//mt.At(2, 2, -1) = -1; mt.At(2, 0, 1) = 1; mt.At(2, 1, 0) = -1;

		//mt.At(0, 0, -1) = 0; mt.At(0, 1, 2) = 0; mt.At(0, 2, 1) = -1; 
		//mt.At(1, 1, -1) = 0; mt.At(1, 0, 2) = 0; mt.At(1, 2, 0) = -1;
		//mt.At(2, 2, -1) = 1; mt.At(2, 0, 1) = 1; mt.At(2, 1, 0) = 1;

		//mt.At(1, 1, -1) = 1; mt.At(0, 1, 1) = 1; mt.At(1, 0, 1) = -1; 

		QuatRep qBasis[4]; qBasis[0] = QuatRep(1,0,0,0); qBasis[1] = QuatRep(0,1,0,0); qBasis[2] = QuatRep(0,0,1,0); qBasis[3] = QuatRep(0,0,0,1);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				QuatRep prod = qBasis[i+1]*qBasis[j+1];
				mt.At(i, j, -1) = prod.s;
				mt.At(i, j, 0) = prod.i;
				mt.At(i, j, 1) = prod.j;
				mt.At(i, j, 2) = prod.k;
			}
		}

		//get new basis where basis vectors square to unitial subspace
		std::vector<SymExp> newBasisEq; newBasisEq.resize(2);
		newBasisEq[0].scalar = (QuatRep(0, 1, 0, 0)*QuatRep(0, 1, 0, 0)).i;
		newBasisEq[1].scalar = (QuatRep(0, 1, 0, 0)*QuatRep(0, 1, 0, 0)).j;
		newBasisEq[0].terms.push_back(Product(2,0));
		newBasisEq[1].terms.push_back(Product(2, 1)); newBasisEq[1].terms.back().MultId(0);
		newBasisEq[0].terms.push_back(Product((QuatRep(0,0,1,0)*QuatRep(0,0,1,0)).i,1,2));
		newBasisEq[1].terms.push_back(Product((QuatRep(0,0,1,0)*QuatRep(0,0,1,0)).j,1,2));
		newBasisEq[0].terms.push_back(Product((QuatRep(0,1,0,0)*QuatRep(0,0,1,0)).i,1,1));
		newBasisEq[1].terms.push_back(Product((QuatRep(0,1,0,0)*QuatRep(0,0,1,0)).j,1,1));
		newBasisEq[0].terms.push_back(Product((QuatRep(0,0,1,0)*QuatRep(0,1,0,0)).i,1,1));
		newBasisEq[1].terms.push_back(Product((QuatRep(0,0,1,0)*QuatRep(0,1,0,0)).j,1,1));

		std::vector<float> initial; initial.push_back(0.35f); initial.push_back(0.47f);
		FlatSymExp test(newBasisEq);
		FlatSymExp testGrad = test.Gradient();
		initial = test.NewtonsMethodSolve(testGrad, initial);
 
		std::vector<SymExp> varExp;
		std::vector<SymExp> fixedExp;
		varExp.push_back(SymExp(Product(1, 0, 1)));
		varExp.push_back(SymExp(Product(1, 1, 1)));
		varExp.push_back(SymExp(Product(1, 2, 1)));
		varExp.push_back(SymExp(Product(1, 3, 1)));
		fixedExp.push_back(SymExp(Product(1, 4, 1)));
		fixedExp.push_back(SymExp(Product(1, 5, 1)));
		fixedExp.push_back(SymExp(Product(1, 6, 1)));
		fixedExp.push_back(SymExp(Product(1, 7, 1)));

			//look at ScalarMinimizing.png to find ideal c for exp(i + c)
			//trying to minimize squared mag of exp(1000i + 1000c) maybe
			//There is always some c such that squared mag of exp(xi + xc) derivative is 0 at x = 0, this is ideal c for complex and hyper complex
			//This means that exp surface for exp(xi) is tangent to 1 + 0i + ...
			//how does this relate to expansion of e^xi, in particular if exp(xi)'|_x=0 is tangent and exp(xj)'|_x=0 is tangent, is exp(x(ai + bj))'|_x=0 tangent for all a,b?
			//this would probably be the ideal form because initially the exp curve isn't increasing, but maybers it increases more down the line, probs not tho

			//still need a decent method for displaying looping spaces + handling looping rotation spaces
			//Just do the association thing normal algebra exp (not mov/rot), where surface point space has a grid, where each grid points to the closest exp

		sqrtSysPre = mt.Mult(varExp, varExp);

		invSys = mt.Mult(varExp, fixedExp);
		invSys[0].scalar -= 1;

		inverseMat.resize((mt.count+1)*(mt.count+1)*(mt.count+1));
		for (int y = 0; y < (mt.count+1); y++) //scalar var case
		{
			inverseMat[(y*(mt.count+1))+y] = 1;
		}
		for (int x = 1; x < (mt.count+1); x++) //variable
		{
			inverseMat[(x*(mt.count+1)*(mt.count+1))+(x*(mt.count+1))+0] = 1;
			for (int y = 0; y < (mt.count+1); y++) //equation
			{
				for (int i = 1; i < mt.count+1; i++) //iterate through fixed components
				{
					inverseMat[(x*(mt.count+1)*(mt.count+1))+(y*(mt.count+1))+i] += mt.At(x-1,i-1,y-1);
				}
			}
		}

		/*invExpSysPre.push_back(SymExp(1));
		invExpSysPre.push_back(SymExp(0));
		invExpSysPre.push_back(SymExp(0));
		invExpSysPre.push_back(SymExp(0));
		varMult = invExpSysPre;
		for (int i = 1; i < 6; i++) //10! ~ 2^22
		{
			varMult = mt.Mult(varMult, varExp);
			for (int j = 0; j < 4; j++)
			{
				varMult[j] *= Product(1.0f/i);
				invExpSysPre[j] += varMult[j];
			}
		}

		invExpDist = 0;
		for (float s = -0.5f; s <= 0.5f; s += 0.1f)
		{
			for (float i = -0.5f; i <= 0.5f; i += 0.1f)
			{
				for (float j = -0.5f; j <= 0.5f; j += 0.1f)
				{
					for (float k = -0.5f; k <= 0.5f; k += 0.1f)
					{
						QuatRep v(s, i, j, k);
						v *= 0.25f/std::sqrtf(v.SqrMag());
						QuatRep q = QuatExp(v);
						q.s -= 1;
						float sqrMag = q.SqrMag();
						if (sqrMag > invExpDist)
							invExpDist = sqrMag;
					}
				}
			}
		}
		*/

		d1 = QuatRep(0, 1, 0, 0);
		d2 = QuatRep(0, 0.2f, 1, 0);
		//d1 = QuatRep(0, -0.7071, 0.7071, 0);
		d1 = AntiScaleDir(d1);
		//d2 = QuatRep(0, 0.7071, 0.7071, 0);
		d2 = AntiScaleDir(d2);
		r1 = (d1*d2 - d2*d1)*0.5f;
		//r1 = QuatRep(0, 0, 0, 1);
		//r1 = QuatRep(0, 0, 0, 0);
		//r1 = AntiScaleDir(r1); pretty sure can't do this
		if (r1.ProjectOffOf(d1).ProjectOffOf(d2).SqrMag() < 0.0001f)
			r1 = QuatRep(0, 0, 0, 0);


		QuatRepDirs::d1 = d1;
		QuatRepDirs::d2 = d2;
		QuatRepDirs::r1 = r1;

		//somehow this isn't 100% like nonbasis change, so something is broken in calculations, especially with regards to rotation

		std::vector<SymExp> moveBase;
		moveBase.push_back(SymExp(Product(d1.s, 0, 1))); moveBase[0].terms.push_back(Product(d2.s, 1, 1)); moveBase[0].Simplify();
		moveBase.push_back(SymExp(Product(d1.i, 0, 1))); moveBase[1].terms.push_back(Product(d2.i, 1, 1)); moveBase[1].Simplify();
		moveBase.push_back(SymExp(Product(d1.j, 0, 1))); moveBase[2].terms.push_back(Product(d2.j, 1, 1)); moveBase[2].Simplify();
		moveBase.push_back(SymExp(Product(d1.k, 0, 1))); moveBase[3].terms.push_back(Product(d2.k, 1, 1)); moveBase[3].Simplify();
		std::vector<SymExp> rotBase;
		rotBase.push_back(SymExp(Product(r1.s, 2, 1))); rotBase[0].Simplify();
		rotBase.push_back(SymExp(Product(r1.i, 2, 1))); rotBase[1].Simplify();
		rotBase.push_back(SymExp(Product(r1.j, 2, 1))); rotBase[2].Simplify();
		rotBase.push_back(SymExp(Product(r1.k, 2, 1))); rotBase[3].Simplify();

		std::vector<SymExp> finalMovSys;
		finalMovSys.push_back(SymExp(1));
		finalMovSys.push_back(SymExp(0));
		finalMovSys.push_back(SymExp(0));
		finalMovSys.push_back(SymExp(0));
		std::vector<SymExp> varMult = finalMovSys;
		for (int i = 1; i < 8; i++) //10! ~ 2^22
		{
			varMult = mt.Mult(varMult, moveBase);
			for (int j = 0; j < 4; j++)
			{
				varMult[j] *= Product(1.0f/i);
				finalMovSys[j] += varMult[j];
			}
		}
		std::vector<SymExp> finalRotSys;
		finalRotSys.push_back(SymExp(1));
		finalRotSys.push_back(SymExp(0));
		finalRotSys.push_back(SymExp(0));
		finalRotSys.push_back(SymExp(0));
		varMult = finalRotSys;
		for (int i = 1; i < 8; i++) //10! ~ 2^22
		{
			varMult = mt.Mult(varMult, rotBase);
			for (int j = 0; j < 4; j++)
			{
				varMult[j] *= Product(1.0f/i);
				finalRotSys[j] += varMult[j];
			}
		}

		//gotta use a different norm Sys.
		//For h2, we need something like current for the rotation dir, but this breaks the movement dir
		//Ideally use normal to exp surface, given evaluated gradient, is easy to calculate, but would need to redo NM for this specific case
		//btw to calc normal, take basis vector, proj remove all known vectors. If basis ends up as 0 pick a different basis until a valid one is found
		//gradient with respect to i at z = sum( (z^n)*i/n! ). So calculating gradient and exp through mult is probs faster than FlatSymExp
		//FlatSymExp calc has complexity that increased polynomially with respect to iteration limit, whereas direct approach is always linear. Likewise direct is inheritantly faster, since iteration 1 they are equal, but Flat gets worse
		//Should implement a new version of NMSolve that takes in functions of vectors to vectors for eval and gradient, instead of SymExps
		

		std::vector<SymExp> normSys;
		normSys.push_back(SymExp(1)); normSys.push_back(SymExp(0)); normSys.push_back(SymExp(0)); normSys.push_back(SymExp(0));
		for (int i = 1; i < 8; i++)
			normSys[0].terms.push_back(Product(1.0f/i, 3, i));
		

		//invGeoExpSysPre = mt.Mult(finalRotSys, finalMovSys); //invGeoExpSysPre = mt.Mult(normSys, invGeoExpSysPre);
		invGeoExpRRSysPre = mt.Mult(finalMovSys, finalRotSys); for (int i = 0; i < 4; i++) invGeoExpRRSysPre[i].Simplify();  invGeoExpRRSysPre = mt.Mult(normSys, invGeoExpRRSysPre);
		
		std::vector<SymExp> offMoveBase = moveBase;
		offMoveBase[1].scalar = -0; 
		offMoveBase[2] = SymExp(0);
		finalMovSys[0] = SymExp(1); finalMovSys[1] = SymExp(0); finalMovSys[2] = SymExp(0); finalMovSys[3] = SymExp(0);

		varMult = finalMovSys;
		for (int i = 1; i < 8; i++) //10! ~ 2^22
		{
			varMult = mt.Mult(varMult, offMoveBase);
			for (int j = 0; j < 4; j++)
			{
				varMult[j] *= Product(1.0f/i);
				for (int i = 0; i < varMult[j].terms.size(); i++)
				{
					Product& term = varMult[j].terms[i];
					int totalPow = 0;
					for (int j = 0; j < term.pows.size(); j++)
						totalPow += term.pows[j];
					if (totalPow > 8 || std::abs(term.coeff) < 0.0000001f)
					{
						varMult[j].terms.erase(varMult[j].terms.begin()+i);
						i--;
					}
				}
				finalMovSys[j] += varMult[j];
			}
		}

		testInvGeoExpSysPre = mt.Mult(normSys,finalMovSys);

		for (int i = 0; i < 4; i++)
		{
			sqrtSysPre[i].Simplify();
			//invSys[i].Simplify();
			//invExpSysPre[i].Simplify();
			//invGeoExpSysPre[i].Simplify();
			invGeoExpRRSysPre[i].Simplify();
			testInvGeoExpSysPre[i].Simplify();
		}

		testInvGeoExpSysPre[3].terms.push_back(Product(0, 1, 1));
		testInvGeoExpSysPre[3].terms.push_back(Product(0, 2, 1));

		sqrtSys = FlatSymExp(sqrtSysPre);
		sqrtSysGrad = sqrtSys.Gradient();
		//invExpSys = FlatSymExp(invExpSysPre);
		//invExpSysGrad = invExpSys.Gradient();
		//invGeoExpSys = FlatSymExp(invGeoExpSysPre);
		//invGeoExpSysGrad = invGeoExpSys.Gradient();
		invGeoExpRRSys = FlatSymExp(invGeoExpRRSysPre);
		invGeoExpRRSysGrad = invGeoExpRRSys.Gradient();
		testInvGeoExpSys = FlatSymExp(testInvGeoExpSysPre);
		testInvGeoExpSysGrad = testInvGeoExpSys.Gradient();
	}
}; thread_local __SetupQuat_ __sq_;

QuatRep QuatRep::Sqrt()
{
	std::vector<float> fixed; fixed.resize(4);

	sqrtSys.expScalar(0) = -s;
	sqrtSys.expScalar(1) = -i;
	sqrtSys.expScalar(2) = -j;
	sqrtSys.expScalar(3) = -k;

	//use fixed as initial guess
	fixed[0] = s;
	fixed[1] = i+0.1f;
	fixed[2] = j;
	fixed[3] = k;
	fixed = sqrtSys.NewtonsMethodSolve(sqrtSysGrad, fixed);

	return QuatRep(fixed[0], fixed[1], fixed[2], fixed[3]);
}

QuatRep QuatRep::Inverse()
{

	//QuatRep a = *this; a.i *= -1; a.j *= -1; a.k *= -1;
	//return a;

	std::vector<float> fixed; fixed.resize(4);
	std::vector<int> ids; ids.resize(4);
	fixed[0] = s; ids[0] = 4; fixed[1] = i; ids[1] = 5; fixed[2] = j; ids[2] = 6; fixed[3] = k; ids[3] = 7; 

	std::vector<SymExp> solveInvSys; solveInvSys.resize(invSys.size());
	for (int i = 0; i < invSys.size(); i++)
		solveInvSys[i] = invSys[i].Eval(ids, fixed);

	//use fixed as initial guess
	ids[0] = 0; ids[1] = 1; ids[2] = 2; ids[3] = 3;
	//fixed = NMnTom(solveInvSys, ids, fixed);

	Vector2D<float> inverseMatEq(5, 4); inverseMatEq.At(4, 0) = 1;
	for (int x = 0; x < 4; x++)
		for (int y = 0; y < 4; y++)
			for (int f = 0; f < 4; f++)
				inverseMatEq.At(x, y) += inverseMat[(x*4*4)+(y*4)+f]*fixed[f];

	fixed = NMSolveLinearSystem(inverseMatEq);

	return QuatRep(fixed[0], fixed[1], fixed[2], fixed[3]);
}

/*QuatRep QuatRep::InvExp()
{
	QuatRep reduce = *this;
	int mult = 1;
	while (QuatRep(reduce.s-1,reduce.i,reduce.j,reduce.k).SqrMag() > invExpDist && mult < 256*256)
	{
		reduce = reduce.Sqrt();
		mult *= 2;
	}
	return reduce.DirectInvExp()*mult;
}*/

/*QuatRep QuatRep::DirectInvExp()
{
	std::vector<float> fixed; fixed.resize(4);

	invExpSys.expScalar(0) = 1-s;
	invExpSys.expScalar(1) = -i;
	invExpSys.expScalar(2) = -j;
	invExpSys.expScalar(3) = -k;

	//use fixed as initial guess
	for (int i = 0; i < 4; i++) fixed[i] = 0;
	fixed = invExpSys.NewtonsMethodSolve(invExpSysGrad, fixed);

	return QuatRep(fixed[0], fixed[1], fixed[2], fixed[3]);
}*/

QuatRep QuatRep::InvGeoExpLR()
{
	return InvGeoExpLR(QuatRep(0, 0, 0, 0));
}

QuatRep QuatRep::InvGeoExpLR(QuatRep initial, int maxCount)
{
	//use fixed as initial guess
	std::vector<float> fixed; fixed.resize(3);
	fixed[0] = initial.i; fixed[1] = initial.j; fixed[2] = initial.k; // fixed[3] = 0;
	fixed.resize(4);

	std::vector<float> target; target.resize(4); target[0] = s; target[1] = i; target[2] = j; target[3] = k;

	fixed = NMnTomFunction(4, QuatGeoExpEval, QuatGeoExpGradient, fixed, target,maxCount);
	fixed[3] = 0;

	return QuatRep(0, fixed[0], fixed[1], fixed[2]);
}

QuatRep QuatRep::InvGeoExp()
{
	return InvGeoExp(QuatRep(0, 0, 0, 0));
}

QuatRep QuatRep::InvGeoExp(QuatRep initial, int maxCount)
{
	if (!(std::abs(initial.s) < 8 && std::abs(initial.i) < 8 && std::abs(initial.j) < 8 && std::abs(initial.k) < 8))
		initial = QuatRep(0, 0, 0, 0);

	invGeoExpRRSys.expScalar(0) = 1-s;
	invGeoExpRRSys.expScalar(1) = -i;
	invGeoExpRRSys.expScalar(2) = -j;
	invGeoExpRRSys.expScalar(3) = -k;

	//use fixed as initial guess
	std::vector<float> fixed; fixed.resize(3);
	fixed[0] = initial.i; fixed[1] = initial.j; fixed[2] = initial.k; // fixed[3] = 0;
	fixed.resize(4);

	std::vector<float> target; target.resize(4); target[0] = s; target[1] = i; target[2] = j; target[3] = k;

	if (initial.i < -0.0f && false)
	{
		fixed.resize(4);
		fixed[0] += 0;
		testInvGeoExpSys.expScalar(0) = 1-s;
		testInvGeoExpSys.expScalar(1) = -i;
		testInvGeoExpSys.expScalar(2) = -j;
		testInvGeoExpSys.expScalar(3) = -k;
		fixed = testInvGeoExpSys.NewtonsMethodSolve(testInvGeoExpSysGrad, fixed, maxCount);
		auto testOut = testInvGeoExpSys.SclEval(fixed);
		fixed[0] -= 0;
	}
	else
		//fixed = invGeoExpRRSys.NewtonsMethodSolve(invGeoExpRRSysGrad, fixed, maxCount);
		fixed = NMnTomFunction(4, QuatGeoExpEval, QuatGeoExpGradient, fixed, target,maxCount);
	fixed[3] = 0;

	//add in some sort of quotient solver later, e.g. generate set of directions, then figure out when going forward in that direction collides with going backwards in that direction.
	//then given a point, find it's closest direction, then modulo by the length of that direction
	//in this case, we would know that the k direction repeats at 3.141
	
	//if (fixed[2] >= 3.14159f)
	//	fixed[2] -= 3.14159f*2.0f;
	//if (fixed[2] <= -3.14159f)
	//	fixed[2] += 3.14159f*2.0f;

	/*float mag = glm::sqrt((fixed[0]*fixed[0])+(fixed[1]*fixed[1]));
	float a = mag;
	if (mag >= 3.14159f/2.0f)
	{
		a = mag - 3.14159f;
		fixed[0] *= a/mag;
		fixed[1] *= a/mag;
	}*/

	return QuatRep(0, fixed[0], fixed[1], fixed[2]);
}

//should precompute rot, but this function only gets called once anyways
QuatRep QuatRep::LocalMove(float x, float y) //tries to use standardized distance based on how R*T (moving after rotation) compares to T*R (quotiented move)
{
	QuatRep decomp = (*this).InvGeoExp();
	QuatRep rot = QuatRotExp(decomp.k);

	QuatRep unstandard = rot*((d1*x)+(d2*y));
	float mag = (x*x)+(y*y);
	if (mag == 0)
		return (*this);
	//we know R' is close to R at least
	//maybe even assume R' = R

	Vector2D<float> grad(3,4);

	QuatRep iGrad = d1*rot;
	QuatRep jGrad = d2*rot;

	//cursed
	grad.At(0, 0) = iGrad.s; grad.At(1, 0) = jGrad.s; grad.At(2, 0) = unstandard.s;
	grad.At(0, 1) = iGrad.i; grad.At(1, 1) = jGrad.i; grad.At(2, 1) = unstandard.i;
	grad.At(0, 2) = iGrad.j; grad.At(1, 2) = jGrad.j; grad.At(2, 2) = unstandard.j;
	grad.At(0, 3) = iGrad.k; grad.At(1, 3) = jGrad.k; grad.At(2, 3) = unstandard.k;

	std::vector<float> standardMove = NMSolveLinearSystem(grad);
	float standardMag = (standardMove[0]*standardMove[0])+(standardMove[1]*standardMove[1]);
	float coeff = glm::sqrt(mag/standardMag); //standardMag is the real mag, so we divide by it then multiply by mag to get movement with the intended mag
	return (*this)*QuatGeoExp(x,y,0);
}

QuatRep QuatRep::MoveAdjust(QuatRep camR)
{
	QuatRep rot = camR;

	QuatRep unstandard = rot*((d1*i)+(d2*j));
	float mag = (i*i)+(j*j);
	if (mag == 0)
		return (*this);
	//we know R' is close to R at least
	//maybe even assume R' = R

	Vector2D<float> grad(3,4);

	QuatRep iGrad = d1*rot;
	QuatRep jGrad = d2*rot;

	//cursed
	grad.At(0, 0) = iGrad.s; grad.At(1, 0) = jGrad.s; grad.At(2, 0) = unstandard.s;
	grad.At(0, 1) = iGrad.i; grad.At(1, 1) = jGrad.i; grad.At(2, 1) = unstandard.i;
	grad.At(0, 2) = iGrad.j; grad.At(1, 2) = jGrad.j; grad.At(2, 2) = unstandard.j;
	grad.At(0, 3) = iGrad.k; grad.At(1, 3) = jGrad.k; grad.At(2, 3) = unstandard.k;

	std::vector<float> standardMove = NMSolveLinearSystem(grad);
	float standardMag = (standardMove[0]*standardMove[0])+(standardMove[1]*standardMove[1]);
	float coeff = glm::sqrt(standardMag/mag); //inverse of move, cause we are dividing by coeff of real move
	return QuatRep(s, i*coeff, j*coeff, k);
}

float QuatRep::SqrMag()
{
	return (s*s)+(i*i)+(j*j)+(k*k);
}

float QuatRep::Dot(QuatRep q)
{
	return (s*q.s)+(i*q.i)+(j*q.j)+(k*q.k);
}

QuatRep QuatRep::ProjectOnTo(QuatRep q)
{
	return q*(Dot(q)/q.Dot(q));
}

QuatRep QuatRep::ProjectOffOf(QuatRep q)
{
	return (*this)-ProjectOnTo(q);
}

std::string QuatRep::ToString()
{
	return std::string("s: ") + std::to_string(s) + " i: " + std::to_string(i) + " j: " + std::to_string(j) + " k: " + std::to_string(k);
}

QuatRep QuatExp(QuatRep quat)
{
	QuatRep outHold;
	QuatRep curHold;
	for (int n = 1; n <= 32; n++)
	{
		curHold *= quat;
		curHold /= n;
		outHold += curHold;
	}

	return outHold;
}
QuatRep QuatExp(float s, float i, float j, float k)
{
	return QuatExp(QuatRep(s, i, j, k));
}

QuatRep QuatGeoExp(float x, float y, float r)
{
	QuatRep quat = r1*r;
	QuatRep rotHold;
	QuatRep outHold;
	QuatRep curHold;

	for (int n = 1; n <= 32; n++)
	{
		curHold *= quat;
		curHold /= n;
		rotHold += curHold;
	}

	quat = (d1*x)+(d2*y);
	curHold = QuatRep();
	for (int n = 1; n <= 32; n++)
	{
		curHold *= quat;
		curHold /= n;
		outHold += curHold;
	}

	return rotHold*outHold;
}
QuatRep QuatGeoExp(QuatRep v)
{
	return QuatGeoExp(v.i,v.j,v.k);
}

void QuatGeoExpEval(std::vector<float>& in, std::vector<float>& out)
{
	in.resize(4); out.resize(4);

	QuatRep moveQ = d1*in[0] + d2*in[1];
	QuatRep rotQ = r1*in[2];

	QuatRep move = QuatMovExp(in[0], in[1]);
	QuatRep rot = QuatRotExp(in[2]);
	
	QuatRep outHold = move*rot;

	out[0] = outHold.s; out[1] = outHold.i;
	out[2] = outHold.j; out[3] = outHold.k;
}

void QuatGeoExpGradient(std::vector<float>& in, std::vector<float>& out)
{
	in.resize(4);
	out.resize(4*4);

	QuatRep moveQ = d1*in[0] + d2*in[1];
	QuatRep rotQ = r1*in[2];
	//der is (zi + iz)/2! + (zzi + ziz + izz)/3! + ...
	QuatRep moveHold[16];
	moveHold[0] = QuatRep();
	for (int i = 1; i < 16; i++)
		moveHold[i] = (moveHold[i-1]*moveQ);

	QuatRep rotHold[16];
	rotHold[0] = QuatRep();
	for (int i = 1; i < 16; i++)
		rotHold[i] = (rotHold[i-1]*rotQ);

	QuatRep move = QuatMovExp(in[0], in[1]);
	QuatRep rot = QuatRotExp(in[2]);
	
	QuatRep iGrad = QuatIDer(VectorRef<QuatRep>(16,moveHold), rot);
	QuatRep jGrad = QuatJDer(VectorRef<QuatRep>(16,moveHold), rot);
	QuatRep kGrad = QuatKDer(VectorRef<QuatRep>(16,rotHold), move);

	out[0] = iGrad.s; out[1] = iGrad.i; out[2] = iGrad.j; out[3] = iGrad.k;
	out[4] = jGrad.s; out[5] = jGrad.i; out[6] = jGrad.j; out[7] = jGrad.k;
	out[8] = kGrad.s; out[9] = kGrad.i; out[10] = kGrad.j; out[11] = kGrad.k;

	//iGrad = iGrad;
	jGrad = jGrad.ProjectOffOf(iGrad); 
	kGrad = kGrad.ProjectOffOf(iGrad).ProjectOffOf(jGrad);

	glm::vec4 norm(1,1,1,1);
	norm -= glm::vec4(iGrad.s,iGrad.i,iGrad.j,iGrad.k)*(glm::dot(norm,glm::vec4(iGrad.s,iGrad.i,iGrad.j,iGrad.k)/iGrad.SqrMag()));
	norm -= glm::vec4(jGrad.s,jGrad.i,jGrad.j,jGrad.k)*(glm::dot(norm,glm::vec4(jGrad.s,jGrad.i,jGrad.j,jGrad.k)/jGrad.SqrMag()));
	norm -= glm::vec4(kGrad.s,kGrad.i,kGrad.j,kGrad.k)*(glm::dot(norm,glm::vec4(kGrad.s,kGrad.i,kGrad.j,kGrad.k)/kGrad.SqrMag()));
	if (glm::dot(norm, norm) < 0.00001f)
	{
		norm = glm::vec4(1.1f, 1, 1, 1);
		norm -= glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)*(glm::dot(norm, glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)/iGrad.SqrMag()));
		norm -= glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)*(glm::dot(norm, glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)/jGrad.SqrMag()));
		norm -= glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)*(glm::dot(norm, glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)/kGrad.SqrMag()));
	}
	if (glm::dot(norm, norm) < 0.00001f)
	{
		norm = glm::vec4(1, 1.1f, 1, 1);
		norm -= glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)*(glm::dot(norm, glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)/iGrad.SqrMag()));
		norm -= glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)*(glm::dot(norm, glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)/jGrad.SqrMag()));
		norm -= glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)*(glm::dot(norm, glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)/kGrad.SqrMag()));
	}
	if (glm::dot(norm, norm) < 0.00001f)
	{
		norm = glm::vec4(1, 1, 1.1f, 1);
		norm -= glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)*(glm::dot(norm, glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)/iGrad.SqrMag()));
		norm -= glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)*(glm::dot(norm, glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)/jGrad.SqrMag()));
		norm -= glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)*(glm::dot(norm, glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)/kGrad.SqrMag()));
	}

	out[12] = norm.x; out[13] = norm.y; out[14] = norm.z; out[15] = norm.w;
	for (int i = 12; i < 16; i++) out[i] = 0;
}

void QuatGeoExpLREval(std::vector<float>& in, std::vector<float>& out)
{
	in.resize(4); out.resize(4);

	QuatRep moveQ = d1*in[0] + d2*in[1];
	QuatRep rotQ = r1*in[2];

	QuatRep move = QuatMovExp(in[0], in[1]);
	QuatRep rot = QuatRotExp(in[2]);

	QuatRep outHold = rot*move;

	out[0] = outHold.s; out[1] = outHold.i;
	out[2] = outHold.j; out[3] = outHold.k;
}

void QuatGeoExpLRGradient(std::vector<float>& in, std::vector<float>& out)
{
	in.resize(4);
	out.resize(4*4);

	QuatRep moveQ = d1*in[0] + d2*in[1];
	QuatRep rotQ = r1*in[2];
	//der is (zi + iz)/2! + (zzi + ziz + izz)/3! + ...
	QuatRep moveHold[16];
	moveHold[0] = QuatRep();
	for (int i = 1; i < 16; i++)
		moveHold[i] = (moveHold[i-1]*moveQ);

	QuatRep rotHold[16];
	rotHold[0] = QuatRep();
	for (int i = 1; i < 16; i++)
		rotHold[i] = (rotHold[i-1]*rotQ);

	QuatRep move = QuatMovExp(in[0], in[1]);
	QuatRep rot = QuatRotExp(in[2]);

	QuatRep iGrad = QuatIDerLR(VectorRef<QuatRep>(16,moveHold), rot);
	QuatRep jGrad = QuatJDerLR(VectorRef<QuatRep>(16,moveHold), rot);
	QuatRep kGrad = QuatKDerLR(VectorRef<QuatRep>(16,rotHold), move);

	out[0] = iGrad.s; out[1] = iGrad.i; out[2] = iGrad.j; out[3] = iGrad.k;
	out[4] = jGrad.s; out[5] = jGrad.i; out[6] = jGrad.j; out[7] = jGrad.k;
	out[8] = kGrad.s; out[9] = kGrad.i; out[10] = kGrad.j; out[11] = kGrad.k;

	//iGrad = iGrad;
	jGrad = jGrad.ProjectOffOf(iGrad); 
	kGrad = kGrad.ProjectOffOf(iGrad).ProjectOffOf(jGrad);

	glm::vec4 norm(1,1,1,1);
	norm -= glm::vec4(iGrad.s,iGrad.i,iGrad.j,iGrad.k)*(glm::dot(norm,glm::vec4(iGrad.s,iGrad.i,iGrad.j,iGrad.k)/iGrad.SqrMag()));
	norm -= glm::vec4(jGrad.s,jGrad.i,jGrad.j,jGrad.k)*(glm::dot(norm,glm::vec4(jGrad.s,jGrad.i,jGrad.j,jGrad.k)/jGrad.SqrMag()));
	norm -= glm::vec4(kGrad.s,kGrad.i,kGrad.j,kGrad.k)*(glm::dot(norm,glm::vec4(kGrad.s,kGrad.i,kGrad.j,kGrad.k)/kGrad.SqrMag()));
	if (glm::dot(norm, norm) < 0.00001f)
	{
		norm = glm::vec4(1.1f, 1, 1, 1);
		norm -= glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)*(glm::dot(norm, glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)/iGrad.SqrMag()));
		norm -= glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)*(glm::dot(norm, glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)/jGrad.SqrMag()));
		norm -= glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)*(glm::dot(norm, glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)/kGrad.SqrMag()));
	}
	if (glm::dot(norm, norm) < 0.00001f)
	{
		norm = glm::vec4(1, 1.1f, 1, 1);
		norm -= glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)*(glm::dot(norm, glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)/iGrad.SqrMag()));
		norm -= glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)*(glm::dot(norm, glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)/jGrad.SqrMag()));
		norm -= glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)*(glm::dot(norm, glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)/kGrad.SqrMag()));
	}
	if (glm::dot(norm, norm) < 0.00001f)
	{
		norm = glm::vec4(1, 1, 1.1f, 1);
		norm -= glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)*(glm::dot(norm, glm::vec4(iGrad.s, iGrad.i, iGrad.j, iGrad.k)/iGrad.SqrMag()));
		norm -= glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)*(glm::dot(norm, glm::vec4(jGrad.s, jGrad.i, jGrad.j, jGrad.k)/jGrad.SqrMag()));
		norm -= glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)*(glm::dot(norm, glm::vec4(kGrad.s, kGrad.i, kGrad.j, kGrad.k)/kGrad.SqrMag()));
	}

	out[12] = norm.x; out[13] = norm.y; out[14] = norm.z; out[15] = norm.w;
	for (int i = 12; i < 16; i++) out[i] = 0;
}

QuatRep QuatMovExp(float x, float y)
{
	QuatRep hold(1,0,0,0);
	QuatRep cur(1,0,0,0);
	for (int i = 1; i < 16; i++)
	{
		cur *= d1*(x/i) + d2*(y/i);
		hold += cur;
	}
	return hold;
}

QuatRep QuatRotExp(float r)
{
	QuatRep hold(1,0,0,0);
	QuatRep cur(1,0,0,0);
	for (int i = 1; i < 16; i++)
	{
		cur *= r1*(r/i);
		hold += cur;
	}
	return hold;
}

QuatRep QuatIDer(VectorRef<QuatRep> moveZ, QuatRep rot)
{
	QuatRep hold(0,0,0,0);
	float denom = 1;
	for (int i = 0; i < moveZ.size(); i++)
	{
		denom *= i+1;
		for (int j = 0; j <= i; j++)
		{
			hold += (moveZ[j]*d1*moveZ[i-j])/denom;
		}
	}
	return hold*rot;
}
QuatRep QuatJDer(VectorRef<QuatRep> moveZ, QuatRep rot)
{
	QuatRep hold(0,0,0,0);
	float denom = 1;
	for (int i = 0; i < moveZ.size(); i++)
	{
		denom *= i+1;
		for (int j = 0; j <= i; j++)
		{
			hold += (moveZ[j]*d2*moveZ[i-j])/denom;
		}
	}
	return hold*rot;
}
QuatRep QuatKDer(VectorRef<QuatRep> rotZ, QuatRep move)
{
	QuatRep hold(0,0,0,0);
	float denom = 1;
	for (int i = 0; i < rotZ.size(); i++)
	{
		denom *= i+1;
		for (int j = 0; j <= i; j++)
		{
			hold += (rotZ[j]*r1*rotZ[i-j])/denom;
		}
	}
	return move*hold;
}

QuatRep QuatIDerLR(VectorRef<QuatRep> moveZ, QuatRep rot)
{
	QuatRep hold(0,0,0,0);
	float denom = 1;
	for (int i = 0; i < moveZ.size(); i++)
	{
		denom *= i+1;
		for (int j = 0; j <= i; j++)
		{
			hold += (moveZ[j]*d1*moveZ[i-j])/denom;
		}
	}
	return rot*hold;
}
QuatRep QuatJDerLR(VectorRef<QuatRep> moveZ, QuatRep rot)
{
	QuatRep hold(0,0,0,0);
	float denom = 1;
	for (int i = 0; i < moveZ.size(); i++)
	{
		denom *= i+1;
		for (int j = 0; j <= i; j++)
		{
			hold += (moveZ[j]*d2*moveZ[i-j])/denom;
		}
	}
	return rot*hold;
}
QuatRep QuatKDerLR(VectorRef<QuatRep> rotZ, QuatRep move)
{
	QuatRep hold(0,0,0,0);
	float denom = 1;
	for (int i = 0; i < rotZ.size(); i++)
	{
		denom *= i+1;
		for (int j = 0; j <= i; j++)
		{
			hold += (rotZ[j]*r1*rotZ[i-j])/denom;
		}
	}
	return hold*move;
}

float TestMag(QuatRep dir, float dist)
{

	float expMag1;
	float expMag2;
	int sign1 = 1;
	int sign2 = -1;

	//calculate mag by projecting to dir, then squaring that and the other components
	QuatRep test = QuatExp(dir*dist);
	float testProj = dir.Dot(test)/dir.Dot(dir);
	expMag1 = (testProj*testProj)+(test.s*test.s)+(test.j*test.j)+(test.k*test.k);
	/*if (expMag1 < 1.0f)
	{
	sign1 = -1;
	expMag1 = 1.0f/expMag1;
	}*/

	test = QuatExp(dir*-dist);
	testProj = dir.Dot(test)/dir.Dot(dir);
	expMag2 = (testProj*testProj)+(test.s*test.s)+(test.j*test.j)+(test.k*test.k);
	/*if (expMag2 < 1.0f)
	{
	sign2 = 1;
	expMag2 = 1.0f/expMag2;
	}*/

	if (expMag2 > expMag1)
	{
		expMag1 = expMag2;
		sign1 = sign2;
	}
	expMag1 = glm::sqrt(expMag1);

	return expMag1;
}

QuatRep AntiScaleDir(QuatRep dir)
{
	QuatRep i2 = dir*dir;
	QuatRep dirProj = dir; dirProj.s = 0;
	

	//square it => project it by checking i coeff versus i coeff of dir => get scaler from i/2
	//e.g. if dir = 2i and square has i coeff = 6i, then basis convert would be 3i - some scalar
	dir.s -= i2.Dot(dirProj)/(2*dirProj.Dot(dirProj));

	if (isnan(dir.s))
		dir.s = 0;
	return dir;

	float dist = 9;
	float coeff = 1;

	for (int i = 0; i < 80; i++)
	{
		float expC = TestMag(dir, dist);
		float expL = TestMag(dir-QuatRep(coeff, 0, 0, 0), dist);
		float expR = TestMag(dir+QuatRep(coeff, 0, 0, 0), dist);

		if (expR < expL && expR < expC)
			dir = dir+QuatRep(coeff, 0, 0, 0);
		else if (expL < expC)
			dir = dir-QuatRep(coeff, 0, 0, 0);
		else
			coeff /= 2;

		if (coeff <= 0.000001)
			break;

		//float test2 = (glm::log(expMag1)/dist)/coeff;
		//dir = dir - QuatRep(sign1*glm::log(expMag1)/(coeff*dist), 0, 0, 0);

	}
	return dir;
}
