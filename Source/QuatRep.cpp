#include "QuatRep.hpp"

#include <SymExp.hpp>
#include <FlatSymExp.hpp>
#include <ArbAlgRep.hpp>

QuatRep QuatRep::operator*(QuatRep quat)
{
	QuatRep hold;
	hold.s = (s*quat.s)+(i*quat.i)+(j*quat.j)-(k*quat.k);
	hold.i = (s*quat.i)+(i*quat.s)+(j*quat.k)-(k*quat.j);
	hold.j = (s*quat.j)-(i*quat.k)+(j*quat.s)+(k*quat.i);
	hold.k = (s*quat.k)-(i*quat.j)+(j*quat.i)+(k*quat.s);
	//hold.s = (s*quat.s)-(i*quat.i)-(j*quat.j)-(k*quat.k);
	//hold.i = (s*quat.i)+(i*quat.s)+(j*quat.k)-(k*quat.j);
	//hold.j = (s*quat.j)-(i*quat.k)+(j*quat.s)+(k*quat.i);
	//hold.k = (s*quat.k)+(i*quat.j)-(j*quat.i)+(k*quat.s);
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
		mt.At(0, 0, -1) = 1; mt.At(0, 1, 2) = -1; mt.At(0, 2, 1) = -1; 
		mt.At(1, 1, -1) = 1; mt.At(1, 0, 2) = 1; mt.At(1, 2, 0) = 1;
		mt.At(2, 2, -1) = -1; mt.At(2, 0, 1) = 1; mt.At(2, 1, 0) = -1;

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
		d2 = QuatRep(0, 0, 1, 0);
		r1 = (d1*d2 - d2*d1)*0.5f;

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

/*QuatRep QuatRep::InvGeoExp()
{
	return InvGeoExp(QuatRep(0, 0, 0, 0));
}*/

/*QuatRep QuatRep::InvGeoExp(QuatRep initial, int maxCount)
{
	if (!(std::abs(initial.s) < 8.0f && std::abs(initial.i) < 8.0f && std::abs(initial.j) < 8.0f && std::abs(initial.k) < 8.0f))
		initial = QuatRep(0, 0, 0, 0);

	invGeoExpSys.expScalar(0) = 1-s;
	invGeoExpSys.expScalar(1) = -i;
	invGeoExpSys.expScalar(2) = -j;
	invGeoExpSys.expScalar(3) = -k;

	//use fixed as initial guess
	std::vector<float> fixed; fixed.resize(3);
	fixed[0] = initial.i; fixed[1] = initial.j; fixed[2] = initial.k;// fixed[3] = 0;
	fixed = invGeoExpSys.NewtonsMethodSolve(invGeoExpSysGrad, fixed, maxCount);

	return QuatRep(0, fixed[0], fixed[1], fixed[2]);
}*/

QuatRep QuatRep::InvGeoExpRR()
{
	return InvGeoExpRR(QuatRep(0, 0, 0, 0));
}

QuatRep QuatRep::InvGeoExpRR(QuatRep initial, int maxCount)
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

QuatRep QuatRep::LocalMove(float x, float y)
{
	QuatRep decomp = InvGeoExpRR(QuatRep(0, 0, 0, 0), 80);
	return QuatGeoExp(decomp.i,decomp.j,0)*QuatGeoExp(x,y,0)*QuatGeoExp(0,0,decomp.k);
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
