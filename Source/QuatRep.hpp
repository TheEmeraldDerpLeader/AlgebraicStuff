#pragma once

#include <string>
#include <vector>
#include <Vector.hpp>

struct QuatRep
{
public:
	float s = 1;
	float i = 0;
	float j = 0;
	float k = 0;

	QuatRep() {}
	QuatRep(float sF, float iF, float jF, float kF) : s(sF), i(iF), j(jF), k(kF) {}

	QuatRep operator*(QuatRep quat);
	QuatRep operator/(QuatRep quat);
	QuatRep operator+(QuatRep quat);
	QuatRep operator-(QuatRep quat);

	QuatRep& operator*=(QuatRep quat);
	QuatRep& operator/=(QuatRep quat);
	QuatRep& operator+=(QuatRep quat);
	QuatRep& operator-=(QuatRep quat);

	QuatRep operator*(float scale);
	QuatRep operator/(float scale);
	QuatRep& operator*=(float scale);
	QuatRep& operator/=(float scale);

	QuatRep Sqrt();
	QuatRep Inverse();
	//QuatRep InvExp();
	//QuatRep DirectInvExp();

	QuatRep InvGeoExpLR();
	QuatRep InvGeoExpLR(QuatRep initial, int maxCount = 80);
	QuatRep InvGeoExp(); //Gives InvGeoExp for rotation on the right side
	QuatRep InvGeoExp(QuatRep initial, int maxCount = 80);

	QuatRep InvTestExp();

	QuatRep LocalMove(float x, float y);
	QuatRep MoveAdjust(QuatRep cam); //converts (*this) as a vector in exp dist to local move dist

	float SqrMag();
	float Dot(QuatRep q);
	QuatRep ProjectOnTo(QuatRep q);
	QuatRep ProjectOffOf(QuatRep q);

	std::string ToString();
};

namespace QuatRepDirs
{
	extern thread_local QuatRep d1;
	extern thread_local QuatRep d2;
	extern thread_local QuatRep r1;
};


QuatRep QuatExp(float s, float i, float j, float k);
QuatRep QuatExp(QuatRep quat);

QuatRep QuatGeoExp(float x, float y, float r);
QuatRep QuatGeoExp(QuatRep v); //don't have float3, so just use last 3 components

void QuatGeoExpEval(std::vector<float>& in, std::vector<float>& out);
void QuatGeoExpGradient(std::vector<float>& in, std::vector<float>& out);
void QuatGeoExpLREval(std::vector<float>& in, std::vector<float>& out);
void QuatGeoExpLRGradient(std::vector<float>& in, std::vector<float>& out);


QuatRep QuatMovExp(float x, float y);
QuatRep QuatRotExp(float r);
QuatRep QuatIDer(VectorRef<QuatRep> moveZ, QuatRep rot);
QuatRep QuatJDer(VectorRef<QuatRep> moveZ, QuatRep rot);
QuatRep QuatKDer(VectorRef<QuatRep> rotZ, QuatRep move);

QuatRep QuatIDerLR(VectorRef<QuatRep> moveZ, QuatRep rot);
QuatRep QuatJDerLR(VectorRef<QuatRep> moveZ, QuatRep rot);
QuatRep QuatKDerLR(VectorRef<QuatRep> rotZ, QuatRep move);

QuatRep AntiScaleDir(QuatRep dir);