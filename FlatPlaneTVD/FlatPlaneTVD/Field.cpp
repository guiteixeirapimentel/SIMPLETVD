#include "Field.h"
#include <fstream>

Field::Field(double L, double H, double W, double deltaSize, double rho, double mu, double Ufarfield, double Vfarfield, double Wfarfield)
	:
	cdd(deltaSize),
	cdA(cdd*cdd),
	cdV(cdd*cdd*cdd),
	cNY(size_t(W/cdd) + 2UL),
	cNX(size_t(L/cdd) + 2UL),
	cNZ(size_t(H/cdd) + 2UL),
	cRHO(rho),
	cMU(mu),
	cD(cdA * mu / cdd)
{
	cPressureField.resize(cNY * cNZ * cNX, 0.0f);

	cPressureCorrField.resize(cNY * cNZ * cNX, 0.0f);

	cUVelField.resize(cNY * cNZ * cNX, Ufarfield);

	cVVelField.resize(cNY * cNZ * cNX, Vfarfield);
	
	cWVelField.resize(cNY * cNZ * cNX, Wfarfield);
	
	cSpUField.resize(cNY * cNZ * cNX, 0.0f);

	cSuUField.resize(cNY * cNZ * cNX, 0.0f);


	cSpVField.resize(cNY * cNZ * cNX, 0.0f);

	cSuVField.resize(cNY * cNZ * cNX, 0.0f);


	cSpWField.resize(cNY * cNZ * cNX, 0.0f);

	cSuWField.resize(cNY * cNZ * cNX, 0.0f);

	caijkU = cSuWField;
	caijkV = cSuWField;
	caijkW = cSuWField;
}

Field::~Field()
{}

Field::Intervalo::Intervalo(size_t v1, size_t v2)
	:
	valor1(v1),
	valor2(v2)
{}

size_t Field::Intervalo::Tamanho() const
{
	return valor2 - valor1;
}

void Field::SetBCFlatPlate()
{
	//Face N# 2
	SetUValueiJ(0, 0.0f);
	SetVValueIj(0, 0.0f);

	SetWValueIJ(0, 0.0f);
	SetWValueIJ(1, 0.0f);

	SetSpUValueInterval({ 2, cNX - 1 }, { 1, cNY - 1 }, { 1, 2 }, -cMU * cdA / (cdd / 2.0f));
	SetSpVValueInterval({ 1, cNX - 1 }, { 2, cNY - 1 }, { 1, 2 }, -cMU * cdA / (cdd / 2.0f));

	// Face N# 1 
	//SetUValueiJ(cNZ - 1, 0.0f);
	//SetVValueIj(cNZ - 1, 0.0f);
	//
	//SetWValueIJ(cNZ - 1, 0.0f);
	//
	//SetSpUValueInterval({ 2, cNX - 1 }, { 1, cNY - 1 }, { cNZ - 2, cNZ - 1 }, -cMU * cdA / (cdd / 2.0f));
	//SetSpVValueInterval({ 1, cNX - 1 }, { 2, cNY - 1 }, { cNZ - 2, cNZ - 1 }, -cMU * cdA / (cdd / 2.0f));

	// Face N# 6
	//SetUValueiK(0, 0.0f);
	//
	//SetVValueIK(0, 0.0f);
	//SetVValueIK(1, 0.0f);
	//
	//SetWValueIk(0, 0.0f);
	//
	//SetSpUValueInterval({ 2, cNX - 1 }, { 1, 2 }, { 1, cNZ - 1 }, -cMU * cdA / (cdd / 2.0f));
	//SetSpWValueInterval({ 1, cNX - 1 }, { 1, 2 }, { 2, cNZ - 1 }, -cMU * cdA / (cdd / 2.0f));

	// Face N# 7
	//SetUValueiK(cNY - 1, 0.0f);
	//
	//SetVValueIK(cNY - 1, 0.0f);
	//
	//SetWValueIk(cNY - 1, 0.0f);
	//
	//SetSpUValueInterval({ 2, cNX - 1 }, { cNY - 2, cNY - 1 }, { 1, cNZ - 1 }, -cMU * cdA / (cdd / 2.0f));
	//SetSpWValueInterval({ 1, cNX - 1 }, { cNY - 2, cNY - 1 }, { 2, cNZ - 1 }, -cMU * cdA / (cdd / 2.0f));
	//
	//// corners
	//SetSpUValueInterval({ 2, cNX - 1 }, { cNY - 2, cNY - 1 }, { cNZ - 2, cNZ - 1 }, -2.0f*cMU * cdA / (cdd / 2.0f));
	//SetSpUValueInterval({ 2, cNX - 1 }, { cNY - 2, cNY - 1 }, { 1, 2 }, -2.0f * cMU * cdA / (cdd / 2.0f));
	//
	//SetSpUValueInterval({ 2, cNX - 1 }, { 1, 2 }, { cNZ - 2, cNZ - 1 }, -2.0f * cMU * cdA / (cdd / 2.0f));
	//SetSpUValueInterval({ 2, cNX - 1 }, { 1, 2 }, { 1, 2 }, -2.0f * cMU * cdA / (cdd / 2.0f));


}

int Field::CreateUMomentumLSCSR(const doubleField3D& uField0, double dt,
	std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val, std::vector<double>& rhs)
{
	int n2 = cNX * cNY * cNZ;        // Number of points in the grid.
		

	ptr.clear();
	if (ptr.capacity() != n2 + 1)
		ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	if (col.capacity() != n2 * 7)
		col.reserve(n2 * 7);

	val.clear();
	if (val.capacity() != n2 * 7)
		val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);
	
	for (int j = 0, index = 0; j < cNY; ++j) 
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++, index++) {

				if (i == 1 || i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

					col.push_back(index);
					val.push_back(1.0f);

					rhs[index] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i];
				}
				else 
				{
					CoefsData aijkU;

					const double Fw = GetFwInternalUMomentum(i, j, k, cRHO) * cdA;
					const double Fe = GetFeInternalUMomentum(i, j, k, cRHO) * cdA;
					const double Fs = GetFsInternalUMomentum(i, j, k, cRHO) * cdA;
					const double Fn = GetFnInternalUMomentum(i, j, k, cRHO) * cdA;
					const double Fb = GetFbInternalUMomentum(i, j, k, cRHO) * cdA;
					const double Ft = GetFtInternalUMomentum(i, j, k, cRHO) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkU.aimjk = cD + fmaxf(Fw, 0.0f);
					aijkU.aipjk = cD + fmaxf(-Fe, 0.0f);

					aijkU.aijmk = cD + fmaxf(Fs, 0.0f);
					aijkU.aijpk = cD + fmaxf(-Fn, 0.0f);

					aijkU.aijkm = cD + fmaxf(Fb, 0.0f);
					aijkU.aijkp = cD + fmaxf(-Ft, 0.0f);

					const double alfaw = (Fw > 0.0f) ? 1.0f : 0.0f;
					const double alfae = (Fe > 0.0f) ? 1.0f : 0.0f;
					const double alfas = (Fs > 0.0f) ? 1.0f : 0.0f;
					const double alfan = (Fn > 0.0f) ? 1.0f : 0.0f;
					const double alfab = (Fb > 0.0f) ? 1.0f : 0.0f;
					const double alfat = (Ft > 0.0f) ? 1.0f : 0.0f;

					const double reneg = ren(i, j, k, cUVelField);
					const double rwneg = rwn(i, j, k, cUVelField);
					const double rnneg = rnn(i, j, k, cUVelField);
					const double rsneg = rsn(i, j, k, cUVelField);
					const double rtneg = rtn(i, j, k, cUVelField);
					const double rbneg = rbn(i, j, k, cUVelField);

					const double repos = rep(i, j, k, cUVelField);
					const double rwpos = rwp(i, j, k, true, cUVelField);
					const double rnpos = rnp(i, j, k, cUVelField);
					const double rspos = rsp(i, j, k, false, cUVelField);
					const double rtpos = rtp(i, j, k, cUVelField);
					const double rbpos = rbp(i, j, k, false, cUVelField);

					const double phiE = cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cUVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cUVelField[((j+1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cUVelField[((j-1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cUVelField[(j * cNZ * cNX) + ((k+1) * cNX) + i];
					const double phiB = cUVelField[(j * cNZ * cNX) + ((k-1) * cNX) + i];

					const double SuDc =
						(0.5f * Fe * (((1.0f - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5f * Fw * (-((1.0f - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5f * Fn * (((1.0f - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5f * Fs * (-((1.0f - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5f * Ft * (((1.0f - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5f * Fb * (-((1.0f - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));


					// SERIA UM ERRO SETAR ESSES COEFS A ZERO (Os valores da fronteira sao utilizados)
					//if (i == 2)
					//	aijkU[J][K][i].aimjk = 0.0f;
					//if (i == NXx - 3)
					//	aijkU[J][K][i].aipjk = 0.0f;
					// ^^ SERIA UM ERRO SETAR ESSES COEFS A ZERO (Os valores da fronteira sao utilizados)
					// FIM
					if (k == 1)
						aijkU.aijkm = 0.0f;
					if (k == cNZ - 2)
						aijkU.aijkp = 0.0f;
					if (j == 1)
						aijkU.aijmk = 0.0f;
					if (j == cNY - 2)
						aijkU.aijpk = 0.0f;

					aijkU.aijk = aijkU.aimjk + aijkU.aipjk + aijkU.aijmk
						+ aijkU.aijpk + aijkU.aijkm + aijkU.aijkp + deltaf - cSpUField[(j * cNZ * cNX) + (k * cNX) + i]
						+ (cRHO * cdV / dt);

					caijkU[(j * cNZ * cNX) + (k * cNX) + i] = aijkU.aijk;

					aijkU.source = SuDc + 
						((cPressureField[(j * cNZ * cNX) + (k * cNX) + i - 1] - cPressureField[(j * cNZ * cNX) + (k * cNX) + i]) * cdA) +
						cSuUField[(j * cNZ * cNX) + (k * cNX) + i] + 
						((cRHO * cdd * cdd * cdd / dt) * uField0[(j * cNZ * cNX) + (k * cNX) + i]);

					//if (isnan(aijkU.aijk) || isnan(aijkU.source))
					//{
					//	int x = 0;
					//}

					// Interior point. Use 5-point finite difference stencil.

					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jkmi
					col.push_back(index + (-1 * cNX) + (0 * cNX * cNZ));
					val.push_back(- aijkU.aijkm);

					// jkim
					col.push_back(index - 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkU.aimjk);

					// jki
					col.push_back(index + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(aijkU.aijk);

					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkU.aipjk);

					// jkpi
					col.push_back(index + (1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkU.aijkp);

					// jmki
					col.push_back(index + (0 * cNX) + (-1 * cNX * cNZ));
					val.push_back(-aijkU.aijmk);

					// jpki
					col.push_back(index + (0 * cNX) + (1 * cNX * cNZ));
					val.push_back(-aijkU.aijpk);

					rhs[index] = aijkU.source;
				}

				ptr.push_back(col.size());
			}
		}		
	}

	return n2;
}

int Field::CreateVMomentumLSCSR(const doubleField3D& vField0, double dt, std::vector<int>& ptr, 
	std::vector<int>& col, std::vector<double>& val, std::vector<double>& rhs)
{
	int n2 = cNX * cNY * cNZ; // Number of points in the grid.
	
	/*
	if (ptr.size() != n2 + 1)
	{
		ptr.clear();
		ptr.reserve(n2 + 1);
		ptr.push_back(0);
	}

	if (col.size() != n2 * 7)
	{
		col.clear();
		col.reserve(n2 * 7); // sete coeficientes
	}

	if (val.size() != n2 * 7)
	{
		val.clear();
		val.reserve(n2 * 7);  // sete coeficientes
	}

	if (rhs.size() != n2)
	{
		rhs.resize(n2);
	}
	*/
	ptr.clear();
	if(ptr.capacity() != n2 +1)
		ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	if (col.capacity() != n2 * 7)
		col.reserve(n2 * 7);

	val.clear();
	if (val.capacity() != n2 * 7)
		val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);

	for (int j = 0, index = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++, index++) {

				if ( i == 0 || i == cNX - 1 || j == 0 || j == 1 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

					col.push_back(index);
					val.push_back(1.0f);

					rhs[index] = cVVelField[(j * cNZ * cNX) + (k * cNX) + i];
				}
				else
				{
					CoefsData aijkV;

					const double Fw = GetFwInternalVMomentum(i, j, k, cRHO) * cdA;
					const double Fe = GetFeInternalVMomentum(i, j, k, cRHO) * cdA;
					const double Fs = GetFsInternalVMomentum(i, j, k, cRHO) * cdA;
					const double Fn = GetFnInternalVMomentum(i, j, k, cRHO) * cdA;
					const double Fb = GetFbInternalVMomentum(i, j, k, cRHO) * cdA;
					const double Ft = GetFtInternalVMomentum(i, j, k, cRHO) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkV.aimjk = cD + fmaxf(Fw, 0.0f);
					aijkV.aipjk = cD + fmaxf(-Fe, 0.0f);

					aijkV.aijmk = cD + fmaxf(Fs, 0.0f);
					aijkV.aijpk = cD + fmaxf(-Fn, 0.0f);

					aijkV.aijkm = cD + fmaxf(Fb, 0.0f);
					aijkV.aijkp = cD + fmaxf(-Ft, 0.0f);

					const double alfaw = (Fw > 0) ? 1.0f : 0.0f;
					const double alfae = (Fe > 0) ? 1.0f : 0.0f;
					const double alfas = (Fs > 0) ? 1.0f : 0.0f;
					const double alfan = (Fn > 0) ? 1.0f : 0.0f;
					const double alfab = (Fb > 0) ? 1.0f : 0.0f;
					const double alfat = (Ft > 0) ? 1.0f : 0.0f;

					const double reneg = ren(i, j, k, cVVelField);
					const double rwneg = rwn(i, j, k, cVVelField);
					const double rnneg = rnn(i, j, k, cVVelField);
					const double rsneg = rsn(i, j, k, cVVelField);
					const double rtneg = rtn(i, j, k, cVVelField);
					const double rbneg = rbn(i, j, k, cVVelField);

					const double repos = rep(i, j, k, cVVelField);
					const double rwpos = rwp(i, j, k, false, cVVelField);
					const double rnpos = rnp(i, j, k, cVVelField);
					const double rspos = rsp(i, j, k, true, cVVelField);
					const double rtpos = rtp(i, j, k, cVVelField);
					const double rbpos = rbp(i, j, k, false, cVVelField);

					const double phiE = cVVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cVVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cVVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cVVelField[((j+1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cVVelField[((j-1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cVVelField[(j * cNZ * cNX) + ((k+1) * cNX) + i];
					const double phiB = cVVelField[(j * cNZ * cNX) + ((k-1) * cNX) + i];


					const double SuDc = (0.5f * Fe * (((1.0f - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5f * Fw * (-((1.0f - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5f * Fn * (((1.0f - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5f * Fs * (-((1.0f - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5f * Ft * (((1.0f - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5f * Fb * (-((1.0f - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));

					//if (K == 1)
					//	(*pOutV)[J][K][I].aijkm = 0.0f;
					//if (K == cNZ - 2)
					//	(*pOutV)[J][K][I].aijkp = 0.0f;
					//if (I == 1)
					//	(*pOutV)[J][K][I].aimjk = 0.0f;
					//if (I == cNX - 2)
					//	(*pOutV)[J][K][I].aipjk = 0.0f;
					if (j == 2)
						aijkV.aijmk = 0.0f;
					if (j == cNY - 2)
						aijkV.aijpk = 0.0f;

					aijkV.aijk = aijkV.aimjk + aijkV.aipjk + aijkV.aijmk
						+ aijkV.aijpk + aijkV.aijkm + aijkV.aijkp + deltaf - cSpVField[(j * cNZ * cNX) + (k * cNX) + i]
						+ (cRHO * cdV / dt);
					
					caijkV[(j * cNZ * cNX) + (k * cNX) + i] = aijkV.aijk;

					aijkV.source = SuDc + 
						((cPressureField[((j-1) * cNZ *cNX) + (k*cNX) + i] - cPressureField[(j * cNZ * cNX) + (k * cNX) + i])*cdA) +
						cSuVField[(j * cNZ * cNX) + (k * cNX) + i] + 
						((cRHO * cdd * cdd * cdd / dt) * vField0[(j * cNZ * cNX) + (k * cNX) + i]);

					//if (isnan(aijkV.aijk) || isnan(aijkV.source))
					//{
					//	int x = 0;
					//}


					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jkmi
					col.push_back(index + (-1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkV.aijkm);

					// jkim
					col.push_back(index - 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkV.aimjk);

					// jki
					col.push_back(index + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(aijkV.aijk);

					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkV.aipjk);

					// jkpi
					col.push_back(index + (1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkV.aijkp);

					// jmki
					col.push_back(index + (0 * cNX) + (-1 * cNX * cNZ));
					val.push_back(-aijkV.aijmk);

					// jpki
					col.push_back(index + (0 * cNX) + (1 * cNX * cNZ));
					val.push_back(-aijkV.aijpk);

					rhs[index] = aijkV.source;
				}

				ptr.push_back(col.size());
			}
		}
	}

	return n2;
}
int Field::CreateWMomentumLSCSR(const doubleField3D& wField0, double dt, std::vector<int>& ptr, 
	std::vector<int>& col, std::vector<double>& val, std::vector<double>& rhs)
{
	int n2 = cNX * cNY * cNZ; // Number of points in the grid.
	/*

	if (ptr.size() != n2 + 1)
	{
		ptr.clear();
		ptr.reserve(n2 + 1);
		ptr.push_back(0);
	}

	if (col.size() != n2 * 7)
	{
		col.clear();
		col.reserve(n2 * 7); // sete coeficientes
	}

	if (val.size() != n2 * 7)
	{
		val.clear();
		val.reserve(n2 * 7);  // sete coeficientes
	}

	if (rhs.size() != n2)
	{
		rhs.resize(n2);
	}
	*/
	ptr.clear();
	if (ptr.capacity() != n2 + 1)
		ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	if (col.capacity() != n2 * 7)
		col.reserve(n2 * 7);

	val.clear();
	if (val.capacity() != n2 * 7)
		val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);

	for (int j = 0, index = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++, index++) {

				if (i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == 1|| k == cNZ - 1)
				{
					// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

					col.push_back(index);
					val.push_back(1.0f);

					rhs[index] = cWVelField[(j * cNZ * cNX) + (k * cNX) + i];
				}
				else
				{
					CoefsData aijkW;

					const double Fw = GetFwInternalWMomentum(i, j, k, cRHO) * cdA;
					const double Fe = GetFeInternalWMomentum(i, j, k, cRHO) * cdA;
					const double Fs = GetFsInternalWMomentum(i, j, k, cRHO) * cdA;
					const double Fn = GetFnInternalWMomentum(i, j, k, cRHO) * cdA;
					const double Fb = GetFbInternalWMomentum(i, j, k, cRHO) * cdA;
					const double Ft = GetFtInternalWMomentum(i, j, k, cRHO) * cdA;

					const double deltaf = Fe - Fw + Fn - Fs + Ft - Fb;

					aijkW.aimjk = cD + fmaxf(Fw, 0.0f);
					aijkW.aipjk = cD + fmaxf(-Fe, 0.0f);

					aijkW.aijmk = cD + fmaxf(Fs, 0.0f);
					aijkW.aijpk = cD + fmaxf(-Fn, 0.0f);

					aijkW.aijkm = cD + fmaxf(Fb, 0.0f);
					aijkW.aijkp = cD + fmaxf(-Ft, 0.0f);

					const double alfaw = (Fw > 0) ? 1.0f : 0.0f;
					const double alfae = (Fe > 0) ? 1.0f : 0.0f;
					const double alfas = (Fs > 0) ? 1.0f : 0.0f;
					const double alfan = (Fn > 0) ? 1.0f : 0.0f;
					const double alfab = (Fb > 0) ? 1.0f : 0.0f;
					const double alfat = (Ft > 0) ? 1.0f : 0.0f;

					const double reneg = ren(i, j, k, cWVelField);
					const double rwneg = rwn(i, j, k, cWVelField);
					const double rnneg = rnn(i, j, k, cWVelField);
					const double rsneg = rsn(i, j, k, cWVelField);
					const double rtneg = rtn(i, j, k, cWVelField);
					const double rbneg = rbn(i, j, k, cWVelField);

					const double repos = rep(i, j, k, cWVelField);
					const double rwpos = rwp(i, j, k, false, cWVelField);
					const double rnpos = rnp(i, j, k, cWVelField);
					const double rspos = rsp(i, j, k, false, cWVelField);
					const double rtpos = rtp(i, j, k, cWVelField);
					const double rbpos = rbp(i, j, k, true, cWVelField);

					const double phiE = cWVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
					const double phiP = cWVelField[(j * cNZ * cNX) + (k * cNX) + i];
					const double phiW = cWVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
					const double phiN = cWVelField[((j+1) * cNZ * cNX) + (k * cNX) + i];
					const double phiS = cWVelField[((j-1) * cNZ * cNX) + (k * cNX) + i];
					const double phiT = cWVelField[(j * cNZ * cNX) + ((k+1) * cNX) + i];
					const double phiB = cWVelField[(j * cNZ * cNX) + ((k-1) * cNX) + i];

					const double SuDc = (0.5f * Fe * (((1.0f - alfae) * Psir(reneg)) - (alfae * Psir(repos))) * (phiE - phiP))
						+ (0.5f * Fw * (-((1.0f - alfaw) * Psir(rwneg)) + (alfaw * Psir(rwpos))) * (phiP - phiW))
						+ (0.5f * Fn * (((1.0f - alfan) * Psir(rnneg)) - (alfan * Psir(rnpos))) * (phiN - phiP))
						+ (0.5f * Fs * (-((1.0f - alfas) * Psir(rsneg)) + (alfas * Psir(rspos))) * (phiP - phiS))
						+ (0.5f * Ft * (((1.0f - alfat) * Psir(rtneg)) - (alfat * Psir(rtpos))) * (phiT - phiP))
						+ (0.5f * Fb * (-((1.0f - alfab) * Psir(rbneg)) + (alfab * Psir(rbpos))) * (phiP - phiB));

					if (k == 2)
						aijkW.aijkm = 0.0f;
					if (k == cNZ - 2)
						aijkW.aijkp = 0.0f;
					//if (J == 1)
					//	(*pOutW)[J][K][I].aijmk = 0.0f;
					//if (J == cNY - 2)
					//	(*pOutW)[J][K][I].aijpk = 0.0f;
					//if (I == 1)
					//	(*pOutW)[J][K][I].aimjk = 0.0f;
					//if (I == cNX - 2)
					//	(*pOutW)[J][K][I].aipjk = 0.0f;

					aijkW.aijk = aijkW.aimjk + aijkW.aipjk + aijkW.aijmk
						+ aijkW.aijpk + aijkW.aijkm + aijkW.aijkp + deltaf - cSpWField[(j * cNZ * cNX) + (k * cNX) + i]
						+ (cRHO * cdV / dt);

					caijkW[(j * cNZ * cNX) + (k * cNX) + i] = aijkW.aijk;

					aijkW.source = SuDc + ((cPressureField[(j * cNZ * cNX) + ((k-1) * cNX) + i] - cPressureField[(j * cNZ * cNX) + (k * cNX) + i]) * cdA) +
						cSuWField[(j * cNZ * cNX) + (k * cNX) + i] + ((cRHO * cdd * cdd * cdd / dt) * wField0[(j * cNZ * cNX) + (k * cNX) + i]);



					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jkmi
					col.push_back(index + (-1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkW.aijkm);

					// jkim
					col.push_back(index - 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkW.aimjk);

					// jki
					col.push_back(index + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(aijkW.aijk);

					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkW.aipjk);

					// jkpi
					col.push_back(index + (1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkW.aijkp);

					// jmki
					col.push_back(index + (0 * cNX) + (-1 * cNX * cNZ));
					val.push_back(-aijkW.aijmk);

					// jpki
					col.push_back(index + (0 * cNX) + (1 * cNX * cNZ));
					val.push_back(-aijkW.aijpk);

					rhs[index] = aijkW.source;
				}

				ptr.push_back(col.size());
			}
		}
	}

	return n2;
}

int Field::CreatePressureCorrectionLSCSR(std::vector<int>& ptr, std::vector<int>& col, std::vector<double>& val,
	std::vector<double>& rhs) const
{
	int n2 = cNX * cNY * cNZ; // Number of points in the grid.
	/*
	if (ptr.size() != n2 + 1)
	{
		ptr.clear();
		ptr.reserve(n2 + 1);
		ptr.push_back(0);
	}

	if (col.size() != n2 * 7)
	{
		col.clear();
		col.reserve(n2 * 7); // sete coeficientes
	}

	if (val.size() != n2 * 7)
	{
		val.clear();
		val.reserve(n2 * 7);  // sete coeficientes
	}

	if (rhs.size() != n2)
	{
		rhs.resize(n2);
	}
	*/

	ptr.clear();
	ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	col.reserve(n2 * 7);

	val.clear();
	val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);
	for (int j = 0, index = 0; j < cNY; ++j)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++, index++) {

				if (i == 0 || i == cNX - 1 || j == 0 || j == cNY - 1 || k == 0 || k == cNZ - 1)
				{
					// Boundary point. Use Dirichlet condition. (seta o valor da variavel diretamente)

					col.push_back(index);
					val.push_back(1.0f);

					rhs[index] = cPressureCorrField[(j * cNZ * cNX) + (k * cNX) + i];
				}
				else if (i == cNX / 2 && j == cNY / 2 && k == cNZ / 2)
				{
					col.push_back(index);
					val.push_back(1.0f);

					rhs[index] = 0.0f;
				}
				else
				{
					CoefsData aijkPCorr;

					aijkPCorr.aimjk = cRHO * cdA * cdA / caijkU[(j * cNZ * cNX) + (k * cNX) + i];
					aijkPCorr.aipjk = cRHO * cdA * cdA / caijkU[(j * cNZ * cNX) + (k * cNX) + i + 1];

					aijkPCorr.aijmk = cRHO * cdA * cdA / caijkV[(j * cNZ * cNX) + (k * cNX) + i];
					aijkPCorr.aijpk = cRHO * cdA * cdA / caijkV[((j+1) * cNZ * cNX) + (k * cNX) + i];

					aijkPCorr.aijkm = cRHO * cdA * cdA / caijkW[(j * cNZ * cNX) + (k * cNX) + i];
					aijkPCorr.aijkp = cRHO * cdA * cdA / caijkW[(j * cNZ * cNX) + ((k+1) * cNX) + i];

					if (i == 1)
						aijkPCorr.aimjk = 0.0f;
					else if (i == cNX - 2)
						aijkPCorr.aipjk = 0.0f;
					if (j == 1)
						aijkPCorr.aijmk = 0.0f;
					else if (j == cNY - 2)
						aijkPCorr.aijpk = 0.0f;
					if (k == 1)
						aijkPCorr.aijkm = 0.0f;
					else if (k == cNZ - 2)
						aijkPCorr.aijkp = 0.0f;


					aijkPCorr.aijk = (aijkPCorr.aimjk + aijkPCorr.aipjk +
						aijkPCorr.aijmk + aijkPCorr.aijpk + aijkPCorr.aijkm +
						aijkPCorr.aijkp);

					aijkPCorr.source = cRHO * cdA * (cUVelField[(j * cNZ * cNX) + (k * cNX) + i] - cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1] +
						cVVelField[(j * cNZ * cNX) + (k * cNX) + i] - cVVelField[((j+1) * cNZ * cNX) + (k * cNX) + i] + 
						cWVelField[(j * cNZ * cNX) + (k * cNX) + i] - cWVelField[(j * cNZ * cNX) + ((k+1) * cNX) + i]);

					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jkmi
					col.push_back(index + (-1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkPCorr.aijkm);

					// jkim
					col.push_back(index - 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkPCorr.aimjk);

					// jki
					col.push_back(index + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(aijkPCorr.aijk);

					// jkip
					col.push_back(index + 1 + (0 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkPCorr.aipjk);

					// jkpi
					col.push_back(index + (1 * cNX) + (0 * cNX * cNZ));
					val.push_back(-aijkPCorr.aijkp);

					// jmki
					col.push_back(index + (0 * cNX) + (-1 * cNX * cNZ));
					val.push_back(-aijkPCorr.aijmk);

					// jpki
					col.push_back(index + (0 * cNX) + (1 * cNX * cNZ));
					val.push_back(-aijkPCorr.aijpk);

					rhs[index] = aijkPCorr.source;
				}

				ptr.push_back(col.size());
			}
		}
	}

	return n2;
}

// Set P value for all 'I' and 'J' in K defined;
void Field::SetPValueIJ(size_t K, double pValue)
{

	for (size_t J = 0; J < cNY; J++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			cPressureField[(J * cNZ * cNX) + (K * cNX) + I] = pValue;
		}
	}
}

// Set P value for all 'J' and 'K' in I defined;
void Field::SetPValueJK(size_t I, double pValue)
{
	for (size_t J = 0; J < cNY; J++)
	{
		for (size_t K = 0; K < cNZ; K++)
		{
			cPressureField[(J * cNZ * cNX) + (K * cNX) + I] = pValue;
		}
	}
}

// Set P value for all 'K' and 'I' in J defined; 
void Field::SetPValueIK(size_t J, double pValue)
{
	for (size_t K = 0; K < cNZ; K++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			cPressureField[(J * cNZ * cNX) + (K * cNX) + I] = pValue;
		}
	}
}

// Set U value for all 'i' and 'J' in K defined;
void Field::SetUValueiJ(size_t K, double uValue)
{
	for (size_t J = 0; J < cNY; J++)
	{
		for (size_t i = 0; i < cNX; i++)
		{
			cUVelField[(J * cNZ * cNX) + (K * cNX) + i] = uValue;
		}
	}
}

// Set U value for all 'J' and 'K' in i defined;
void Field::SetUValueJK(size_t i, double uValue)
{
	for (size_t J = 0; J < cNY; J++)
	{
		for (size_t K = 0; K < cNZ; K++)
		{
			cUVelField[(J * cNZ * cNX) + (K * cNX) + i] = uValue;
		}
	}
}

// Set U value for all 'K' and 'i' in J defined; 
void Field::SetUValueiK(size_t J, double uValue)
{
	for (size_t K = 0; K < cNZ; K++)
	{
		for (size_t i = 0; i < cNX; i++)
		{
			cUVelField[(J * cNZ * cNX) + (K * cNX) + i] = uValue;
		}
	}
}

// Set V value for all 'I' and 'K' in j defined;
void Field::SetVValueIK(size_t j, double uValue)
{
	for (size_t K = 0; K < cNZ; K++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			cVVelField[(j * cNZ * cNX) + (K * cNX) + I] = uValue;
		}
	}
}

// Set V value for all 'j' and 'I' in K defined;
void Field::SetVValueIj(size_t K, double vValue)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			cVVelField[(j * cNZ * cNX) + (K * cNX) + I] = vValue;
		}
	}
}

// Set V value for all 'K' and 'j' in I defined; 
void Field::SetVValuejK(size_t I, double vValue)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t K = 0; K < cNZ; K++)
		{
			cVVelField[(j * cNZ * cNX) + (K * cNX) + I] = vValue;
		}
	}
}

// Set W value for all 'I' and 'J' in k defined;
void Field::SetWValueIJ(size_t k, double wValue)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = wValue;
		}
	}
}

// Set W value for all 'J' and 'k' in I defined;
void Field::SetWValueJk(size_t I, double wValue)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = wValue;
		}
	}
}

// Set W value for all 'k' and 'I' in J defined; 
void Field::SetWValueIk(size_t J, double wValue)
{
	for (size_t k = 0; k < cNZ; k++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			cWVelField[(J * cNZ * cNX) + (k * cNX) + I] = wValue;
		}
	}
}

void Field::SetPValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& KK, double pValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cPressureField[(j * cNZ * cNX) + (k * cNX) + i] = pValue;
			}
		}
	}
}

void Field::SetUValueInterval(const Intervalo& ii, const Intervalo& JJ, const Intervalo& KK, double uValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			for (size_t i = ii.valor1; i < ii.valor2; i++)
			{
				cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = uValue;
			}
		}
	}
}

void Field::SetVValueInterval(const Intervalo& II, const Intervalo& jj, const Intervalo& KK, double vValue)
{
	for (size_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cVVelField[(j * cNZ * cNX) + (k * cNX) + i] = vValue;
			}
		}
	}
}

void Field::SetWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double wValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cWVelField[(j * cNZ * cNX) + (k * cNX) + i] = wValue;
			}
		}
	}
}

void Field::ExtrapolateForwardPJK(size_t I)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cPressureField[(j * cNZ * cNX) + (k * cNX) + I] = cPressureField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void Field::ExtrapolateForwardUJK(size_t i)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
		}
	}
}

void Field::ExtrapolateForwardVjK(size_t I)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void Field::ExtrapolateForwardWJk(size_t I)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void Field::ExtrapolateBackwardPJK(size_t I)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cPressureField[(j * cNZ * cNX) + (k * cNX) + I] = cPressureField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void Field::ExtrapolateBackwardUJK(size_t i)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
		}
	}
}

void Field::ExtrapolateBackwardVjK(size_t I)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void Field::ExtrapolateBackwardWJk(size_t I)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void Field::ExtrapolateBackwardWeightedUJK(size_t i, double weight)
{
	for (size_t j = 1; j < cNY - 1; j++)
	{
		for (size_t k = 1; k < cNZ - 1; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1] * weight;
		}
	}
}

void Field::ExtrapolateBackwardWeightedVjK(size_t I, double weight)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I - 1] * weight;
		}
	}
}

void Field::ExtrapolateBackwardWeightedWJk(size_t I, double weight)
{
	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I - 1] * weight;
		}
	}
}

void Field::ExtrapolateForwardIntervalPJK(size_t I, const Intervalo& JJ, const Intervalo& KK)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			cPressureField[(j * cNZ * cNX) + (k * cNX) + I] = cPressureField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void Field::ExtrapolateForwardIntervalUJK(size_t i, const Intervalo& JJ, const Intervalo& KK)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i + 1];
		}
	}
}

void Field::ExtrapolateForwardIntervalVjK(size_t I, const Intervalo& jj, const Intervalo& KK)
{
	for (size_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void Field::ExtrapolateForwardIntervalWJk(size_t I, const Intervalo& JJ, const Intervalo& kk)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I + 1];
		}
	}
}

void Field::ExtrapolateBackwardIntervalPJK(size_t I, const Intervalo& JJ, const Intervalo& KK)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			cPressureField[(j * cNZ * cNX) + (k * cNX) + I] = cPressureField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void Field::ExtrapolateBackwardIntervalUJK(size_t i, const Intervalo& JJ, const Intervalo& KK)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			cUVelField[(j * cNZ * cNX) + (k * cNX) + i] = cUVelField[(j * cNZ * cNX) + (k * cNX) + i - 1];
		}
	}
}

void Field::ExtrapolateBackwardIntervalVjK(size_t I, const Intervalo& jj, const Intervalo& KK)
{
	for (size_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (size_t k = KK.valor1; k < KK.valor2; k++)
		{
			cVVelField[(j * cNZ * cNX) + (k * cNX) + I] = cVVelField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void Field::ExtrapolateBackwardIntervalWJk(size_t I, const Intervalo& JJ, const Intervalo& kk)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			cWVelField[(j * cNZ * cNX) + (k * cNX) + I] = cWVelField[(j * cNZ * cNX) + (k * cNX) + I - 1];
		}
	}
}

void Field::ExtrapolateForwardIntervalPJI(size_t K, const Intervalo& II, const Intervalo& JJ)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t i = II.valor1; i < II.valor2; i++)
		{
			cPressureField[(j * cNZ * cNX) + (K * cNX) + i] = cPressureField[(j * cNZ * cNX) + ((K+1) * cNX) + i];
		}
	}
}

void Field::ExtrapolateForwardIntervalUJI(size_t K, const Intervalo& ii, const Intervalo& JJ)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t i = ii.valor1; i < ii.valor2; i++)
		{
			cUVelField[(j * cNZ * cNX) + ((K ) * cNX) + i] = cUVelField[(j * cNZ * cNX) + ((K + 1) * cNX) + i];
		}
	}
}
void Field::ExtrapolateForwardIntervalVJI(size_t K, const Intervalo& II, const Intervalo& jj)
{
	for (size_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (size_t i = II.valor1; i < II.valor2; i++)
		{
			cVVelField[(j * cNZ * cNX) + ((K) * cNX) + i] = cVVelField[(j * cNZ * cNX) + ((K + 1) * cNX) + i];
		}
	}
}
void Field::ExtrapolateForwardIntervalWJI(size_t k, const Intervalo& II, const Intervalo& JJ)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t i = II.valor1; i < II.valor2; i++)
		{
			cVVelField[(j * cNZ * cNX) + ((k) * cNX) + i] = cVVelField[(j * cNZ * cNX) + ((k + 1) * cNX) + i];
		}
	}
}

void Field::ExtrapolateBackwardIntervalPJI(size_t K, const Intervalo& II, const Intervalo& JJ)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t i = II.valor1; i < II.valor2; i++)
		{
			cPressureField[(j * cNZ * cNX) + ((K) * cNX) + i] = cPressureField[(j * cNZ * cNX) + ((K - 1) * cNX) + i];
		}
	}
}
void Field::ExtrapolateBackwardIntervalUJI(size_t K, const Intervalo& ii, const Intervalo& JJ)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t i = ii.valor1; i < ii.valor2; i++)
		{
			cUVelField[(j * cNZ * cNX) + ((K + 1) * cNX) + i] = cUVelField[(j * cNZ * cNX) + ((K - 1) * cNX) + i];
		}
	}
}
void Field::ExtrapolateBackwardIntervalVJI(size_t K, const Intervalo& II, const Intervalo& jj)
{
	for (size_t j = jj.valor1; j < jj.valor2; j++)
	{
		for (size_t i = II.valor1; i < II.valor2; i++)
		{
			cVVelField[(j * cNZ * cNX) + ((K) * cNX) + i] = cVVelField[(j * cNZ * cNX) + ((K - 1) * cNX) + i];
		}
	}
}
void Field::ExtrapolateBackwardIntervalWJI(size_t k, const Intervalo& II, const Intervalo& JJ)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t i = II.valor1; i < II.valor2; i++)
		{
			cVVelField[(j * cNZ * cNX) + ((k) * cNX) + i] = cVVelField[(j * cNZ * cNX) + ((k - 1) * cNX) + i];
		}
	}
}

void Field::SaveIJCutFieldToCSVFile(size_t K, const std::string& baseFileName) const
{
	std::ofstream pFile(baseFileName + "P.csv");

	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			pFile << cPressureField[(j * cNZ * cNX) + ((K) * cNX) + I] << ";";
		}
		pFile << "\n";
	}

	pFile.close();

	std::ofstream uFile(baseFileName + "U.csv");

	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			uFile << cUVelField[(j * cNZ * cNX) + ((K)*cNX) + I] << ";";
		}
		uFile << "\n";
	}

	uFile.close();

	std::ofstream vFile(baseFileName + "V.csv");

	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			vFile << cVVelField[(j * cNZ * cNX) + ((K)*cNX) + I] << ";";
		}
		vFile << "\n";
	}

	vFile.close();

	std::ofstream wFile(baseFileName + "W.csv");

	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			wFile << cWVelField[(j * cNZ * cNX) + ((K)*cNX) + I] << ";";
		}
		wFile << "\n";
	}

	wFile.close();

	std::ofstream uspFile(baseFileName + "USP.csv");

	for (size_t J = 0; J < cNY; J++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			uspFile << cSpUField[(J * cNZ * cNX) + ((K)*cNX) + I] << ";";
		}
		uspFile << "\n";
	}

	uspFile.close();
}

void Field::SaveIKCutFieldToCSVFile(size_t J, const std::string& baseFileName) const
{
	std::ofstream pFile(baseFileName + "P.csv");

	for (size_t K = 0; K < cNZ; K++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			pFile << cPressureField[(J * cNZ * cNX) + ((K)*cNX) + I] << ";";
		}
		pFile << "\n";
	}

	pFile.close();

	//std::ofstream pcorrFile(baseFileName + "Pcorr.csv");
	//
	//for (size_t K = 0; K < cPressureCorrField[0].size(); K++)
	//{
	//	for (size_t I = 0; I < cPressureCorrField[0][0].size(); I++)
	//	{
	//		pcorrFile << cPressureCorrField[J][K][I] << ";";
	//	}
	//	pcorrFile << "\n";
	//}
	//
	//pcorrFile.close();

	std::ofstream uFile(baseFileName + "U.csv");

	for (size_t K = 0; K < cNZ; K++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			uFile << cUVelField[(J* cNZ * cNX) + ((K)*cNX) + I] << ";";
		}
		uFile << "\n";
	}

	uFile.close();

	std::ofstream vFile(baseFileName + "V.csv");

	for (size_t K = 0; K < cNZ; K++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			vFile << cVVelField[(J * cNZ * cNX) + ((K)*cNX) + I] << ";";
		}
		vFile << "\n";
	}

	vFile.close();

	std::ofstream wFile(baseFileName + "W.csv");

	for (size_t K = 0; K < cNZ; K++)
	{
		for (size_t I = 0; I < cNX; I++)
		{
			wFile << cWVelField[(J * cNZ * cNX) + ((K)*cNX) + I] << ";";
		}
		wFile << "\n";
	}

	wFile.close();

	//std::ofstream uspFile(baseFileName + "USP.csv");
	//
	//for (size_t K = 0; K < cSpUField[0].size(); K++)
	//{
	//	for (size_t I = 0; I < cSpUField[0][0].size(); I++)
	//	{
	//		uspFile << cSpUField[J][K][I] << ";";
	//	}
	//	uspFile << "\n";
	//}
	//
	//uspFile.close();
}

void Field::SaveJKCutFieldToCSVFile(size_t I, const std::string& baseFileName) const
{
	std::ofstream pFile(baseFileName + "P.csv");

	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			pFile << cPressureField[(j * cNZ * cNX) + ((k)*cNX) + I] << ";";
		}
		pFile << "\n";
	}
	pFile.close();

	std::ofstream uFile(baseFileName + "U.csv");

	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			uFile << cUVelField[(j * cNZ * cNX) + ((k)*cNX) + I] << ";";
		}
		uFile << "\n";
	}
	uFile.close();

	std::ofstream vFile(baseFileName + "V.csv");

	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			vFile << cVVelField[(j * cNZ * cNX) + ((k) * cNX) + I] << ";";
		}
		vFile << "\n";
	}
	vFile.close();

	std::ofstream wFile(baseFileName + "W.csv");

	for (size_t j = 0; j < cNY; j++)
	{
		for (size_t k = 0; k < cNZ; k++)
		{
			wFile << cWVelField[(j * cNZ * cNX) + ((k)*cNX) + I] << ";";
		}
		wFile << "\n";
	}
	wFile.close();

	//std::ofstream uspFile(baseFileName + "USP.csv");
	//for (size_t j = 0; j < cSpUField.size(); j++)
	//{
	//	for (size_t K = 0; K < cSpUField[0].size(); K++)
	//	{
	//		uspFile << cSpUField[j][K][I] << ";";
	//	}
	//
	//	uspFile << "\n";
	//}
	//
	//
	//uspFile.close();
}

double Field::GetFwInternalUMomentum(size_t i, size_t J, size_t K, double rho) const
{
	return (rho / 2.0f) * (cUVelField[(J * cNZ * cNX) + ((K)*cNX) + i] + cUVelField[(J * cNZ * cNX) + ((K)*cNX) + i - 1]);
}
double Field::GetFeInternalUMomentum(size_t i, size_t J, size_t K, double rho) const
{
	return (rho / 2.0f) * (cUVelField[(J * cNZ * cNX) + ((K)*cNX) + i + 1] + cUVelField[(J * cNZ * cNX) + ((K)*cNX) + i]);
}
double Field::GetFsInternalUMomentum(size_t i, size_t J, size_t K, double rho) const
{
	return (rho / 2.0f) * (cVVelField[(J * cNZ * cNX) + ((K)*cNX) + i] + cVVelField[(J * cNZ * cNX) + ((K)*cNX) + i - 1]);
}
double Field::GetFnInternalUMomentum(size_t i, size_t J, size_t K, double rho) const
{
	return (rho / 2.0f) * (cVVelField[((J+1) * cNZ * cNX) + ((K)*cNX) + i] + cVVelField[((J+1) * cNZ * cNX) + ((K)*cNX) + i - 1]);
}
double Field::GetFbInternalUMomentum(size_t i, size_t J, size_t K, double rho) const
{
	return (rho / 2.0f) * (cWVelField[(J * cNZ * cNX) + ((K)*cNX) + i] + cWVelField[(J * cNZ * cNX) + ((K)*cNX) + i - 1]);
}
double Field::GetFtInternalUMomentum(size_t i, size_t J, size_t K, double rho) const
{
	return (rho / 2.0f) * (cWVelField[(J * cNZ * cNX) + ((K+1)*cNX) + i] + cWVelField[(J * cNZ * cNX) + ((K + 1)*cNX) + i - 1]);
}

double Field::GetFwInternalVMomentum(size_t I, size_t j, size_t K, double rho) const
{
	return (rho / 2.0f) * (cUVelField[(j * cNZ * cNX) + ((K)*cNX) + I] + cUVelField[((j-1) * cNZ * cNX) + ((K)*cNX) + I]);
}
double Field::GetFeInternalVMomentum(size_t I, size_t j, size_t K, double rho) const
{
	return (rho / 2.0f) * (cUVelField[(j * cNZ * cNX) + ((K)*cNX) + I + 1] + cUVelField[((j-1) * cNZ * cNX) + ((K)*cNX) + I + 1]);
}
double Field::GetFsInternalVMomentum(size_t I, size_t j, size_t K, double rho) const
{
	return (rho / 2.0f) * (cVVelField[(j * cNZ * cNX) + ((K)*cNX) + I] + cVVelField[((j-1) * cNZ * cNX) + ((K)*cNX) + I]);
}
double Field::GetFnInternalVMomentum(size_t I, size_t j, size_t K, double rho) const
{
	return (rho / 2.0f) * (cVVelField[(j * cNZ * cNX) + ((K)*cNX) + I] + cVVelField[((j+1) * cNZ * cNX) + ((K)*cNX) + I]);
}
double Field::GetFbInternalVMomentum(size_t I, size_t j, size_t K, double rho) const
{
	return (rho / 2.0f) * (cWVelField[(j * cNZ * cNX) + ((K)*cNX) + I] + cWVelField[((j-1) * cNZ * cNX) + ((K)*cNX) + I]);
}
double Field::GetFtInternalVMomentum(size_t I, size_t j, size_t K, double rho) const
{
	return (rho / 2.0f) * (cWVelField[(j * cNZ * cNX) + ((K+1)*cNX) + I] + cWVelField[((j-1) * cNZ * cNX) + ((K+1)*cNX) + I]);
}

double Field::GetFwInternalWMomentum(size_t I, size_t J, size_t k, double rho) const
{
	return (rho / 2.0f) * (cUVelField[(J * cNZ * cNX) + ((k)*cNX) + I] + cUVelField[(J * cNZ * cNX) + ((k-1)*cNX) + I]);
}
double Field::GetFeInternalWMomentum(size_t I, size_t J, size_t k, double rho) const
{
	return (rho / 2.0f) * (cUVelField[(J * cNZ * cNX) + ((k)*cNX) + I + 1] + cUVelField[(J * cNZ * cNX) + ((k-1)*cNX) + I + 1]);
}
double Field::GetFsInternalWMomentum(size_t I, size_t J, size_t k, double rho) const
{
	return (rho / 2.0f) * (cVVelField[(J * cNZ * cNX) + ((k)*cNX) + I] + cVVelField[(J * cNZ * cNX) + ((k-1)*cNX) + I]);
}
double Field::GetFnInternalWMomentum(size_t I, size_t J, size_t k, double rho) const
{
	return (rho / 2.0f) * (cVVelField[((J+1) * cNZ * cNX) + ((k)*cNX) + I] + cVVelField[((J+1) * cNZ * cNX) + ((k-1)*cNX) + I]);
}
double Field::GetFbInternalWMomentum(size_t I, size_t J, size_t k, double rho) const
{
	return (rho / 2.0f) * (cWVelField[(J * cNZ * cNX) + ((k)*cNX) + I] + cWVelField[(J * cNZ * cNX) + ((k - 1)*cNX) + I]);
}
double Field::GetFtInternalWMomentum(size_t I, size_t J, size_t k, double rho) const
{
	return (rho / 2.0f) * (cWVelField[(J * cNZ * cNX) + ((k)*cNX) + I] + cWVelField[(J * cNZ * cNX) + ((k+1)*cNX) + I]);
}

void Field::SetSpUValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cSpUField[(j * cNZ * cNX) + ((k)*cNX) + i] = SpValue;
			}
		}
	}
}
void Field::SetSuUValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cSuUField[(j * cNZ * cNX) + ((k)*cNX) + i] = SuValue;
			}
		}
	}
}

void Field::SetSpVValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cSpVField[(j * cNZ * cNX) + ((k)*cNX) + i] = SpValue;
			}
		}
	}
}
void Field::SetSuVValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cSuVField[(j * cNZ * cNX) + ((k)*cNX) + i] = SuValue;
			}
		}
	}
}

void Field::SetSpWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cSpWField[(j * cNZ * cNX) + ((k)*cNX) + i] = SpValue;
			}
		}
	}
}
void Field::SetSuWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue)
{
	for (size_t j = JJ.valor1; j < JJ.valor2; j++)
	{
		for (size_t k = kk.valor1; k < kk.valor2; k++)
		{
			for (size_t i = II.valor1; i < II.valor2; i++)
			{
				cSuWField[(j * cNZ * cNX) + ((k)*cNX) + i] = SuValue;
			}
		}
	}
}

double Field::Psir(double r) const
{
	return (r + (r * r)) / (1.0f + (r * r));
}

double Field::rep(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiW = field[(J * cNZ * cNX) + ((K)*cNX) + I - 1];
	const double phiE = field[(J * cNZ * cNX) + ((K)*cNX) + I + 1];

	if (phiE == phiP)
		return 1e15f;
	
	return (phiP - phiW) / (phiE - phiP);
}
double Field::ren(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiE = field[(J * cNZ * cNX) + ((K)*cNX) + I + 1];
	double phiEE;

	if (I == cNX - 2)
		phiEE = (2.0f * phiE) - phiP;
	else
		phiEE = field[(J * cNZ * cNX) + ((K)*cNX) + I + 2];

	if (phiE == phiP)
		return 0.0f;

	return (phiEE - phiE) / (phiE - phiP);
}
	 
double Field::rwp(size_t I, size_t J, size_t K, bool ufield, const doubleField3D& field) const
{
	const double phiW = field[(J * cNZ * cNX) + ((K)*cNX) + I - 1];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	const int iInicial = (ufield) ? 2 : 1;

	double phiWW;
	if (I == iInicial)
		phiWW = (2.0f * phiW) - phiP;
	else
		phiWW = field[(J * cNZ * cNX) + ((K)*cNX) + I - 2];

	if (phiP == phiW)
		return 0.0f;

	return (phiW - phiWW) / (phiP - phiW);
}
double Field::rwn(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiE = field[(J * cNZ * cNX) + ((K)*cNX) + I + 1];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiW = field[(J * cNZ * cNX) + ((K)*cNX) + I - 1];

	if (phiP == phiW)
		return 0.0f;

	return (phiE - phiP) / (phiP - phiW);
}
	 
double Field::rnp(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiS = field[((J-1) * cNZ * cNX) + ((K)*cNX) + I];
	const double phiN = field[((J+1) * cNZ * cNX) + ((K)*cNX) + I];

	if (phiN == phiP)
		return 0.0f;

	return (phiP - phiS) / (phiN - phiP);
}
double Field::rnn(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiN = field[((J+1) * cNZ * cNX) + ((K)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	double phiNN;
	if (J == cNY - 2)
		phiNN = (2.0f * phiN) - phiP;
	else
		phiNN = field[((J+2) * cNZ * cNX) + ((K)*cNX) + I];

	if (phiN == phiP)
		return 0.0f;

	return (phiNN - phiN) / (phiN - phiP);
}

double Field::rsp(size_t I, size_t J, size_t K, bool vfield, const doubleField3D& field) const
{
	const double phiS = field[((J-1) * cNZ * cNX) + ((K)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	const int jInicial = (vfield) ? 2 : 1;

	double phiSS;
	if (J == jInicial)
		phiSS = (phiS * 2.0f) - phiP;
	else
		phiSS = field[((J-2) * cNZ * cNX) + ((K)*cNX) + I];

	if (phiP == phiS)
		return 0.0f;

	return (phiS - phiSS) / (phiP - phiS);
}
double Field::rsn(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiN = field[((J+1) * cNZ * cNX) + ((K)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiS = field[((J-1) * cNZ * cNX) + ((K)*cNX) + I];

	if (phiP == phiS)
		return 0.0f;

	return (phiN - phiP) / (phiP - phiS);
}
	
double Field::rbp(size_t I, size_t J, size_t K, bool wfield, const doubleField3D& field) const
{
	const double phiB = field[(J * cNZ * cNX) + ((K-1)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	const int kInicial = wfield ?  2 : 1;

	double phiBB;
	if (K == kInicial)
		phiBB = (2.0f * phiB) - phiP;
	else
		phiBB = field[(J * cNZ * cNX) + ((K-2)*cNX) + I];

	if (phiP == phiB)
		return 0.0f;

	return (phiB - phiBB) / (phiP - phiB);
}
double Field::rbn(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiT = field[(J * cNZ * cNX) + ((K + 1)*cNX) + I];
	const double phiB = field[(J * cNZ * cNX) + ((K-1)*cNX) + I];

	if (phiP == phiB)
		return 0.0f;
	
	return (phiT - phiP) / (phiP - phiB);
}

double Field::rtp(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];
	const double phiT = field[(J * cNZ * cNX) + ((K+1)*cNX) + I];
	const double phiB = field[(J * cNZ * cNX) + ((K-1)*cNX) + I];

	if (phiT == phiP)
		return 0.0f;
	
	return (phiP - phiB) / (phiT - phiP);
}
double Field::rtn(size_t I, size_t J, size_t K, const doubleField3D& field) const
{
	const double phiT = field[(J * cNZ * cNX) + ((K+1)*cNX) + I];
	const double phiP = field[(J * cNZ * cNX) + ((K)*cNX) + I];

	double phiTT;
	if (K == cNZ - 2)
		phiTT = (2.0f * phiT) - phiP;
	else
		phiTT = field[(J * cNZ * cNX) + ((K+2)*cNX) + I];

	if (phiT == phiP)
		return 0.0f;

	return (phiTT - phiT) / (phiT - phiP);
}