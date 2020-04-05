#pragma once
#include <vector>
#include <string>

typedef std::vector<double> doubleField3D;

struct CoefsData
{
public:
	double aijk;
	double aipjk;
	double aimjk;
	double aijpk;
	double aijmk;
	double aijkp;
	double aijkm;

	double source;
};

class Field
{
public:
	class Intervalo
	{
	public:
		Intervalo(size_t v1, size_t v2);
		size_t valor1;
		size_t valor2;

		size_t Tamanho() const;
	};
public:
	Field(double L, double H, double W, double deltaSize, double rho, double mu, double Ufarfield, double Vfarfield, double Wfarfield);
	~Field();

	inline size_t GetNX() { return cNX; }
	inline size_t GetNY() { return cNY; }
	inline size_t GetNZ() { return cNZ; }
	
	// SIMPLE FLAT PLATE AT ZERO DEEGRE TO THE FLOW
	void SetBCFlatPlate();
	
	int CreateUMomentumLSCSR(
		const doubleField3D& uvel0, 
		double dt, 
		std::vector<int>& ptr, 
		std::vector<int>& col, 
		std::vector<double>& val, 
		std::vector<double>& rhs
	);
	int CreateVMomentumLSCSR(
		const doubleField3D& vvel0,
		double dt,
		std::vector<int>& ptr, 
		std::vector<int>& col, 
		std::vector<double>& val, 
		std::vector<double>& rhs
	);
	int CreateWMomentumLSCSR(
		const doubleField3D& Wvel0,
		double dt,
		std::vector<int>& ptr, 
		std::vector<int>& col, 
		std::vector<double>& val, 
		std::vector<double>& rhs
	);

	int CreatePressureCorrectionLSCSR(
		std::vector<int>& ptr, 
		std::vector<int>& col, 
		std::vector<double>& val, 
		std::vector<double>& rhs
	) const;

	// Set P value for all 'I' and 'J' in K defined;
	void SetPValueIJ(size_t K, double pValue);

	// Set P value for all 'J' and 'K' in I defined;
	void SetPValueJK(size_t I, double pValue);

	// Set P value for all 'K' and 'I' in J defined; 
	void SetPValueIK(size_t J, double pValue);

	// Set U value for all 'i' and 'J' in K defined;
	void SetUValueiJ(size_t K, double uValue);

	// Set U value for all 'J' and 'K' in i defined;
	void SetUValueJK(size_t i, double uValue);

	// Set U value for all 'K' and 'i' in J defined; 
	void SetUValueiK(size_t J, double uValue);

	// Set V value for all 'I' and 'K' in j defined;
	void SetVValueIK(size_t j, double vValue);

	// Set V value for all 'j' and 'I' in K defined;
	void SetVValueIj(size_t K, double vValue);

	// Set V value for all 'K' and 'j' in I defined; 
	void SetVValuejK(size_t I, double vValue);

	// Set W value for all 'I' and 'J' in k defined;
	void SetWValueIJ(size_t k, double wValue);

	// Set W value for all 'J' and 'k' in I defined;
	void SetWValueJk(size_t I, double wValue);

	// Set W value for all 'k' and 'I' in J defined; 
	void SetWValueIk(size_t J, double wValue);

	void SetPValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& KK, double pValue);
	void SetUValueInterval(const Intervalo& ii, const Intervalo& JJ, const Intervalo& KK, double uValue);
	void SetVValueInterval(const Intervalo& II, const Intervalo& jj, const Intervalo& KK, double vValue);
	void SetWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double wValue);

	void SetSpUValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue);
	void SetSuUValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue);

	void SetSpVValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue);
	void SetSuVValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue);

	void SetSpWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SpValue);
	void SetSuWValueInterval(const Intervalo& II, const Intervalo& JJ, const Intervalo& kk, double SuValue);

	void ExtrapolateForwardPJK(size_t I);
	void ExtrapolateForwardUJK(size_t i);
	void ExtrapolateForwardVjK(size_t I);
	void ExtrapolateForwardWJk(size_t I);

	void ExtrapolateBackwardPJK(size_t I);
	void ExtrapolateBackwardUJK(size_t i);
	void ExtrapolateBackwardVjK(size_t I);
	void ExtrapolateBackwardWJk(size_t I);

	void ExtrapolateBackwardWeightedUJK(size_t i, double weight);
	void ExtrapolateBackwardWeightedVjK(size_t I, double weight);
	void ExtrapolateBackwardWeightedWJk(size_t I, double weight);

	void ExtrapolateForwardIntervalPJK(size_t I, const Intervalo& JJ, const Intervalo& KK);
	void ExtrapolateForwardIntervalUJK(size_t i, const Intervalo& JJ, const Intervalo& KK);
	void ExtrapolateForwardIntervalVjK(size_t I, const Intervalo& jj, const Intervalo& KK);
	void ExtrapolateForwardIntervalWJk(size_t I, const Intervalo& JJ, const Intervalo& kk);

	void ExtrapolateBackwardIntervalPJK(size_t I, const Intervalo& JJ, const Intervalo& KK);
	void ExtrapolateBackwardIntervalUJK(size_t i, const Intervalo& JJ, const Intervalo& KK);
	void ExtrapolateBackwardIntervalVjK(size_t I, const Intervalo& jj, const Intervalo& KK);
	void ExtrapolateBackwardIntervalWJk(size_t I, const Intervalo& JJ, const Intervalo& kk);

	void ExtrapolateForwardIntervalPJI(size_t K, const Intervalo& II, const Intervalo& JJ);
	void ExtrapolateForwardIntervalUJI(size_t K, const Intervalo& ii, const Intervalo& JJ);
	void ExtrapolateForwardIntervalVJI(size_t K, const Intervalo& II, const Intervalo& jj);
	void ExtrapolateForwardIntervalWJI(size_t k, const Intervalo& II, const Intervalo& JJ);

	void ExtrapolateBackwardIntervalPJI(size_t K, const Intervalo& II, const Intervalo& JJ);
	void ExtrapolateBackwardIntervalUJI(size_t K, const Intervalo& ii, const Intervalo& JJ);
	void ExtrapolateBackwardIntervalVJI(size_t K, const Intervalo& II, const Intervalo& jj);
	void ExtrapolateBackwardIntervalWJI(size_t k, const Intervalo& II, const Intervalo& JJ);

	void SaveIJCutFieldToCSVFile(size_t K, const std::string& baseFileName) const;
	void SaveIKCutFieldToCSVFile(size_t J, const std::string& baseFileName) const;
	void SaveJKCutFieldToCSVFile(size_t I, const std::string& baseFileName) const;

	double GetFwInternalUMomentum(size_t i, size_t J, size_t K, double rho) const;
	double GetFeInternalUMomentum(size_t i, size_t J, size_t K, double rho) const;
	double GetFsInternalUMomentum(size_t i, size_t J, size_t K, double rho) const;
	double GetFnInternalUMomentum(size_t i, size_t J, size_t K, double rho) const;
	double GetFbInternalUMomentum(size_t i, size_t J, size_t K, double rho) const;
	double GetFtInternalUMomentum(size_t i, size_t J, size_t K, double rho) const;

	double GetFwInternalVMomentum(size_t I, size_t j, size_t K, double rho) const;
	double GetFeInternalVMomentum(size_t I, size_t j, size_t K, double rho) const;
	double GetFsInternalVMomentum(size_t I, size_t j, size_t K, double rho) const;
	double GetFnInternalVMomentum(size_t I, size_t j, size_t K, double rho) const;
	double GetFbInternalVMomentum(size_t I, size_t j, size_t K, double rho) const;
	double GetFtInternalVMomentum(size_t I, size_t j, size_t K, double rho) const;

	double GetFwInternalWMomentum(size_t I, size_t J, size_t k, double rho) const;
	double GetFeInternalWMomentum(size_t I, size_t J, size_t k, double rho) const;
	double GetFsInternalWMomentum(size_t I, size_t J, size_t k, double rho) const;
	double GetFnInternalWMomentum(size_t I, size_t J, size_t k, double rho) const;
	double GetFbInternalWMomentum(size_t I, size_t J, size_t k, double rho) const;
	double GetFtInternalWMomentum(size_t I, size_t J, size_t k, double rho) const;

private:
	double Psir(double r) const;

	double rep(size_t I, size_t J, size_t K, const doubleField3D& field) const;
	double ren(size_t I, size_t J, size_t K, const doubleField3D& field) const;
					 					  	
	double rwp(size_t I, size_t J, size_t K, bool ufield, const doubleField3D& field) const;
	double rwn(size_t I, size_t J, size_t K, const doubleField3D& field) const;
					 					  		  
	double rnp(size_t I, size_t J, size_t K, const doubleField3D& field) const;
	double rnn(size_t I, size_t J, size_t K, const doubleField3D& field) const;
					 					  		 
	double rsp(size_t I, size_t J, size_t K, bool vfield, const doubleField3D& field) const;
	double rsn(size_t I, size_t J, size_t K, const doubleField3D& field) const;
					 					  		 
	double rbp(size_t I, size_t J, size_t K, bool wfield, const doubleField3D& field) const;
	double rbn(size_t I, size_t J, size_t K, const doubleField3D& field) const;
					 					 		
	double rtp(size_t I, size_t J, size_t K, const doubleField3D& field) const;
	double rtn(size_t I, size_t J, size_t K, const doubleField3D& field) const;

public:
	double cdd;
	double cdA;
	double cdV;

	size_t cNX;
	size_t cNY;
	size_t cNZ;

	double cRHO;
	double cMU;
	double cD;

	doubleField3D cPressureField;
	doubleField3D cPressureCorrField;

	doubleField3D cUVelField;
	doubleField3D cVVelField;
	doubleField3D cWVelField;

	doubleField3D cSpUField;
	doubleField3D cSuUField;

	doubleField3D cSpVField;
	doubleField3D cSuVField;

	doubleField3D cSpWField;
	doubleField3D cSuWField;

	doubleField3D caijkU;
	doubleField3D caijkV;
	doubleField3D caijkW;
};