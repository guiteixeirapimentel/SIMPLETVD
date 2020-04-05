#include <iostream>
#include <ctime>

#include "Field.h"
#include "LSSolver.h"

typedef std::numeric_limits<double> doublelimits;

int main()
{
	constexpr double underRelaxPFactor = 0.3;
	constexpr double underRelaxUFactor = 0.5;
	constexpr double underRelaxVFactor = 0.5;
	constexpr double underRelaxWFactor = 0.5; 
	
	constexpr double dt = 0.0005;
	constexpr double maxTime = 1.0;
	constexpr int NITMAX = int(maxTime / dt) + 1;

	constexpr int NITMAXMOMENTUMEQTS = 500;
	constexpr int NITMAXPRESSCORREQTS = 15000;

	constexpr int NITMAXSIMPLE = 1000;

	constexpr double L = 1.0;
	constexpr double H = 0.25;
	constexpr double W = 0.25;
	
	constexpr double dd = 0.003;

	constexpr double rho = 1.0;
	constexpr double mu = 1e-5;

	constexpr double ufreestream = 1.0;
	constexpr double vfreestream = 0.0;
	constexpr double wfreestream = 0.0;

	// FOR DEBUGGING PURPOSES

	constexpr size_t NX = size_t(L / dd) + 2;
	constexpr size_t NY = size_t(W / dd) + 2;
	constexpr size_t NZ = size_t(H / dd) + 2;

	constexpr size_t NN = NX * NY * NZ;

	if (NN > 500000)
	{
		std::cout << "Atencao: o valor de pontos eh de " << NN << " sera necessario aprox. " << 0.75 * double(NN) / 500000 << "GB de memoria" << std::endl;
		std::cout << "Continuar? (y/n): ";
		char b;
		std::cin >> b;
		if (b != 'y')
			return 0;
		else
			std::cout << "Continuando....\n";
	}

	constexpr double ReL = rho * ufreestream * L / mu;
	constexpr double ReD = rho * ufreestream * H / mu;

	const double delta = 5.0 * L / sqrt(ReL);

	constexpr double Peclet = rho * ufreestream / (mu / dd);

	constexpr double CFL = 3.0 * ufreestream * dt / dd;

	// FOR DEBUGGING PURPOSES
	
	Field field(L, H, W, dd, rho, mu, ufreestream, vfreestream, wfreestream);	

	const double refPressure = 0.0;
	const int iRefPressure = (int)field.cNX / 2;
	const int jRefPressure = (int)field.cNY / 2;
	const int kRefPressure = (int)field.cNZ / 2;

	field.cPressureField[(jRefPressure * NX * NZ) + (kRefPressure * NX) + iRefPressure] = refPressure;
	
	field.SetBCFlatPlate();

#ifndef _DEBUG
	field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/initlong");
	field.SaveIJCutFieldToCSVFile(1, "OUT_DEBUG/initlat");
#endif
	
	Field fieldnm = field;

	Field field0 = field;

	std::vector<int>   ptr, col;
	std::vector<double> val, rhs;
	
	LSSolver linearSystemSolver;

	double tempoPassado = 0.0;

	for (int itTime = 0; itTime < NITMAX; itTime++)
	{
		time_t tSimple = time(nullptr);

		for (int itSimple = 0; itSimple < NITMAXSIMPLE; itSimple++)
		{
			//field.ExtrapolateBackwardIntervalUJI(NZ - 1, { 0, NX }, { 0, NY });

			{
				double uSumOut = 0.0f;
				double uSumIn = 0.0f;
				for (size_t j = 1; j < NY - 1; j++)
					for (size_t k = 1; k < NZ - 1; k++)
					{
						uSumOut += field.cUVelField[(j * NZ * NX) + ((k)*NX) + NX - 1];
						uSumIn += field.cUVelField[(j * NZ * NX) + ((k)*NX) + 1];
					}
			
				const double weight = uSumIn / uSumOut;
			
				field.ExtrapolateBackwardWeightedUJK(NX - 1, weight);
				field.ExtrapolateBackwardVjK(NX - 1);
				field.ExtrapolateBackwardWJk(NX - 1);
			}

			time_t t1 = time(nullptr);			

			double resSimple = 0.0;

			// Solve Momentum equations		

			for (int itMomentum = 0; itMomentum < NITMAXMOMENTUMEQTS; itMomentum++)
			{
				bool convergiuTudo = true;

				double resU = 0.0;

				doubleField3D Ulast = field.cUVelField;

				double resV = 0.0;

				doubleField3D Vlast = field.cVVelField;

				double resW = 0.0;

				doubleField3D Wlast = field.cWVelField;

				int nn = field.CreateUMomentumLSCSR(field0.cUVelField, dt, ptr, col, val, rhs);

				double errULS = 0.0;
				int nitULS = 0;
				// solve LS U
				field.cUVelField = linearSystemSolver.SolveSparseCRS(nn, ptr, col, val, rhs, errULS, nitULS);

				if (isnan(errULS))
				{
					std::cout << "Erro ls U - nao foi possivel resolver ls\n";
					int x = 0;
					std::cin >> x;
				}

				std::cout << "Solved U LS it " << nitULS << " residual " << errULS << std::endl;

				if ((itMomentum % 1) == 0)
				{
					double maxResUNonLinear = 0.0;
					double resUNonLinearAbs = 0.0;
					for (size_t j = 0; j < NY; j++)
						for (size_t k = 0; k < NZ; k++)
							for (size_t i = 0; i < NX; i++)
							{

								maxResUNonLinear = fmax(fabs((Ulast[(j * NX * NZ) + (k * NX) + i] - field.cUVelField[(j * NX * NZ) + (k * NX) + i])
									/ field.cUVelField[(j * NX * NZ) + (k * NX) + i]), maxResUNonLinear);
								resUNonLinearAbs = fmax(fabs((Ulast[(j * NX * NZ) + (k * NX) + i] - field.cUVelField[(j * NX * NZ) + (k * NX) + i])),
									resUNonLinearAbs);
							}

					std::cout << "max res U abs" << resUNonLinearAbs << std::endl;

					if (maxResUNonLinear < 1e-5 || resUNonLinearAbs < 1e-5)
					{
						std::cout << "Convergiu Eq. nao linear Momentum U it: " << itMomentum << "abs res: " << resUNonLinearAbs << std::endl;
					}
					else
					{
						convergiuTudo = false;
					}
				}

				nn = field.CreateVMomentumLSCSR(field0.cVVelField, dt, ptr, col, val, rhs);

				double errVLS = 0.0f;
				int nitVLS = 0;

				field.cVVelField = linearSystemSolver.SolveSparseCRS(nn, ptr, col, val, rhs, errVLS, nitVLS);

				if (isnan(errVLS))
				{
					std::cout << "Erro ls V - nao foi possivel resolver ls\n";
					int x = 0;
					std::cin >> x;
				}

				std::cout << "Solved V LS it " << nitVLS << " residual " << errVLS << std::endl;

				if ((itMomentum % 1) == 0)
				{
					double maxResVNonLinear = 0.0;
					double maxResVNonLinearAbs = 0.0;
					for (size_t j = 0; j < NY; j++)
						for (size_t k = 0; k < NZ; k++)
							for (size_t i = 0; i < NX; i++)
							{
								maxResVNonLinear = fmax(fabs((Vlast[(j * NX * NZ) + (k * NX) + i] - field.cVVelField[(j * NX * NZ) + (k * NX) + i]) / 
									field.cVVelField[(j * NX * NZ) + (k * NX) + i]), maxResVNonLinear);

								maxResVNonLinearAbs = fmax(fabs((Vlast[(j * NX * NZ) + (k * NX) + i] - field.cVVelField[(j * NX * NZ) + (k * NX) + i])), 
									maxResVNonLinearAbs);
							}


					std::cout << "max res V abs" << maxResVNonLinearAbs << std::endl;

					if (maxResVNonLinear < 1e-5 || maxResVNonLinearAbs < 1e-5)
					{
						std::cout << "Convergiu Eq. nao linear Momentum V it: " << itMomentum << " abs " << maxResVNonLinearAbs << std::endl;

					}
					else
					{
						convergiuTudo = false;
					}
				}

				nn = field.CreateWMomentumLSCSR(field0.cWVelField, dt, ptr, col, val, rhs);

				double errWLS = 0.0f;
				int nitWLS = 0;

				field.cWVelField = linearSystemSolver.SolveSparseCRS(nn, ptr, col, val, rhs, errWLS, nitWLS);

				if (isnan(errWLS))
				{
					std::cout << "Erro ls W - nao foi possivel resolver ls\n";
					int x = 0;
					std::cin >> x;
				}

				std::cout << "Solved W LS it " << nitWLS << " residual " << errWLS << std::endl;

				if ((itMomentum % 1) == 0)
				{
					double maxResWNonLinear = 0.0;
					double maxResWNonLinearAbs = 0.0;
					for (size_t j = 0; j < NY; j++)
						for (size_t k = 0; k < NZ; k++)
							for (size_t i = 0; i < NX; i++)
							{
								maxResWNonLinear = fmax(fabs((Wlast[(j * NX * NZ) + (k * NX) + i] - field.cWVelField[(j * NX * NZ) + (k * NX) + i])
									/ field.cWVelField[(j * NX * NZ) + (k * NX) + i]), maxResWNonLinear);
								maxResWNonLinearAbs = fmax(fabs((Wlast[(j * NX * NZ) + (k * NX) + i] - field.cWVelField[(j * NX * NZ) + (k * NX) + i])),
									maxResWNonLinearAbs);
							}

					std::cout << "max res W abs" << maxResWNonLinearAbs << std::endl;

					if (maxResWNonLinear < 1e-5 || maxResWNonLinearAbs < 1e-5)
					{
						std::cout << "Convergiu Eq. nao linear Momentum W it: " << itMomentum << " abs " << maxResWNonLinearAbs << std::endl;

					}
					else
					{
						convergiuTudo = false;

					}
				}

				if (convergiuTudo)
					break;
			}


#ifdef _DEBUG
			field.SaveIKCutFieldToCSVFile(NY / 2, "OUT_DEBUG/posmomlong");
			field.SaveJKCutFieldToCSVFile(NX / 2, "OUT_DEBUG/posmomtrans");
#endif

			//field.SaveIKCutFieldToCSVFile(NYy / 2, "posMom");
			//field.SaveIJCutFieldToCSVFile(NZz / 2, "posMomTrans");
			int x = 0;
			//std::cin >> x;				   

			// Solve Pressure Correction Equations
			{
				int nn = field.CreatePressureCorrectionLSCSR(ptr, col, val, rhs);

				double errPCLS = 0.0;
				int nitPCLS = 0;

				field.cPressureCorrField = linearSystemSolver.SolveSparseCRS(nn, ptr, col, val, rhs, errPCLS, nitPCLS);

				if (isnan(errPCLS))
				{
					std::cout << "Erro ls PC - nao foi possivel resolver ls\n";
					int x = 0;
					std::cin >> x;
				}

				std::cout << "Solved pressure corr LS it " << nitPCLS << " residual " << errPCLS << std::endl;
				
			}


			// Correct
			// pressure
			for (int index = 0; index < NY * NX * NZ; index++)
			{
				field.cPressureField[index] = field.cPressureField[index] +
					((field.cPressureCorrField[index]) * underRelaxPFactor);
			}

			// velocities
			// u
			for (int J = 1; J < NY - 1; J++)
			{
				for (int K = 1; K < NZ - 1; K++)
				{
					for (int i = 2; i < NX - 1; i++)
					{
						field.cUVelField[(J * NZ * NX) + (K * NX) + i] += (dd * dd / field.caijkU[(J * NZ * NX) + (K * NX) + i]) *
							(field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i - 1] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);

						field.cUVelField[(J * NZ * NX) + (K * NX) + i] = (field.cUVelField[(J * NZ * NX) + (K * NX) + i] * underRelaxUFactor) +
							((1.0f - underRelaxUFactor) * fieldnm.cUVelField[(J * NZ * NX) + (K * NX) + i]);
					}
				}
			}
			//v
			for (int J = 2; J < NY - 1; J++)
			{
				for (int K = 1; K < NZ - 1; K++)
				{
					for (int i = 1; i < NX - 1; i++)
					{
						field.cVVelField[(J * NZ * NX) + (K * NX) + i] += (dd* dd/ field.caijkV[(J * NZ * NX) + (K * NX) + i]) *
							(field.cPressureCorrField[((J-1) * NZ * NX) + (K * NX) + i] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);

						field.cVVelField[(J * NZ * NX) + (K * NX) + i] = (field.cVVelField[(J * NZ * NX) + (K * NX) + i] * underRelaxVFactor) +
							((1.0f - underRelaxVFactor) * fieldnm.cVVelField[(J * NZ * NX) + (K * NX) + i]);
					}
				}
			}
			//w
			for (int J = 1; J < NY - 1; J++)
			{
				for (int K = 2; K < NZ - 1; K++)
				{
					for (int i = 1; i < NX - 1; i++)
					{
						field.cWVelField[(J * NZ * NX) + (K * NX) + i] += (dd* dd/ field.caijkW[(J * NZ * NX) + (K * NX) + i]) *
							(field.cPressureCorrField[(J * NZ * NX) + ((K-1) * NX) + i] - field.cPressureCorrField[(J * NZ * NX) + (K * NX) + i]);

						field.cWVelField[(J * NZ * NX) + (K * NX) + i] = (field.cWVelField[(J * NZ * NX) + (K * NX) + i] * underRelaxWFactor) +
							((1.0f - underRelaxWFactor) * fieldnm.cWVelField[(J * NZ * NX) + (K * NX) + i]);
						
					}
				}
			}

			double ures = 0.0;
			double vres = 0.0;
			double wres = 0.0;
			double pres = 0.0;

			// calc residuals velocities
			for (size_t j = 0; j < NY; j++)
				for (size_t k = 0; k < NZ; k++)
					for (size_t i = 0; i < NX; i++)
					{						
						ures = fmax(fabs(field.cUVelField[(j * NZ * NX) + (k * NX) + i] - fieldnm.cUVelField[(j * NZ * NX) + (k * NX) + i]), ures);
						vres = fmax(fabs(field.cVVelField[(j * NZ * NX) + (k * NX) + i] - fieldnm.cVVelField[(j * NZ * NX) + (k * NX) + i]), vres);
						wres = fmax(fabs(field.cWVelField[(j * NZ * NX) + (k * NX) + i] - fieldnm.cWVelField[(j * NZ * NX) + (k * NX) + i]), wres);
					}

			for (size_t j = 0; j < NY; j++)
				for (size_t k = 0; k < NZ; k++)
					for (size_t i = 0; i < NX; i++)
					{
						pres = fmax(fabs(field.cPressureField[(j * NZ * NX) + (k * NX) + i] - 
							fieldnm.cPressureField[(j * NZ * NX) + (k * NX) + i]), pres);
					}

			fieldnm = field;

			for(int j = 1; j < field.cNY - 1; j++)
				for(int k = 1; k < field.cNZ - 1; k++)
					for (int i = 1; i < field.cNX - 1; i++)
					{
						resSimple = fmax(
							fabs(
								 (
									field.cUVelField[(j * NZ * NX) + (k * NX) + i] - field.cUVelField[(j * NZ * NX) + (k * NX) + i + 1] +
							field.cVVelField[(j * NZ * NX) + (k * NX) + i] - field.cVVelField[((j+1) * NZ * NX) + (k * NX) + i] + field.cWVelField[(j * NZ * NX) + (k * NX) + i] -
									field.cWVelField[(j * NZ * NX) + ((k+1) * NX) + i]
									)
							), resSimple);
					}

			std::cout << "##########################################\n";
			std::cout << "##########################################\n";

			std::cout << "resSimple: " << resSimple << " ures "
				      << ures << " vres " << vres << " wres " << wres << " pres " << pres << std::endl;

			if ((resSimple < 1e-6 && pres < 1e-4) /*|| (ures < 1e-6 && vres < 1e-6 && wres < 1e-6 && pres < 1e-6)*/)
			{
				std::cout << "Simple converged; time: " << tempoPassado << std::endl;
				std::cout << "Tempo para convergir simple: " << time(nullptr) - tSimple << " segundos \n";

				std::cout << "Tempo aprox restante para conclusao : " << (NITMAX - itTime) * (time(nullptr) - tSimple) << " segs" << std::endl;

				std::cout << "##########################################\n";
				std::cout << "##########################################\n";

				field.SaveIKCutFieldToCSVFile(NY / 2, "turbFlatPlate/corteLong" /*+ std::to_string(tempoPassado)*/);
				field.SaveJKCutFieldToCSVFile(NX / 2, "turbFlatPlate/corteTrans"/* + std::to_string(tempoPassado)*/);

				break;
			}



			std::cout << "Last simple iteration spent " << time(nullptr) - t1 << " seconds\n";
		}

		fieldnm = field;
		field0 = field;

		tempoPassado += dt;
	}

	return 0;
}