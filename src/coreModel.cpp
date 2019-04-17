#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List csFs(std::vector<double> uIn,
                std::vector<double> vIn,
                std::vector<double> wIn,
                const double& ZIn,
                const double& ustarIn,
                const double& LinvIn,
                const double& ZoIn,
                const double& bwIn,
                const double& sUustarIn,
                const double& sVustarIn,
                const double& kvIn,
                const double& C0In,
                const double& alphaIn,
                const double& MaxFetchIn,
                int outSteps)
	{
	Rcpp::RNGScope scope;

	Rcpp::NumericVector Ri(1003);
	std::vector<double> xOut;
	std::vector<double> yOut;
	std::vector<double> zOut;
	std::vector<double> wTDOut;
	std::vector<double> TimeOut;
	std::vector<int> IDOut;

	const double ustar2 = ustarIn*ustarIn;
	const double ustar4 = ustar2*ustar2;
	const double sigmaU2 = sUustarIn*sUustarIn*ustar2;
	const double sigmaV22inv = 0.5/(sVustarIn*sVustarIn*ustar2);
	const double bw2 = bwIn*bwIn;
	const double ukv = ustarIn/kvIn;
	const double cu3kv = C0In*ustar2*ukv;
	const double psiMZo = 4.8*ZoIn*LinvIn;
	const double dpsiMdz = 4.8*LinvIn;
	const double sigmaW2 = ustar2*bw2;
	const double s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);
	const double alpha2sW2 = alphaIn*2.0*sigmaW2;
	const double bisqdT = std::sqrt(alpha2sW2);
	double Time, u, v, w, x, y, z, U, dUdz, bsquare, deltaT, deltaXx, deltaXy, deltaXz, fracZ;


	const int N = uIn.size();
	//const int frac = (int)floor(N/10);
	int i;
  int j = 0;
	for(int ID = 0; ID < N; ID++){
		x = 0.0;
		y = 0.0;
		z = ZIn;
		u = uIn[ID];
		v = vIn[ID];
		w = wIn[ID];
		Time = 0.0;
		Ri = Rcpp::rnorm(1003);
		i = 0;
		while((z < 1000.0) && (abs(x) < MaxFetchIn)){
			U = ukv*(std::log(z/ZoIn) + z*dpsiMdz - psiMZo);
			dUdz = ukv*(1.0/z + dpsiMdz);
			bsquare = cu3kv*(1.0/z + 5.0*LinvIn);
			deltaT = -alpha2sW2/bsquare;

			u += (s2inv*bsquare*(sigmaW2*(u - U) + ustar2*w) + w*dUdz)*deltaT + Ri[i]*bisqdT;
			i += 1;
			v += bsquare*v*sigmaV22inv*deltaT + Ri[i]*bisqdT;
			i += 1;
			w += s2inv*bsquare*(ustar2*(u - U) + sigmaU2*w)*deltaT + Ri[i]*bisqdT;
			i += 1;
			if(i >= 1000){
				Ri = Rcpp::rnorm(1003);
				i = 0;
			}

			deltaXx = u*deltaT;
			deltaXy = v*deltaT;
			deltaXz = w*deltaT;

			if((z+deltaXz)<ZoIn){
				fracZ = (ZoIn - z)/deltaXz;
				x += deltaXx*fracZ;
				y += deltaXy*fracZ;
				Time += deltaT*fracZ;

				u = 2.0*U - u;
				v = -v;
				w = -w;

				x += u*deltaT*(1.0 - fracZ);
				y -= deltaXy*(1.0 - fracZ);
				z += 2.0*(ZoIn - z) - deltaXz;
				Time += deltaT*(1.0 - fracZ);
			} else {
				x += deltaXx;
				y += deltaXy;
				z += deltaXz;
				Time += deltaT;
			}
    if (j % outSteps == 0){
			IDOut.push_back(ID);
			xOut.push_back(x);
			yOut.push_back(y);
			zOut.push_back(z);
			TimeOut.push_back(Time);
			wTDOut.push_back(w);
    }
    j += 1;
		}
		// if((ID % frac) == 0){
		// 	Rcpp::Rcout << "\t" << (10*ID)/frac << "%" << std::endl;
		// }
		Rcpp::checkUserInterrupt();
	}

	return Rcpp::List::create(
		_["Traj_IDOut"] = IDOut,
		_["TimeOut"] = TimeOut,
		_["xOut"] = xOut,
		_["yOut"] = yOut,
		_["zOut"] = zOut,
		_["wTDOut"] = wTDOut);
}


// [[Rcpp::export]]
Rcpp::List csFi(std::vector<double> uIn,
                std::vector<double> vIn,
                std::vector<double> wIn,
                const double& ZIn,
                const double& ustarIn,
                const double& LinvIn,
                const double& ZoIn,
                const double& bwIn,
                const double& sUustarIn,
                const double& sVustarIn,
                const double& kvIn,
                const double& C0In,
                const double& alphaIn,
                const double& MaxFetchIn,
                int outSteps)
	{
	Rcpp::RNGScope scope;


	Rcpp::NumericVector Ri(1003);
	std::vector<double> xOut;
	std::vector<double> yOut;
	std::vector<double> zOut;
	std::vector<double> wTDOut;
	std::vector<double> TimeOut;
	std::vector<int> IDOut;


	const double ustar2 = ustarIn*ustarIn;
	const double ustar4 = ustar2*ustar2;
	const double sigmaU2 = sUustarIn*sUustarIn*ustar2;
	const double sigmaV22inv = 0.5/(sVustarIn*sVustarIn*ustar2);
	const double bw2 = bwIn*bwIn;
	const double x2Zo = std::sqrt(1.0-16.0*ZoIn*LinvIn);
	const double xZo = std::sqrt(x2Zo);
	const double ukv = ustarIn/kvIn;
	const double cu3kv = C0In*ustar2*ukv;
	const double psiMZo = std::log(8.0/(1.0+2.0*xZo+x2Zo)/(1.0+x2Zo)) +
	  2.0*std::atan(xZo) - PI*0.5;
	const double alpha2 = 2.0*alphaIn;
	double Time, u, v, w, x, y, z, zL, xi, xi2, psiM, phiE, dsW2dz, sigmaW2;
	double s2inv, U, dUdz, powW, bsquare, deltaT, bisqdT;
	double deltaXx, deltaXy, deltaXz, fracZ;


	const int N = uIn.size();
	// const int frac = (int)floor(uIn.size()/10);
	int i;
  int j = 0;

	for(int ID = 0; ID < N; ID++){
		x = 0.0;
		y = 0.0;
		z = ZIn;
		u = uIn[ID];
		v = vIn[ID];
		w = wIn[ID];
		Time = 0.0;
		Ri = Rcpp::rnorm(1003);
		i = 0;

		while((z < 1000.0) && (abs(x) < MaxFetchIn)){
			zL = z*LinvIn;
			xi2 = std::sqrt(1.0-16.0*zL);
			xi = std::sqrt(xi2);
			psiM = std::log(8.0/(1.0+2.0*xi+xi2)/(1.0+xi2)) + 2.0*std::atan(xi) - PI*0.5;
			powW = std::pow(1.0 - 3.0*zL,1.0/3.0);
			phiE = (bw2*bw2*(1.0 - 3.0*zL) + 1.0/powW)/((bw2*bw2 + 1.0)*std::sqrt(std::sqrt(1.0 - 6.0*zL)));
			dsW2dz = -2.0*LinvIn/powW*bw2*ustar2;

			sigmaW2 = ustar2*powW*powW*bw2;
			s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);

			U = ukv*(std::log(z/ZoIn) + psiM - psiMZo);
			dUdz = ukv/z/xi;

			bsquare = cu3kv*phiE/z;
			deltaT = -alpha2*sigmaW2/bsquare;
			bisqdT = std::sqrt(alpha2*sigmaW2);

			u += (s2inv*bsquare*(sigmaW2*(u - U) + ustar2*w) + w*dUdz)*deltaT + Ri[i]*bisqdT;
			i += 1;
			v += bsquare*v*sigmaV22inv*deltaT + Ri[i]*bisqdT;
			i += 1;
			w += (s2inv*bsquare*(ustar2*(u - U) + sigmaU2*w) + dsW2dz*(0.5 + s2inv*(ustar2*(u - U)*w + sigmaU2*w*w)))*deltaT + Ri[i]*bisqdT;
			i += 1;
			if(i >= 1000){
				Ri = Rcpp::rnorm(1003);
				i = 0;
			}

			deltaXx = u*deltaT;
			deltaXy = v*deltaT;
			deltaXz = w*deltaT;

			if((z+deltaXz)<ZoIn){
				fracZ = (ZoIn - z)/deltaXz;
				x += deltaXx*fracZ;
				y += deltaXy*fracZ;
				Time += deltaT*fracZ;

				u = 2.0*U - u;
				v = -v;
				w = -w;

				x += u*deltaT*(1.0 - fracZ);
				y -= deltaXy*(1.0 - fracZ);
				z += 2.0*(ZoIn - z) - deltaXz;
				Time += deltaT*(1.0 - fracZ);
			} else {
				x += deltaXx;
				y += deltaXy;
				z += deltaXz;
				Time += deltaT;
			}
      if (j % outSteps == 0){
  			IDOut.push_back(ID);
  			xOut.push_back(x);
  			yOut.push_back(y);
  			zOut.push_back(z);
  			TimeOut.push_back(Time);
  			wTDOut.push_back(w);
      }
      j += 1;
		}
		// if((ID % frac) == 0){
		// 	Rcpp::Rcout << "\t" << (10*ID)/frac << "%" << std::endl;
		// }
	}

	return Rcpp::List::create(
		_["Traj_IDOut"] = IDOut,
		_["TimeOut"] = TimeOut,
		_["xOut"] = xOut,
		_["yOut"] = yOut,
		_["zOut"] = zOut,
		_["wTDOut"] = wTDOut);
}


