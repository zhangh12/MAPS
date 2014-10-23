
#include <R.h>
#include <R_ext/Print.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <float.h>
#include <fcntl.h>
#include <omp.h>

using namespace std;

struct SCORE{
	
	double val;
	int id;
	
	bool operator < (const SCORE& sc) const{
		return (val > sc.val);
	}
	
	SCORE(const double &v, const int &i){
		val = v;
		id = i;
	}
	
};

void rnorm(drand48_data &buf, int n, vector<double> &v1, vector<double> &v2){
	
	v1 = vector<double> (n, .0);
	v2 = vector<double> (n, .0);
	
	for(int i = 0; i < n; ++i){
		double x1, x2, w;
		double r1, r2;
		do{
			drand48_r(&buf, &r1);
			drand48_r(&buf, &r2);
			x1 = 2.0 * r1 - 1.0;
			x2 = 2.0 * r2 - 1.0;
			w = x1 * x1 + x2 * x2;
		}while(w >= 1.0);
		
		w = sqrt((-2.0 * log(w)) / w);
		v1[i] = x1 * w;
		v2[i] = x2 * w;
	}
	
}

SCORE score0(.0, -1);

extern "C" {
	
void eval_pval_opt(double * const input_obs_score_D, double *const input_obs_score_Z,  
double * const input_inv_VD, double * const input_inv_VZ, 
double * const input_rho, double * const input_kappa, 
const int * input_np, const int * input_nperm, 
const int * input_nrho, const int * input_nkappa, 
const int * input_seed, const int * input_nthread, 
double *pval, int *obs_rank){
	
	int np = *input_np;
	int nperm = *input_nperm;
	int nrho = *input_nrho;
	int nkappa = *input_nkappa;
	int seed = *input_seed;
	int nthread = *input_nthread;
	
	if(nthread <= 0 || nthread > sysconf( _SC_NPROCESSORS_ONLN )){
		nthread = sysconf( _SC_NPROCESSORS_ONLN );
	}
	
	vector<double> obs_score_D(np, .0);
	vector<double> obs_score_Z(np, .0);
	vector<vector<double> > inv_VD(np, vector<double> (np, .0));
	vector<vector<double> > inv_VZ(np, vector<double> (np, .0));
	vector<double> rho(nrho, .0);
	vector<double> kappa(nkappa, .0);
	
	int k = -1;
	for(int i = 0; i < np; ++i){
		obs_score_D[i] = input_obs_score_D[i];
		obs_score_Z[i] = input_obs_score_Z[i];
		for(int j = 0; j < np; ++j){
			++k;
			inv_VD[i][j] = input_inv_VD[k];
			inv_VZ[i][j] = input_inv_VZ[k];
		}
	}
	
	for(int i = 0; i < nrho; ++i){
		rho[i] = input_rho[i];
	}
	
	for(int i = 0; i < nkappa; ++i){
		kappa[i] = input_kappa[i];
	}
	
	vector<double> x1 (nperm+1, .0);
	vector<double> x2 (nperm+1, .0);
	vector<double> x3 (nperm+1, .0);
	
	for(int i = 0; i < np; ++i){
		x1[0] += obs_score_D[i] * obs_score_D[i];
		x2[0] += obs_score_Z[i] * obs_score_Z[i];
		x3[0] += obs_score_D[i] * obs_score_Z[i];
	}
	
	drand48_data buf;
	
	//cout << "nthread = " << nthread << endl;
	Rprintf("nthread = %d\n", nthread);
	R_FlushConsole();
	R_ProcessEvents();
	
	#pragma omp parallel num_threads(nthread) private(buf)
	{
		if(seed < 0){
			srand48_r(time(NULL), &buf);
		}else{
			srand48_r(seed, &buf);
		}
		
		#pragma omp for
		for(int k = 0; k < nperm; ++k){
			vector<double> v1, v2;
			rnorm(buf, np, v1, v2);
			
			vector<double> u1(np, .0);
			vector<double> u2(np, .0);
			
			for(int i = 0; i < np; ++i){
				for(int j = 0; j < np; ++j){
					u1[i] += inv_VD[i][j] * v1[j];
					u2[i] += inv_VZ[i][j] * v2[j];
				}
			}
			
			for(int i = 0; i < np; ++i){
				x1[k+1] += u1[i] * u1[i];
				x2[k+1] += u2[i] * u2[i];
				x3[k+1] += u1[i] * u2[i];
			}
			
		}
		
	}
	
	vector<int> rank_ref (nperm+1, nperm+2);
	vector<SCORE> score_tmp (nperm+1, score0);
	
	int l = -1;
	for(int i = 0; i < nkappa; ++i){
		for(int j = 0; j < nrho; ++j){
			++l;
			
			#pragma omp parallel num_threads(nthread)
			{
				#pragma omp for
				for(int k = 0; k < nperm+1; ++k){
					score_tmp[k].val = kappa[i] * x1[k] + (1.0-kappa[i]) * x2[k] + 2*rho[j] * sqrt(kappa[i]*(1.0-kappa[i])) * x3[k];
					score_tmp[k].id = k;
				}
			}
			
			sort(score_tmp.begin(), score_tmp.end());
			
			#pragma omp parallel num_threads(nthread)
			{
				#pragma omp for
				for(int k = 0; k < nperm +1; ++k){
					int id = score_tmp[k].id;
					if(id == 0){
						obs_rank[l] = k;
					}
					
					if(rank_ref[id] > k){
						rank_ref[id] = k;
					}
				}
			}
			
		}
	}
	
	
	*pval = .0;
	for(int k = 0; k < nperm + 1; ++k){
		if(rank_ref[k] <= rank_ref[0]){
			*pval += 1.0;
		}
	}
	
	*pval /= nperm+1;
	
}

}


