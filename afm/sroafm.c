// Sr4RuO2 Model - AFM

#define N    6
#define OBT  3
#define LN  (2*OBT)
#define LDA  LN
#define LDVL LN
#define LDVR LN

#define USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <lapack.h>
#include <omp.h>
#include <memory.h>
#include <sys/stat.h>

lapack_int lwork = -1;
const int ks = 128;
const double t1 = -4.000;
const double dt =  0.000;
double n0;
double U;
double J;
char *input;
char *output;
char *runtime;

struct Model {
	double n[2][OBT];
	double m[OBT];
	double e[OBT];
	double ntot;
	double mtot;
	double etot;
	double c;
	double mu;
	int itr;
};	
typedef struct Model Model;

void BuildH(int uord, Model *md, double k1, double k2, lapack_complex_double *h) { // Build Hamiltonian matrix
	double t2, t3, t4, ld0, e[OBT];
	double n_sum, m_sum;
	double a1, a2, a3;

	t2  = 0.375*t1;
	t3  = 1.250*t1;
	t4  = 0.125*t1;
	ld0 = 0.200*t1;

	e[0] = (2*t1*(cos(k1)+cos(k2)) + 4*t2*cos(k1)*cos(k2));
	e[1] = (dt + 2*t3*cos(k2) + 2*t4*cos(k1));
	e[2] = (dt + 2*t3*cos(k1) + 2*t4*cos(k2));
	a3   = (4*ld0*sin(k1)*sin(k2));

	for(int i=0; i<OBT; i++) {
		n_sum = 0;
		m_sum = 0;

		for(int j=0; j<OBT; j++) {
			if(j != i) {
				n_sum += md->n[uord][j];
				m_sum += md->m[j];
			}
		}

		a1 = e[i] + U*md->n[uord][i] + (2*U-5*J)*n_sum;
		a2 = -((double)uord-0.5) * (U*md->m[i] + J*m_sum)/2; // uord : up = 1, dn = 0 
		h[7*i]   = a1;
		h[7*i+3] = a2;
	}

	h[8] = h[13] = a3;

	k1 += M_PI;
	k2 += M_PI;

	e[0] = (2*t1*(cos(k1)+cos(k2)) + 4*t2*cos(k1)*cos(k2));
	e[1] = (dt + 2*t3*cos(k2) + 2*t4*cos(k1));
	e[2] = (dt + 2*t3*cos(k1) + 2*t4*cos(k2));
	a3   = (4*ld0*sin(k1)*sin(k2));

	for(int i=0; i<OBT; i++) {
		n_sum = 0;
		m_sum = 0;

		for(int j=0; j<OBT; j++) {
			if(j != i) {
				n_sum += md->n[uord][j];
				m_sum += md->m[j];
			}
		}

		a1 = e[i] + U*md->n[uord][i] + (2*U-5*J)*n_sum;
		a2 = -((double)uord-0.5) * (U*md->m[i] + J*m_sum)/2; // uord : up = 1, dn = 0 
		h[7*i+21] = a1;
		h[7*i+18] = a2;
	}

	h[29] = h[34] = a3;

	/*
	for(int i=0; i<LDA*LN; i++) {
		if(i%LN == 0) printf("\n");
		printf("%f\t", creal(h[i]));
	}
	printf("\n");
	*/
}

void OptCalcEigen() { // Optimize CalcEigen
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double h[LDA*LN] = {0, }, w[LN], vl[LDVL*LN], vr[LDVR*LN], work;
	char jobvl = 'N', jobvr = 'V';
	double rwork[2*LN];

	LAPACK_zgeev(&jobvl, &jobvr, &ln, h, &lda, w, vl, &ldvl, vr, &ldvr, &work, &lwork, rwork, &info);
	if(info != 0){
		printf("OptCalcEigen FAIL\n");
		exit(1);
	}

	lwork = sizeof(work);
}

void CalcEigen(int uord, Model *md, double k1, double k2, lapack_complex_double *w, lapack_complex_double *vr) { // Calculate Eigenproblem (Using LAPACK_zgeev)
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double h[LDA*LN] = {0, }, w0[LN], vr0[LDVR*LN], vl[LDVL*LN], *work;
	char jobvl = 'N', jobvr = 'V';
	double rwork[2*LN];

	work = (lapack_complex_double*)malloc(sizeof(lapack_complex_double) * lwork);
	BuildH(uord, md, k1, k2, h);
	
	LAPACK_zgeev(&jobvl, &jobvr, &ln, h, &lda, w0, vl, &ldvl, vr0, &ldvr, work, &lwork, rwork, &info);
	if(info != 0){
		printf("CalcEigen FAIL\n");
		exit(1);
	}

	free(work);

	// Sorting
	int MAX_idx, idx1[2], idx2[2], x1 = 0, x2 = 0;
	double MAX;

	for(int i=0; i<LN; i++) {
		MAX = 0;
		MAX_idx = 0;

		for(int j=0; j<LDVR; j++) {
			if( MAX < pow(creal(vr0[LN*i+j]), 2)+pow(cimag(vr0[LN*i+j]), 2)) {
				MAX = pow(creal(vr0[LN*i+j]), 2)+pow(cimag(vr0[LN*i+j]), 2);
				MAX_idx = j;
			}
		}
		if(MAX_idx < OBT) {
			if((pow(creal(vr0[LN*i+1]), 2)+pow(cimag(vr0[LN*i+1]), 2)) + (pow(creal(vr0[LN*i+2]), 2)+pow(cimag(vr0[LN*i+2]), 2)) < 1e-6) {
				w[2] = w0[i];
				for(int j=0; j<LDVR; j++) vr[LN*2+j] = vr0[LN*i+j];
			}
			else {
				idx1[x1] = i;
				x1++;
			}
		}
		else {
			if((pow(creal(vr0[LN*i+4]), 2)+pow(cimag(vr0[LN*i+4]), 2)) + (pow(creal(vr0[LN*i+5]), 2)+pow(cimag(vr0[LN*i+5]), 2)) < 1e-6) {
				w[5] = w0[i];
				for(int j=0; j<LDVR; j++) vr[LN*5+j] = vr0[LN*i+j];
			}
			else {
				idx2[x2] = i;
				x2++;
			}
		}
	}

	if(creal(w0[idx1[0]]) < creal(w0[idx1[1]])) {
		w[0] = w0[idx1[0]];
		w[1] = w0[idx1[1]];
		for(int i=0; i<LDVR; i++) {
			vr[LN*0+i] = vr0[LN*idx1[0]+i];
			vr[LN*1+i] = vr0[LN*idx1[1]+i];
		}
	}
	else {
		w[0] = w0[idx1[1]];
		w[1] = w0[idx1[0]];
		for(int i=0; i<LDVR; i++) {
			vr[LN*0+i] = vr0[LN*idx1[1]+i];
			vr[LN*1+i] = vr0[LN*idx1[0]+i];
		}
	}

	if(creal(w0[idx2[0]]) < creal(w0[idx2[1]])) {
		w[3] = w0[idx2[0]];
		w[4] = w0[idx2[1]];
		for(int i=0; i<LDVR; i++) {
			vr[LN*3+i] = vr0[LN*idx2[0]+i];
			vr[LN*4+i] = vr0[LN*idx2[1]+i];
		}
	}
	else {
		w[3] = w0[idx2[1]];
		w[4] = w0[idx2[0]];
		for(int i=0; i<LDVR; i++) {
			vr[LN*3+i] = vr0[LN*idx2[1]+i];
			vr[LN*4+i] = vr0[LN*idx2[0]+i];
		}
	}
}

void PrintBand(Model *md) { // Print band structure data
	lapack_complex_double w[LN], vr[LDVR*LN];
	FILE *fp;
	char buf[2048];
	int p = 0;
	double k1, k2;
	double *o = (double*)malloc(sizeof(double) * 2*LN);

	sprintf(buf, "%s/band.txt", output);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("PrintBand fopen ERROR\n");
		exit(1);
	}

	fprintf(fp, "#p\ta1_up\ta1_dn\tb1_up\tb1_dn\tc1_up\tc1_dn\ta2_up\ta2_dn\tb2_up\tb2_dn\tc2_up\tc2_dn\tNaN\n");
	for(int i=0; i<ks; i++) { // r(0, 0) ~ M(pi, 0)
		k1 = M_PI*i/(double)ks;
		k2 = 0;

		CalcEigen(1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j] = creal(w[j]);
		CalcEigen(0, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j+1] = creal(w[j]);
		
		fprintf(fp, "%d\t", p);
		for(int j=0; j<2*LN; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
		p++;
	}

	for(int i=0; i<ks; i++) { // M(pi, 0) ~ X(pi, pi)
		k1 = M_PI;
		k2 = M_PI*i/(double)ks;

		CalcEigen(1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j] = creal(w[j]);
		CalcEigen(0, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j+1] = creal(w[j]);
		
		fprintf(fp, "%d\t", p);
		for(int j=0; j<2*LN; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
		p++;
	}

	for(int i=0; i<=ks; i++) { // X(pi, pi) ~ r(0, 0)
		k1 = M_PI + M_PI*i/(double)ks;
		k2 = M_PI + M_PI*i/(double)ks;

		CalcEigen(1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j] = creal(w[j]);
		CalcEigen(0, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j+1] = creal(w[j]);
		
		fprintf(fp, "%d\t", p);
		for(int j=0; j<2*LN; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
		p++;
	}

	free(o);
	fclose(fp);
}

void PrintSurface(Model *md) { // Print Fermi surface data
	lapack_complex_double w[LN], vr[LDVR*LN];
	FILE *fp;
	char buf[2048];
	double k1, k2;
	double *o = (double*)malloc(sizeof(double) * 2*LN);

	sprintf(buf, "%s/surface.txt", output);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("PrintSurface fopen ERROR\n");
		exit(1);
	}

	fprintf(fp, "#p1\tp2\ta1_up\ta1_dn\tb1_up\tb1_dn\tc1_up\tc1_dn\ta2_up\ta2_dn\tb2_up\tb2_dn\tc2_up\tc2_dn\tNaN\n");
	for(int i=0; i<ks*ks; i++) {
		k1 = -M_PI + 2*M_PI*(i/ks)/(double)ks;
		k2 = -M_PI + 2*M_PI*(i%ks)/(double)ks;

		CalcEigen(1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) {
			if(creal(w[j]) < md->mu && creal(w[j]) > md->mu - 0.5) o[2*j] = 1;
			else o[2*j] = 0;
		}

		CalcEigen(0, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) {
			if(creal(w[j]) < md->mu && creal(w[j]) > md->mu - 0.5) o[2*j+1] = 1; 
			else o[2*j+1] = 0;
		}

		fprintf(fp, "%f\t%f\t", k1, k2);
		for(int j=0; j<2*LN; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
	}
	
	free(o);
	fclose(fp);
}

void CalcNME(Model *md, double n[2][OBT], double m[OBT], double e[OBT]) { // Calculate n, m, energy
	lapack_complex_double w[LN], vr[LDVR*LN];
	double k1, k2;
	double n1up_sum[OBT] = {0, }, n1dn_sum[OBT] = {0, }, e1up_sum[OBT] = {0, }, e1dn_sum[OBT] = {0, };
	double n2up_sum[OBT] = {0, }, n2dn_sum[OBT] = {0, }, e2up_sum[OBT] = {0, }, e2dn_sum[OBT] = {0, };

	for(int i=0; i<ks*ks; i++) {
		k1 = -M_PI + 2*M_PI*(i/ks)/(double)ks;
		k2 = -M_PI + 2*M_PI*(i%ks)/(double)ks;

		CalcEigen(1, md, k1, k2, w, vr);
		for(int j=0; j<LDVR*LN; j++) {
			if(creal(w[j/LN]) < md->mu) {
				if(j%LN < OBT) {
					n1up_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2));
					e1up_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2)) * creal(w[j/LN]);
				}
				else {
					n2up_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2));
					e2up_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2)) * creal(w[j/LN]);
				}
			}
		}

		CalcEigen(0, md, k1, k2, w, vr);
		for(int j=0; j<LDVR*LN; j++) {
			if(creal(w[j/LN]) < md->mu) {
				if(j%LN < OBT) {
					n1dn_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2));
					e1dn_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2)) * creal(w[j/LN]);
				}
				else {
					n2dn_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2));
					e2dn_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2)) * creal(w[j/LN]);
				}
			}
		}
	}

	for(int i=0; i<OBT; i++) {
		n[1][i] = (n1up_sum[i] + n2up_sum[i])/(2*ks*ks);
		n[0][i] = (n1dn_sum[i] + n2dn_sum[i])/(2*ks*ks);
		printf("%f\t%f\t%f\t%f\n", n1up_sum[i], n1dn_sum[i], n2up_sum[i], n2dn_sum[i]);
	}
	printf("\n");

	for(int i=0; i<OBT; i++) {
		m[i] = ((n1up_sum[i]+n2up_sum[i]) - (n1dn_sum[i]+n2dn_sum[i]))/(2*ks*ks);
		e[i] = (e1up_sum[i] + e1dn_sum[i] + e2up_sum[i] + e2dn_sum[i])/(2*ks*ks);
	}
}

void FindM(Model *md) { // Find m converged
	FILE *fp;
	char buf[2048];
	int itr, v1, v2;
	double n[2][OBT], m[OBT], e[OBT], itv;
	double m_cvg[4] = {-100, -100, -100, -100};

	md->mu = -20;

	sprintf(buf, "%s/%s.txt", input, runtime);
	printf("#%s\n", buf);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("FindM fopen ERROR\n");
		exit(1);
	}

	printf("#%7s%16s%16s%16s%16s\n", "itr", "mu", "ntot", "mtot", "etot");
	fprintf(fp, "#%7s%16s%16s%16s%16s\n", "itr", "mu", "ntot", "mtot", "etot");
	for(itr=1; itr<50; itr++) {
		md->mu -= 1;
		itv = fabs(md->mu);
		v1 = 1;
		v2 = 1;

		while(itv > 1e-6) {
			md->ntot = 0;
			md->mtot = 0;
			md->etot = 0;

			CalcNME(md, n, m, e);
			for(int j=0; j<OBT; j++) {
				md->ntot += n[1][j] + n[0][j];
				md->mtot += m[j];
				md->etot += e[j];
			}
			//printf("%f\t%f\n", md->mu, md->ntot);
			if(fabs(md->ntot-n0) < 1e-3) break;

			if(md->ntot < n0) {
				if(v1 == 0) itv /= 2;
				md->mu += itv;
				v1 = 1;
				v2 = 0;
			}
			else {
				if(v2 == 0) itv /= 2;
				md->mu -= itv;
				v1 = 0;
				v2 = 1;
			}
		}
		printf("%8d%16.6f%16.6f%16.6f%16.6f\n", itr, md->mu, md->ntot, md->mtot, md->etot);
		fprintf(fp, "%8d%16.6f%16.6f%16.6f%16.6f\n", itr, md->mu, md->ntot, md->mtot, md->etot);

		for(int j=0; j<OBT; j++) {
			md->n[1][j] = n[1][j];
			md->n[0][j] = n[0][j];
			md->m[j] = m[j];
			md->e[j] = e[j];
			md->itr = itr;
		}

		m_cvg[itr%4] = fabs(md->mtot);
		if(fabs((m_cvg[0]+m_cvg[1]+m_cvg[2]+m_cvg[3])/4 - m_cvg[itr%4]) < 1e-3) break;
	}

	fclose(fp);
}

void CalcC(Model *md) { // Calculate C
	double n2_sum = 0, m2_sum = 0, nn_sum = 0, mm_sum = 0;

	for(int i=0; i<OBT; i++) {
		n2_sum += (md->n[1][i] + md->n[0][i])*(md->n[1][i] + md->n[0][i]);
		m2_sum += md->m[i]*md->m[i];

		for(int j=0; j<OBT; j++) {
			if(j != i) {
				nn_sum += (md->n[1][j] + md->n[0][j])*(md->n[1][i] + md->n[0][i]);
				mm_sum += md->m[j]*md->m[i];
			}
		}
	}

	md->c = -N*U*(n2_sum-m2_sum/4) - 2*N*(U-2*J)*nn_sum + N*J*(nn_sum+mm_sum/4);
}

int main(int argc, char *argv[]) {
	if(argc != 4) {
		printf("Usage : %s <n value> <U value> <J/U value>\n", argv[0]);
		exit(1);
	}

	n0 = atof(argv[1]);
	U  = atof(argv[2]);
	J  = atof(argv[2]) * atof(argv[3]);

	Model md;
	for(int i=0; i<OBT; i++) { 
		md.n[1][i] = n0/3-1.8 + 0.3*i;
		md.n[0][i] = n0/3+1.2 - 0.2*i;
		md.m[i] = n0/12 + 0.01*i;
	}

	OptCalcEigen();

	input = (char*)malloc(sizeof(char) * 1024);
	sprintf(input, "/home/9yelin9/sro/afm/data/n%.1f_U%.1f_J%.3f_k%d", n0, U, J, ks);
	mkdir(input, 0777);

	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	runtime = (char*)malloc(sizeof(char) * 1024);
	sprintf(runtime, "%d%d%d%d", tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min);

	FindM(&md);
	CalcC(&md);

	output = (char*)malloc(sizeof(char) * 1024);
	sprintf(output, "%s/nt%f_mt%f_et%f_c%f_mu%f_itr%d_v%s", input, md.ntot, md.mtot, md.etot, md.c, md.mu, md.itr, runtime);
	mkdir(output, 0777);

	PrintBand(&md);
	PrintSurface(&md);

	free(input);
	free(output);
	free(runtime);

	printf("\n");

	return 0;
}
