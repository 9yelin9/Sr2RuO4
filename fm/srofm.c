// Sr4RuO2 Model - FM

#define N    6
#define OBT  3
#define LN   OBT
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
double JoverU;
double U;
double J;
char *input;
char *output;
char *runtime;

struct Model {
	double n[OBT];
	double m[OBT];
	double e[OBT];
	double ntot;
	double mtot;
	double etot;
	double c;
	double mu;
	double itr;
};	
typedef struct Model Model;

void BuildH(int uord, Model *md, double k1, double k2, lapack_complex_double *h) { // Build Hamiltonian matrix
	double t2, t3, t4, ld0, *e;
	double n_sum, m_sum;
	double a1, a2, a3;

	t2  = 0.375*t1;
	t3  = 1.250*t1;
	t4  = 0.125*t1;
	ld0 = 0.200*t1;

	e    = (double*)malloc(sizeof(double) * OBT);
	e[0] = (2*t1*(cos(k1)+cos(k2)) + 4*t2*cos(k1)*cos(k2));
	e[1] = (dt + 2*t3*cos(k2) + 2*t4*cos(k1));
	e[2] = (dt + 2*t3*cos(k1) + 2*t4*cos(k2));
	a3   = (4*ld0*sin(k1)*sin(k2));

	// <Hamiltonian matrix>
	//   	xy		yz		zx
	// xy  	a1+a2	0		0
	// yz	0	    a1+a2	a3		+ c
	// zx	0	    a3		a1+a2

	for(int i=0; i<OBT; i++) {
		n_sum = 0;
		m_sum = 0;

		for(int j=0; j<OBT; j++) {
			if(j != i) {
				n_sum += md->n[j];
				m_sum += md->m[j];
			}
		}

		a1 = e[i] + U*md->n[i] + (2*U-5*J)*n_sum;
		a2 = -((double)uord/2) * (U*md->m[i] + J*m_sum); // uord : up = 1, dn = -1 
		h[i*4] = a1 + a2;
	}

	h[5] = h[7] = a3;

	free(e);
}

void OptCalcEigen() { // Optimize CalcEigen
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double h[LDA*LN] = {0, }, w[LN], vl[LDVL*LN], vr[LDVR*LN], work;
	char jobvl = 'N', jobvr = 'V';
	double rwork[2*LN];

	LAPACK_zgeev(&jobvl, &jobvr, &ln, h, &lda, w, vl, &ldvl, vr, &ldvr, &work, &lwork, rwork, &info);
	if(info != 0){
		printf("OptCalEigen FAIL\n");
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
		printf("CalEigen FAIL\n");
		exit(1);
	}

	free(work);

	// Sorting
	int idx[2], x = 0;

	for(int i=0; i<LN; i++) {
		if((pow(creal(vr0[LN*i+1]), 2)+pow(cimag(vr0[LN*i+1]), 2)) + (pow(creal(vr0[LN*i+2]), 2)+pow(cimag(vr0[LN*i+2]), 2)) < 1e-6) {
			w[2] = w0[i];
			for(int j=0; j<LN; j++) vr[LN*2+j] = vr0[LN*i+j];
		}
		else {
			idx[x] = i;
			x++;
		}
	}

	if(creal(w0[idx[0]]) < creal(w0[idx[1]])) {
		w[0] = w0[idx[0]];
		w[1] = w0[idx[1]];
		for(int i=0; i<LN; i++) {
			vr[LN*0+i] = vr0[LN*idx[0]+i];
			vr[LN*1+i] = vr0[LN*idx[1]+i];
		}
	}
	else {
		w[0] = w0[idx[1]];
		w[1] = w0[idx[0]];
		for(int i=0; i<LN; i++) {
			vr[LN*0+i] = vr0[LN*idx[1]+i];
			vr[LN*1+i] = vr0[LN*idx[0]+i];
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

	fprintf(fp, "#p\ta_up\ta_dn\tb_up\tb_dn\tc_up\tc_dn\tNaN\n");
	for(int i=0; i<ks; i++) { // r(0, 0) ~ M(pi, 0)
		k1 = M_PI*i/(double)ks;
		k2 = 0;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j] = creal(w[j]);
		CalcEigen(-1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j+1] = creal(w[j]);
		
		fprintf(fp, "%d\t", p);
		for(int j=0; j<2*LN; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
		p++;
	}

	for(int i=0; i<ks; i++) { // M(pi, 0) ~ X(pi, pi)
		k1 = M_PI;
		k2 = M_PI*i/(double)ks;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j] = creal(w[j]);
		CalcEigen(-1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j+1] = creal(w[j]);
		
		fprintf(fp, "%d\t", p);
		for(int j=0; j<2*LN; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
		p++;
	}

	for(int i=0; i<=ks; i++) { // X(pi, pi) ~ r(0, 0)
		k1 = M_PI + M_PI*i/(double)ks;
		k2 = M_PI + M_PI*i/(double)ks;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) o[2*j] = creal(w[j]);
		CalcEigen(-1, md, k1, k2, w, vr);
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

	fprintf(fp, "#p1\tp2\ta_up\ta_dn\tb_up\tb_dn\tc_up\tc_dn\tNaN\n");
	for(int i=0; i<ks*ks; i++) {
		k1 = -M_PI + 2*M_PI*(i/ks)/(double)ks;
		k2 = -M_PI + 2*M_PI*(i%ks)/(double)ks;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<LN; j++) {
			if(creal(w[j]) < md->mu && creal(w[j]) > md->mu - 0.5) o[2*j] = 1;
			else o[2*j] = 0;
		}

		CalcEigen(-1, md, k1, k2, w, vr);
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

void CalcNME(Model *md, double *n, double *m, double *e) { // Calculate n, m, energy
	lapack_complex_double w[LN], vr[LDVR*LN];
	double k1, k2;
	double nup_sum[OBT] = {0, }, ndn_sum[OBT] = {0, }, eup_sum[OBT] = {0, }, edn_sum[OBT] = {0, };

	for(int i=0; i<ks*ks; i++) {
		k1 = -M_PI + 2*M_PI*(i/ks)/(double)ks;
		k2 = -M_PI + 2*M_PI*(i%ks)/(double)ks;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<LDVR*LN; j++) {
			if(creal(w[j/LN]) < md->mu) {
				nup_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2));
				eup_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2)) * creal(w[j/LN]);
			}
		}

		CalcEigen(-1, md, k1, k2, w, vr);
		for(int j=0; j<LDVR*LN; j++) {
			if(creal(w[j/LN]) < md->mu) {
				ndn_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2));
				edn_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2)) * creal(w[j/LN]);
			}
		}
	}

	for(int i=0; i<OBT; i++) {
		n[i] = (nup_sum[i] + ndn_sum[i])/(ks*ks);
		m[i] = (nup_sum[i] - ndn_sum[i])/(ks*ks*2);
		e[i] = (eup_sum[i] + edn_sum[i])/(ks*ks);
	}
}

void FindM(Model *md) { // Find m converged
	FILE *fp;
	char buf[2048];
	int itr;
	double n[OBT], m[OBT], e[OBT], itv;
	double m_cvg[3] = {-100, -100, -100};

	for(int i=0; i<OBT; i++) md->m[i] += 0.01*U; 

	sprintf(buf, "%s/%s.txt", input, runtime);
	printf("#%s\n", buf);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("FindM fopen ERROR\n");
		exit(1);
	}

	printf("#%7s%16s%16s%16s%16s\n", "itr", "mu", "ntot", "mtot", "etot");
	fprintf(fp, "#%7s%16s%16s%16s%16s\n", "itr", "mu", "ntot", "mtot", "etot");
	for(itr=1; itr<100; itr++) {
		itv = 1;
		md->mu = -20;

		while(itv > 1e-6) {
			md->ntot = 0;
			md->mtot = 0;
			md->etot = 0;

			CalcNME(md, n, m, e);
			for(int j=0; j<OBT; j++) {
				md->ntot += n[j];
				md->mtot += m[j];
				md->etot += e[j];
			}
			if(fabs(md->ntot-n0) < 1e-3) break;

			if(md->ntot > n0-itv) {
				md->mu -= itv;
				itv *= 0.1;
			}
			md->mu += itv;
		}
		printf("%8d%16.6f%16.6f%16.6f%16.6f\n", itr, md->mu, md->ntot, md->mtot, md->etot);
		fprintf(fp, "%8d%16.6f%16.6f%16.6f%16.6f\n", itr, md->mu, md->ntot, md->mtot, md->etot);

		for(int j=0; j<OBT; j++) {
			md->n[j] = n[j];
			md->m[j] = m[j];
			md->e[j] = e[j];
			md->itr = (double)itr;
		}

		m_cvg[itr%3] = md->mtot;
		if(fabs((m_cvg[0]+m_cvg[1]+m_cvg[2])/3 - m_cvg[itr%3]) < 1e-6) break;
	}

	fclose(fp);
}

void CalcC(Model *md) { // Calculate C
	double n2_sum = 0, m2_sum = 0, nn_sum = 0, mm_sum = 0;

	for(int i=0; i<OBT; i++) {
		n2_sum += md->n[i]*md->n[i];
		m2_sum += md->m[i]*md->m[i];

		for(int j=0; j<OBT; j++) {
			if(j != i) {
				nn_sum += md->n[j]*md->n[i];
				mm_sum += md->m[j]*md->m[i];
			}
		}
	}

	md->c = -N*U*(n2_sum-m2_sum/4) - 2*N*(U-2*J)*nn_sum + N*J*(nn_sum+mm_sum/4);
}

void PrintPhase(double U_start, double U_stop) { // Print magnetic phase data
	Model md;

	for(int i=0; i<OBT; i++) { 
		md.n[i] = n0/3;
		md.m[i] = n0/12;
	}

	for(U=U_start; U<U_stop; U+=1) {
		J = U * JoverU;

		input = (char*)malloc(sizeof(char) * 1024);
		sprintf(input, "/home/9yelin9/sro/fm/data/n%.1fU%.1fJ%.3fk%.1f", n0, U, J, (double)ks);
		mkdir(input, 0777);

		time_t t = time(NULL);
		struct tm tm = *localtime(&t);
		runtime = (char*)malloc(sizeof(char) * 1024);
		sprintf(runtime, "%d%d%d%d", tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min);

		FindM(&md);
		CalcC(&md);

		output = (char*)malloc(sizeof(char) * 1024);
		sprintf(output, "%s/n%fm%fe%fc%fmu%fitr%.1f_%s", input, md.ntot, md.mtot, md.etot, md.c, md.mu, md.itr, runtime);
		mkdir(output, 0777);

		PrintBand(&md);
		PrintSurface(&md);

		free(input);
		free(output);
		free(runtime);

		printf("\n");
	}
}

int main(int argc, char *argv[]) {
	if(argc != 5) {
		printf("Usage : %s <n value> <J/U value> <U_start value> <U_stop value>\n", argv[0]);
		exit(1);
	}

	n0     = atof(argv[1]);
	JoverU = atof(argv[2]);
	double U_start = atof(argv[3]);
	double U_stop  = atof(argv[4]);

	OptCalcEigen();
	PrintPhase(U_start, U_stop);  

	return 0;
}
