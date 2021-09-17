// Multi-Orbital Hubbard Model - SRO Model(FM)

#define N    6
#define OBT  3
#define LN   OBT
#define LDA  LN
#define LDVL LN
#define LDVR LN

#define USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <lapack.h>
#include <omp.h>

const int ks = 128; // k size
lapack_int lwork = -1;

struct Model {
	double U;
	double V;
	double J;
	double n[OBT];
	double m[OBT];
	double e[OBT];
	double ntot;
	double mtot;
	double etot;
	double mu;
};	
typedef struct Model Model;

void BuildH(int uord, Model *md, double k1, double k2, lapack_complex_double *h) { // Build Hamiltonian matrix
	double t1 = -4.000, t2, t3, t4, tb0, sp = 0.000, *e;
	double n_sum = 0, m_sum = 0;
	double tb, c1, c2;

	t2  = 0.375*t1;
	t3  = 1.250*t1;
	t4  = 0.125*t1;
	tb0 = 0.200*t1;

	e    = (double*)malloc(sizeof(double) * OBT);
	e[0] = (2*t1*(cos(k1)+cos(k2)) + 4*t2*cos(k1)*cos(k2));
	e[1] = (sp + 2*t3*cos(k2) + 2*t4*cos(k1));
	e[2] = (sp + 2*t3*cos(k1) + 2*t4*cos(k2));
	tb   = (4*tb0*sin(k1)*sin(k2));

	//	   <Hamiltonian matrix>
	//   	xy		yz		zx
	// xy  	c1+c2	0		0
	// yz	0	    c1+c2	tb		+ c0
	// zx	0	    tb		c1+c2

	// c1, c2
	for(int i=0; i<OBT; i++) {
		for(int j=0; j<OBT; j++) {
			if(j != i) {
				n_sum += md->n[j];
				m_sum += md->m[j];
			}
		}
		c1 = e[i] + md->U*md->n[i] + (2*md->V-md->J)*n_sum;
		c2 = -uord * (md->U*md->m[i] + md->J*m_sum)/2; // uord : up = 1, dn = -1 
		h[i*4] = c1 + c2;
	}

	// tb
	h[5] = h[7] = tb;

	free(e);
}

void OptCalcEigen() { // Optimize CalcEigen
	lapack_int ln = LN, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
	lapack_complex_double a[LDA*LN] = {0, }, w[LN], vl[LDVL*LN], vr[LDVR*LN], work;
	char jobvl = 'N', jobvr = 'V';
	double rwork[2*LN];

	LAPACK_zgeev(&jobvl, &jobvr, &ln, a, &lda, w, vl, &ldvl, vr, &ldvr, &work, &lwork, rwork, &info);
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

	for(int i=0; i<OBT; i++) {
		if((pow(creal(vr0[OBT*i+1]), 2)+pow(cimag(vr0[OBT*i+1]), 2)) + (pow(creal(vr0[OBT*i+2]), 2)+pow(cimag(vr0[OBT*i+2]), 2)) < 1e-6) {
			w[2] = w0[i];
			for(int j=0; j<OBT; j++) vr[OBT*2+j] = vr0[OBT*i+j];
		}
		else {
			idx[x] = i;
			x++;
		}
	}

	if(creal(w0[idx[0]]) < creal(w0[idx[1]])) {
		w[0] = w0[idx[0]];
		w[1] = w0[idx[1]];
		for(int i=0; i<OBT; i++) {
			vr[OBT*0+i] = vr0[OBT*idx[0]+i];
			vr[OBT*1+i] = vr0[OBT*idx[1]+i];
		}
	}
	else {
		w[0] = w0[idx[1]];
		w[1] = w0[idx[0]];
		for(int i=0; i<OBT; i++) {
			vr[OBT*0+i] = vr0[OBT*idx[1]+i];
			vr[OBT*1+i] = vr0[OBT*idx[0]+i];
		}
	}
}

void PrintBand(Model *md, int itr) { // Print band structure data
	lapack_complex_double w[LN], vr[LDVR*LN];
	FILE *fp;
	char buf[100];
	int p = 0;
	double k1, k2;
	double *o = (double*)malloc(sizeof(double) * 2*OBT);

	sprintf(buf, "data/ks%.1fU%.2fV%.2fJ%.2fn%.2fm%.2fmu%.2fitr%.1fband.txt", (double)ks, md->U, md->V, md->J, md->ntot, md->mtot, md->mu, (double)itr);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("fopen ERROR\n");
		exit(1);
	}

	fprintf(fp, "#p\ta_up\ta_dn\tb_up\tb_dn\tc_up\tc_dn\tNaN\n");
	for(int i=0; i<ks; i++) { // r(0, 0) ~ M(pi, 0)
		k1 = M_PI*i/(double)ks;
		k2 = 0;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<OBT; j++) o[2*j] = creal(w[j]);
		CalcEigen(-1, md, k1, k2, w, vr);
		for(int j=0; j<OBT; j++) o[2*j+1] = creal(w[j]);
		
		fprintf(fp, "%d\t", p);
		for(int j=0; j<2*OBT; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
		p++;
	}

	for(int i=0; i<ks; i++) { // M(pi, 0) ~ X(pi, pi)
		k1 = M_PI;
		k2 = M_PI*i/(double)ks;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<OBT; j++) o[2*j] = creal(w[j]);
		CalcEigen(-1, md, k1, k2, w, vr);
		for(int j=0; j<OBT; j++) o[2*j+1] = creal(w[j]);
		
		fprintf(fp, "%d\t", p);
		for(int j=0; j<2*OBT; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
		p++;
	}

	for(int i=0; i<=ks; i++) { // X(pi, pi) ~ r(0, 0)
		k1 = M_PI + M_PI*i/(double)ks;
		k2 = M_PI + M_PI*i/(double)ks;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<OBT; j++) o[2*j] = creal(w[j]);
		CalcEigen(-1, md, k1, k2, w, vr);
		for(int j=0; j<OBT; j++) o[2*j+1] = creal(w[j]);
		
		fprintf(fp, "%d\t", p);
		for(int j=0; j<2*OBT; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
		p++;
	}

	free(o);
}

void PrintSurface(Model *md, int itr) { // Print Fermi surface data
	lapack_complex_double w[LN], vr[LDVR*LN];
	FILE *fp;
	char buf[100];
	double k1, k2;
	double *o = (double*)malloc(sizeof(double) * 2*OBT);

	sprintf(buf, "data/ks%.1fU%.2fV%.2fJ%.2fn%.2fm%.2fmu%.2fitr%.1fsurface.txt", (double)ks, md->U, md->V, md->J, md->ntot, md->mtot, md->mu, (double)itr);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("fopen ERROR\n");
		exit(1);
	}

	fprintf(fp, "#k1\tk2\ta_up\ta_dn\tb_up\tb_dn\tc_up\tc_dn\tNaN\n");
	for(int i=0; i<ks*ks; i++) {
		k1 = -M_PI + 2*M_PI*(i/ks)/(double)ks;
		k2 = -M_PI + 2*M_PI*(i%ks)/(double)ks;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<OBT; j++) {
			if(creal(w[j]) < md->mu && creal(w[j]) > md->mu - 0.5) o[2*j] = 1;
			else o[2*j] = 0;
		}

		CalcEigen(-1, md, k1, k2, w, vr);
		for(int j=0; j<OBT; j++) {
			if(creal(w[j]) < md->mu && creal(w[j]) > md->mu - 0.5) o[2*j+1] = 1; 
			else o[2*j+1] = 0;
		}

		fprintf(fp, "%f\t%f\t", k1, k2);
		for(int j=0; j<2*OBT; j++) fprintf(fp, "%f\t", o[j]);
		fprintf(fp, "\n");
	}
	
	free(o);
}

void CalcNME(Model *md, double *n, double *m, double *e) { // Calculate n, m, energy
	lapack_complex_double w[LN], vr[LDVR*LN];
	double k1, k2;
	double nup_sum[OBT] = {0, }, ndn_sum[OBT] = {0, }, eup_sum[OBT] = {0, }, edn_sum[OBT] = {0, };
	double c0, n2_sum = 0, m2_sum = 0, nn_sum = 0, mm_sum = 0;

	for(int i=0; i<ks*ks; i++) {
		k1 = -M_PI + 2*M_PI*(i/ks)/(double)ks;
		k2 = -M_PI + 2*M_PI*(i%ks)/(double)ks;

		CalcEigen( 1, md, k1, k2, w, vr);
		for(int j=0; j<LDA*LN; j++) {
			if(creal(w[j/OBT]) < md->mu) {
				nup_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2));
				eup_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2)) * creal(w[j/OBT]);
			}
		}

		CalcEigen(-1, md, k1, k2, w, vr);
		for(int j=0; j<LDA*LN; j++) {
			if(creal(w[j/OBT]) < md->mu) {
				ndn_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2));
				edn_sum[j%OBT] += (pow(creal(vr[j]), 2) + pow(cimag(vr[j]), 2)) * creal(w[j/OBT]);
			}
		}
	}

	for(int i=0; i<OBT; i++) {
		n[i] = (nup_sum[i] + ndn_sum[i])/(ks*ks);
		m[i] = (nup_sum[i] - ndn_sum[i])/(ks*ks*2);
		e[i] = (eup_sum[i] + edn_sum[i])/(ks*ks);

		n2_sum += n[i]*n[i];
		m2_sum += m[i]*m[i];

		for(int j=0; j<OBT; j++) {
			if(j != i) {
				nn_sum += n[j]*n[i];
				mm_sum += m[j]*m[i];
			}
		}
	}
	c0 = -N*md->U*(n2_sum-m2_sum/4) - 2*N*md->V*nn_sum + N*md->J*(nn_sum+mm_sum/4);
	for(int i=0; i<OBT; i++) e[i] += c0;
}

void FindM(Model *md, double ntot_target) { // Find m converged
	int itr;
	double n[OBT], m[OBT], e[OBT], itv;
	double m_cvg[3] = {-100, -100, -100};

	printf("# ks = %.1f U  = %.2f V = %.2f J = %.2f ntot_target = %.1f\n", (double)ks, md->U, md->V, md->J, ntot_target);
	for(int i=0; i<OBT; i++) { 
		md->n[i] = ntot_target/3;
		md->m[i] = ntot_target/12;
	}

	for(itr=1; itr<100; itr++) {
		printf("#%15s%16s%16s%16s\n", "mu", "ntot", "mtot", "etot");
		itv = 0.1;
		md->mu = 5;

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
			printf("%16.6f%16.6f%16.6f%16.6f\n", md->mu, md->ntot, md->mtot, md->etot);
			if(fabs(md->ntot-ntot_target) < 1e-3) break;

			if(md->ntot > ntot_target-itv) {
				md->mu -= itv;
				itv *= 0.1;
			}
			md->mu += itv;
		}
		printf("\n# itr %d : mu = %f ntot = %f mtot = %f etot = %f\n\n", itr, md->mu, md->ntot, md->mtot, md->etot);

		for(int j=0; j<OBT; j++) {
			md->n[j] = n[j];
			md->m[j] = m[j];
			md->e[j] = e[j];
		}

		m_cvg[itr%3] = md->mtot;
		if(fabs((m_cvg[0]+m_cvg[1]+m_cvg[2])/3 - m_cvg[itr%3]) < 1e-6) break;
	}
	PrintBand(md, itr);
	PrintSurface(md, itr);
}

int main(int argc, char *argv[]) {
	if(argc != 1) {
		printf("md->Usage : %s\n", argv[0]);
		exit(1);
	}

	OptCalcEigen();

	Model md;
	md.U = md.V = md.J = 0;

	/*
	md.U = 7.0;
	md.V = 3.5;
	md.J = 1.75;
	*/

	FindM(&md, 4.0);

	return 0;
}
