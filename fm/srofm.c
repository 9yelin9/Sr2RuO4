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
char *dirname;
char *subdirname;
char *runtime;

struct Model {
	double n0;
	double U;
	double J;
	double n[OBT];
	double m[OBT];
	double e[OBT];
	double ntot;
	double mtot;
	double etot;
	double mu;
	double itr;
};	
typedef struct Model Model;

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

void BuildH(int uord, Model *md, double k1, double k2, lapack_complex_double *h) { // Build Hamiltonian matrix
	double t2, t3, t4, ld0, *e;
	double n_sum, m_sum;
	double c1, c2, c3;

	t2  = 0.375*t1;
	t3  = 1.250*t1;
	t4  = 0.125*t1;
	ld0 = 0.200*t1;

	e    = (double*)malloc(sizeof(double) * OBT);
	e[0] = (2*t1*(cos(k1)+cos(k2)) + 4*t2*cos(k1)*cos(k2));
	e[1] = (dt + 2*t3*cos(k2) + 2*t4*cos(k1));
	e[2] = (dt + 2*t3*cos(k1) + 2*t4*cos(k2));
	c3   = (4*ld0*sin(k1)*sin(k2));

	//	   <Hamiltonian matrix>
	//   	xy		yz		zx
	// xy  	c1+c2	0		0
	// yz	0	    c1+c2	c3		+ c0
	// zx	0	    c3		c1+c2

	// c1, c2
	for(int i=0; i<OBT; i++) {
		n_sum = 0;
		m_sum = 0;

		for(int j=0; j<OBT; j++) {
			if(j != i) {
				n_sum += md->n[j];
				m_sum += md->m[j];
			}
		}

		c1 = e[i] + md->U*md->n[i] + (2*md->U-5*md->J)*n_sum;
		c2 = -((double)uord/2) * (md->U*md->m[i] + md->J*m_sum); // uord : up = 1, dn = -1 
		h[i*4] = c1 + c2;
	}

	// c3
	h[5] = h[7] = c3;

	free(e);
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

void PrintBand(Model *md) { // Print band structure data
	lapack_complex_double w[LN], vr[LDVR*LN];
	FILE *fp;
	char buf[2048];
	int p = 0;
	double k1, k2;
	double *o = (double*)malloc(sizeof(double) * 2*OBT);

	sprintf(buf, "%s/band.txt", subdirname);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("band fopen ERROR\n");
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
	fclose(fp);
}

void PrintSurface(Model *md) { // Print Fermi surface data
	lapack_complex_double w[LN], vr[LDVR*LN];
	FILE *fp;
	char buf[2048];
	double k1, k2;
	double *o = (double*)malloc(sizeof(double) * 2*OBT);

	sprintf(buf, "%s/surface.txt", subdirname);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("surface fopen ERROR\n");
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
	fclose(fp);
}

void PrintPhase(Model *md) { // Print magnetic phase data
	FILE *fp;
	char buf[2048];
	double c0, n2_sum = 0, m2_sum = 0, nn_sum = 0, mm_sum = 0;

	sprintf(buf, "%s/phase.txt", subdirname);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("phase fopen ERROR\n");
		exit(1);
	}

	// c0
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
	c0 = -N*md->U*(n2_sum-m2_sum/4) - 2*N*(md->U-2*md->J)*nn_sum + N*md->J*(nn_sum+mm_sum/4);

	fprintf(fp, "#-U/t1\tmtot\tetot\tc0\n%f\t%f\t%f\t%f\n", -md->U/t1, md->mtot, md->etot, c0);

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
	}
}

void FindM(Model *md) { // Find m converged
	FILE *fp;
	char buf[2048];
	int itr;
	double n[OBT], m[OBT], e[OBT], itv;
	double m_cvg[3] = {-100, -100, -100};

	sprintf(buf, "%s/%s.txt", dirname, runtime);
	if((fp = fopen(buf, "w")) == NULL) {
		printf("fopen ERROR\n");
		exit(1);
	}

	printf("# n = %.1f U  = %.1f J = %.3f k = %.1f\n", md->n0, md->U, md->J, (double)ks);
	fprintf(fp, "# n = %.1f U  = %.1f J = %.3f k = %.1f\n", md->n0, md->U, md->J, (double)ks);
	for(int i=0; i<OBT; i++) { 
		md->n[i] = md->n0/3;
		md->m[i] = md->n0/12;
	}

	printf("#%7s%16s%16s%16s%16s\n", "itr", "mu", "ntot", "mtot", "etot");
	fprintf(fp, "#%7s%16s%16s%16s%16s\n", "itr", "mu", "ntot", "mtot", "etot");
	for(itr=1; itr<100; itr++) {
		itv = 1;
		md->mu = 0;

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
			//printf("%f\t%f\t%f\t%f\t%f\n", md->mu, md->ntot, m[0], m[1], m[2]);
			if(fabs(md->ntot-md->n0) < 1e-3) break;

			if(md->ntot > md->n0-itv) {
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

int main(int argc, char *argv[]) {
	if(argc != 4) {
		printf("Usage : %s <n_value> <U_value> <J/U_ratio>\n", argv[0]);
		exit(1);
	}

	Model md;
	md.n0 = atof(argv[1]);
	md.U  = atof(argv[2]);
	md.J  = atof(argv[2]) * atof(argv[3]);

	dirname = (char*)malloc(sizeof(char) * 1024);
	sprintf(dirname, "/home/9yelin9/mom/fm/data/n%.1fU%.1fJ%.3fk%.1f", md.n0, md.U, md.J, (double)ks);
	mkdir(dirname, 0777);

	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	runtime = (char*)malloc(sizeof(char) * 1024);
	sprintf(runtime, "%d%d%d%d", tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min);

	//////////////////
	OptCalcEigen();
	FindM(&md);
	//////////////////

	subdirname = (char*)malloc(sizeof(char) * 1024);
	sprintf(subdirname, "%s/n%fm%fmu%fitr%.1f_%s", dirname, md.ntot, md.mtot, md.mu, md.itr, runtime);
	mkdir(subdirname, 0777);

	//////////////////
	PrintBand(&md);
	PrintSurface(&md);
	PrintPhase(&md);
	//////////////////

	free(dirname);
	free(subdirname);
	free(runtime);

	return 0;
}
