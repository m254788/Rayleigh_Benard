#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

__device__ void circshift(int* out, int* in, int numrows, int numcols, int rowshift, int colshift) {
	for (int i = 0; i < numrows; i++){
		int ii = (i + rowshift)%numrows;
		if (ii<0) ii=numrows+ii;
		for (int j = 0; j < numcols; j++){
			int jj = (j + colshift) % numcols;
			if (jj<0) jj = numcols + jj;
			out[ii*numcols + jj] = in[i*numcols + j];
		}	
	}
}

__device__ void columnSum(double* sum, double* mat, int rows, int cols){
	for(int j = 0; j < cols; j++){
		double colsum = 0;
		for(int i = 0; i < rows; i++){
			colsum += mat[i*cols+j];
		}
		sum[j] = colsum;
	}
}

void write_vector_to_file(double* vector, int dim, int evolution) {
	std::ofstream myfile;
	myfile.open("evolution"+std::to_string(evolution)+".txt");
	for(int i = 0; i < dim; i++){
		myfile << vector[i] << "\n";
	}
}

int main(int argc, char* argv[]) {
	int ly = 51;
	int aspect_ratio = 2;
	int lx = ly*aspect_ratio;
	int nnodes = lx*ly;
	double delta_x = 1.0/(ly-2);
	double Pr = 1.0;
	double Ra = 1000000;
	double gr = 0.001;
	double buoyancy [] = {0, gr};
	
	double Thot = 1.0;
	int Tcold = 0;
	double T0 = (Thot+Tcold)/2;
	
	double delta_t = sqrt(gr*delta_x);
	double nu = (sqrt(Pr/Ra)*delta_t)/(delta_x*delta_x);
	double k = sqrt(1.0/(Pr*Ra))*delta_t/(delta_x*delta_x);
	double omegaNS = 1.0/(3*nu + 0.5);
	double omegaT = 1.0/(3*k + 0.5);
	
	int maxT = 10000;
	int Vis_ts = 100;
	int Vis_ind = 0;
	
	double tNS[] = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};
	int cxNS[] = {0,1,0,-1,0,1,-1,-1,1};
	int cyNS[] = {0,0,1,0,-1,1,1,-1,-1};
	int oppNS[] = {0,3,4,1,2,7,8,5,6};
	
	double tT[] = {1./3,1./6,1./6,1./6,1./6};

	int cxT[] = {0,1,0,-1,0};
	int cyT[] = {0,0,1,0,-1};
	int oppT[] = {0,3,4,1,2};
	
	//identify top and bottom nodes
	int* top_nodes = new int [lx];
	int* bottom_nodes = new int [lx];
	for(int i = 0; i < lx; i++){
		bottom_nodes[i] = i;
		top_nodes[i] = nnodes-lx+i;
	}
	
	//initialize stuff
	double* fIn = new double [9*nnodes]; //initialize fIn
	for(int s = 0; s < 9; s++){
		for(int n = 0; n < nnodes; n++){
			fIn[s*nnodes+n] = tNS[s];
		}
	}

	double* tIn = new double [5*nnodes]; //initialize tIn
	for(int s = 0; s < 5; s++){
		for(int n = 0; n < nnodes; n++){
			tIn[s*nnodes+n] = tT[s]*Tcold;
		}
	}

	for(int s = 0; s < 5; s++){ //go along bottom nodes, set equal to Thot
		for(int b = 0; b < lx; b++){ 
			tIn[s*nnodes+b] = tT[s]*Thot;
		}
	}
	
	for(int s = 0; s < 5; s++){ //create asymmetry
		tIn[s*nnodes+int(1.5*lx)] = tT[s]*(1.1*Thot);
	}
	
	//create stream target matrices
	
	int* stmNS = new int [9*nnodes];
	int* stmT = new int [5*nnodes];
	int* ind = new int [nnodes];
	int* tInd = new int [nnodes];
	for (int i = 0; i < nnodes; i++){
		ind[i] = i;
	}
	for (int i = 0; i < 9; i++){
		circshift(tInd, ind, ly, lx, -cyNS[i], -cxNS[i]);
		for(int n = 0; n < nnodes; n++){
			stmNS[i*nnodes+n] = tInd[n];
		}
		if (i < 5){
			circshift(tInd, ind, ly, lx, -cyT[i], -cxT[i]);
			for(int n = 0; n < nnodes; n++){
				stmT[i*nnodes+n] = tInd[n];
			}
		}
	}
	
	delete[] ind;
	delete[] tInd;
	double* rho = new double [nnodes];
	double* T = new double [nnodes];
	double* ux = new double [nnodes];
	double* uy = new double [nnodes];
	double* fEq = new double [9*nnodes];
	double* force = new double [9*nnodes];
	double* fOut = new double [9*nnodes];
	double cu;
	double sumx;
	double sumy;
	double* tEq = new double [5*nnodes];
	double* tOut = new double [5*nnodes];
	double* allData = new double [nnodes*3];
	
	
	double *rho_d, *T_d, *ux_d, *uy_d, *fEq_d, *force_d, *fOut_d, cu_d, sumx_d, sumy_d, *tEq_d, *tOut_d, *all_Data_d;
	
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	
	for(int cycle=0; cycle < maxT; cycle++) {
		
	}
	
	
}
