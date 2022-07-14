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
	//give to device as arguments in kernel
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
	
	//host variables
	int maxT = 10000;
	int Vis_ts = 100;
	int Vis_ind = 0;
	
	//device needs to know these values
	double tNS[] = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};
	int cxNS[] = {0,1,0,-1,0,1,-1,-1,1};
	int cyNS[] = {0,0,1,0,-1,1,1,-1,-1};
	int oppNS[] = {0,3,4,1,2,7,8,5,6};

	double *tNS_d;
	int *cxNS_d, *cyNS_d, *oppNS_d;
	cudaMalloc((void**)&tNS_d,9*sizeof(double));
	cudaMalloc((void**)&cxNS_d, 9*sizeof(int));
	cudaMalloc((void**)&cyNS_d, 9*sizeof(int));
	cudaMalloc((void**)&oppNS_d, 9*sizeof(int));
	cudaMemcpy(tNS_d, tNS, cudaMemcpyHostToDevice);
	cudaMemcpy(cxNS_d, cxNS, cudaMemcpyHostToDevice);
	cudaMemcpy(cyNS_d, cyNS, cudaMemcpyHostToDevice);
	cudaMemcpy(oppNS_d, oppNS, cudaMemcpyHostToDevice);
	
	double tT[] = {1./3,1./6,1./6,1./6,1./6};
	int cxT[] = {0,1,0,-1,0};
	int cyT[] = {0,0,1,0,-1};
	int oppT[] = {0,3,4,1,2};
	
	double *tT_d;
	int *cxT_d, *cyT_d, oppT_d;
	cudaMalloc((void**)&tT_d,5*sizeof(double));
	cudaMalloc((void**)&cxT_d,5*sizeof(int));
	cudaMalloc((void**)&cyT_d,5*sizeof(int));
	cudaMalloc((void**)&oppT_d,5*sizeof(int));
	cudaMemcpy(tT_d, tT, cudaMemcpyHostToDevice);
	cudaMemcpy(cxT_d, cxT, cudaMemcpyHostToDevice);
	cudaMemcpy(cyT_d, cyT, cudaMemcpyHostToDevice);
	cudaMemcpy(oppT_d, oppT, cudaMemcpyHostToDevice);
	
	
	
	
	//identify top and bottom nodes
	int* top_nodes = new int [lx];
	int* bottom_nodes = new int [lx];
	for(int i = 0; i < lx; i++){
		bottom_nodes[i] = i;
		top_nodes[i] = nnodes-lx+i;
	}
	
	int *top_nodes_d, *bottom_nodes_d;
	cudaMalloc((void**)&top_nodes_d,lx*sizeof(int));
	cudaMalloc((void**)&bottom_nodes_d,lx*sizeof(int));
	cudaMemcpy(top_nodes_d,top_nodes,lx*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(bottom_nodes_d,bottom_nodes,lx*sizeof(int),cudaMemcpyHostToDevice);
	
	//initialize stuff
	double* fEven = new double [9*nnodes]; //initialize fIn
	for(int s = 0; s < 9; s++){
		for(int n = 0; n < nnodes; n++){
			fEven[s*nnodes+n] = tNS[s];
		}
	}

	double* tEven = new double [5*nnodes]; //initialize tIn
	for(int s = 0; s < 5; s++){
		for(int n = 0; n < nnodes; n++){
			tEven[s*nnodes+n] = tT[s]*Tcold;
		}
	}

	for(int s = 0; s < 5; s++){ //go along bottom nodes, set equal to Thot
		for(int b = 0; b < lx; b++){ 
			tEven[s*nnodes+b] = tT[s]*Thot;
		}
	}
	
	for(int s = 0; s < 5; s++){ //create asymmetry
		tEven[s*nnodes+int(1.5*lx)] = tT[s]*(1.1*Thot);
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
	
	double *stmNS_d, *stmT_d;
	cudaMalloc(&stmNS, 9*nnodes*sizeof(double));
	cudaMalloc(&stmT, 9*nnodes*sizeof(double));
	cudaMemcpy(stmNS_d, stmNS, 9*nnodes*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(stmT_d, stmT, 5*nnodes*sizeof(double), cudaMemcpyHostToDevice);
	
	double* allData = new double [3*nnodes];

	double *rho_d, *T_d, *ux_d, *uy_d, *fEven_d, *fOdd_d, *fEq_d, *force_d, *tEven_d, *tOdd_d, *tEq_d, *allData_d;
	
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&T_d, nnodes*sizeof(double));
	cudaMalloc(&ux_d, nnodes*sizeof(double));
	cudaMalloc(&uy_d, nnodes*sizeof(double));
	cudaMalloc(&fEven_d, 9*nnodes*sizeof(double));
	cudaMalloc(&fEq_d, 9*nnodes*sizeof(double));
	cudaMalloc(&force_d, 9*nnodes*sizeof(double));
	cudaMalloc(&fOdd_d, 9*nnodes*sizeof(double));
	cudaMalloc(&tEq_d, 5*nnodes*sizeof(double));
	cudaMalloc(&tEven_d, 5*nnodes*sizeof(double));
	cudaMalloc(&tOdd_d, 5*nnodes*sizeof(double));
	cudaMalloc(&allData_d, 3*nnodes*sizeof(double));
	
	cudaMemcpy(fEven_d, fEven, 9*nnodes*sizeof(double));
	cudaMemcpy(tEven_d, tEven, 5*nnodes*sizeof(double));
	
	
	
	
	
}
