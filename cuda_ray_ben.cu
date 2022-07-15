#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

void circshift(int* out, int* in, int numrows, int numcols, int rowshift, int colshift) {
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
/*
__device__ void columnSum(double* sum, double* mat, int rows, int cols){
	for(int j = 0; j < cols; j++){
		double colsum = 0;
		for(int i = 0; i < rows; i++){
			colsum += mat[i*cols+j];
		}
		sum[j] = colsum;
	}
}
*/

void printMat(double* mat, int rows, int cols){
        for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                        std::cout << mat[i*cols+j]<<" ";
                }
                std::cout<<"\n";
        }
	std::cout<<"\n";
}
void write_vector_to_file(double* vector, int dim, int evolution) {
	std::ofstream myfile;
	myfile.open("evolution"+std::to_string(evolution)+".txt");
	for(int i = 0; i < dim; i++){
		myfile << vector[i] << "\n";
	}
}

__global__ void timestep(double* fIn, double* fOut,double* fEq,double* force, double* rho, double* T, double* ux, double* uy, double* tIn, double* tOut, double *tEq, int lx, int ly, int* cxNS, int* cyNS,int* cxT, int* cyT, double* tNS,double* tT, int Thot, int Tcold, double omegaNS, double omegaT, int* oppNS, int* stmNS, int* stmT){
	int nid = threadIdx.x + blockIdx.x*blockDim.x;
	int nnodes = lx*ly;
	if (nid < nnodes){
		bool bottom_node = (nid<lx);
		bool top_node = (nid>=(nnodes-lx));

		rho[nid] = 0;
		for(int i = 0; i < 9; i++){
			rho[nid] += fIn[i*nnodes + nid];
		}
		T[nid] = 0;
		for(int i = 0; i < 5; i++){
			T[nid] += tIn[i*nnodes + nid];
		}

		ux[nid] = (fIn[1*nnodes+nid]-fIn[3*nnodes+nid]+fIn[5*nnodes+nid]-fIn[6*nnodes+nid]-fIn[7*nnodes+nid]+fIn[8*nnodes+nid])/rho[nid];
		uy[nid] = (fIn[2*nnodes+nid]-fIn[4*nnodes+nid]+fIn[5*nnodes+nid]+fIn[6*nnodes+nid]-fIn[7*nnodes+nid]-fIn[8*nnodes+nid])/rho[nid];
		//collision
		for(int spd = 0; spd < 9; spd++){
			double cu = 3*(cxNS[spd]*ux[nid] + cyNS[spd]*uy[nid]);
			fEq[spd*nnodes + nid] = tNS[spd]*rho[nid]*(1+cu+(0.5)*cu*cu - (1.5)*(ux[nid]*ux[nid]+uy[nid]*uy[nid]));
			force[spd*nnodes + nid] = 3*tNS[spd]*rho[nid]*(T[nid]-((Thot+Tcold)/2))*(cyNS[spd]*0.001)/(Thot-Tcold);
			fOut[spd*nnodes + nid] = fIn[spd*nnodes+nid]-omegaNS*(fIn[spd*nnodes+nid]-fEq[spd*nnodes+nid])+force[spd*nnodes+nid];
		}
		
		for(int i = 0; i < 5; i++){
			double cu = 3*(cxT[i]*ux[nid] + cyT[i]*uy[nid]);
			tEq[i*nnodes+nid] = T[nid]*tT[i]*(1+cu);
			tOut[i*nnodes+nid] = tIn[i*nnodes+nid]-omegaT*(tIn[i*nnodes+nid]-tEq[i*nnodes+nid]);
		}

		//micro boundary fluid
		if (bottom_node || top_node){
			for(int i = 0; i < 9; i++){
				fOut[i*nnodes+nid] = fIn[oppNS[i]*nnodes + nid];
			}
		}

		//streaming
		for(int i = 0; i < 9; i++){
			fIn[i*nnodes+stmNS[i*nnodes+nid]] = fOut[i*nnodes+nid];
		}
		for(int i = 0; i < 5; i++){
			tIn[i*nnodes+(stmNS[i*nnodes+nid])] = tOut[i*nnodes+nid];
		}

		//micro boundary temp
		if(top_node){
			tIn[4*nnodes+nid] = Tcold-tIn[0*nnodes+nid]-tIn[1*nnodes+nid]-tIn[2*nnodes+nid]-tIn[3*nnodes+nid];
		}
		if(bottom_node){
			tIn[2*nnodes+nid] = Thot-tIn[0*nnodes+nid]-tIn[1*nnodes+nid]-tIn[3*nnodes+nid]-tIn[4*nnodes+nid];
		}
	
	}

}


int main(int argc, char* argv[]) {
	//give to device as arguments in kernel
	int ly = 3;//kernel arg
	int aspect_ratio = 2;
	int lx = ly*aspect_ratio; //kernel arg
	int nnodes = lx*ly; //kernel arg
	double delta_x = 1.0/(ly-2);
	double Pr = 1.0;
	double Ra = 1000000;
	double gr = 0.001;
	double Thot = 1.0; //kernel arg
	int Tcold = 0; //kernel arg
	//double T0 = (Thot+Tcold)/2;
	double delta_t = sqrt(gr*delta_x);
	double nu = (sqrt(Pr/Ra)*delta_t)/(delta_x*delta_x);
	double k = sqrt(1.0/(Pr*Ra))*delta_t/(delta_x*delta_x);
	double omegaNS = 1.0/(3*nu + 0.5);//kernel arg
	double omegaT = 1.0/(3*k + 0.5);//kernel arg
	
	/*host variables
	int maxT = 10000;
	int Vis_ts = 100;
	int Vis_ind = 0;
	*/
	//device needs to know these values
	double tNS[] = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};//kernel arg
	int cxNS[] = {0,1,0,-1,0,1,-1,-1,1};//kernel arg
	int cyNS[] = {0,0,1,0,-1,1,1,-1,-1};//kernel arg
	int oppNS[] = {0,3,4,1,2,7,8,5,6};//kernel arg

	
	double *tNS_d;
	int *cxNS_d, *cyNS_d, *oppNS_d;
	cudaMalloc((void**)&tNS_d,9*sizeof(double));
	cudaMalloc((void**)&cxNS_d, 9*sizeof(int));
	cudaMalloc((void**)&cyNS_d, 9*sizeof(int));
	cudaMalloc((void**)&oppNS_d, 9*sizeof(int));
	cudaMemcpy(tNS_d, tNS, 9*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(cxNS_d, cxNS, 9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(cyNS_d, cyNS, 9*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(oppNS_d, oppNS, 9*sizeof(int),cudaMemcpyHostToDevice);
	

	double tT[] = {1./3,1./6,1./6,1./6,1./6};//kernel arg
	int cxT[] = {0,1,0,-1,0};//kernel arg
	int cyT[] = {0,0,1,0,-1};//kernel arg
	//int oppT[] = {0,3,4,1,2};//kernel arg
	
	double *tT_d;
	int *cxT_d, *cyT_d;
	cudaMalloc(&tT_d,5*sizeof(double));
	cudaMalloc(&cxT_d,5*sizeof(int));
	cudaMalloc(&cyT_d,5*sizeof(int));
	
	cudaMemcpy(tT_d, tT, 5*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(cxT_d, cxT, 5*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(cyT_d, cyT, 5*sizeof(int),cudaMemcpyHostToDevice);
	
	
	
	
	
	/*identify top and bottom nodes
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
	*/

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
	
	
	int *stmNS_d, *stmT_d;
	cudaMalloc(&stmNS_d, 9*nnodes*sizeof(int));
	cudaMalloc(&stmT_d, 9*nnodes*sizeof(int));
	cudaMemcpy(stmNS_d, stmNS, 9*nnodes*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(stmT_d, stmT, 5*nnodes*sizeof(int), cudaMemcpyHostToDevice);
	
	
	double* T = new double[nnodes];
	double* ux = new double[nnodes];
	double* uy = new double[nnodes];
	double *rho_d, *T_d, *ux_d, *uy_d, *fEven_d, *fOdd_d, *fEq_d, *force_d, *tEven_d, *tOdd_d, *tEq_d;
	
	cudaMalloc(&rho_d, nnodes*sizeof(double));
	cudaMalloc(&T_d, nnodes*sizeof(double));
	cudaMalloc(&ux_d, nnodes*sizeof(double));
	cudaMalloc(&uy_d, nnodes*sizeof(double));
	cudaMalloc(&fEven_d, 9*nnodes*sizeof(double));
	cudaMalloc(&fOdd_d, 9*nnodes*sizeof(double));
	cudaMalloc(&fEq_d, 9*nnodes*sizeof(double));
	cudaMalloc(&force_d, 9*nnodes*sizeof(double));
	cudaMalloc(&tEven_d, 5*nnodes*sizeof(double));
	cudaMalloc(&tOdd_d, 5*nnodes*sizeof(double));
	cudaMalloc(&tEq_d, 5*nnodes*sizeof(double));

	
	cudaMemcpy(fEven_d, fEven, 9*nnodes*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(tEven_d, tEven, 5*nnodes*sizeof(double), cudaMemcpyHostToDevice);
	double* test = new double[9*nnodes];	
	timestep<<<((nnodes+127)/128),128>>>(fEven_d,fOdd_d,fEq_d,force_d, rho_d, T_d, ux_d, uy_d, tEven_d, tOdd_d, tEq_d, lx, ly, cxNS_d, cyNS_d,cxT_d,cyT_d, tNS_d,tT_d, Thot, Tcold, omegaNS, omegaT, oppNS_d, stmNS_d,stmT_d);
	cudaMemcpy(test,fEven_d,9*nnodes*sizeof(double),cudaMemcpyDeviceToHost);
	printMat(test,9,nnodes);

	cudaFree(rho_d);
	cudaFree(T_d);
	cudaFree(ux_d);
	cudaFree(uy_d);
	cudaFree(fEven_d);
	cudaFree(fOdd_d);
	cudaFree(fEq_d);
	cudaFree(force_d);
	cudaFree(tEven_d);
	cudaFree(tOdd_d);
	cudaFree(tEq_d);

	delete[] fEven;
	delete[] tEven;
	delete[] stmNS;
	delete[] stmT;	
	delete[] T;
	delete[] ux;
	delete[] uy;
	delete[] test;
}
