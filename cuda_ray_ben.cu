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

__global__ void timestep(double* fIn, double* fOut,double* fTemp, double* fEq,double* force, double* rho, double* T, double* ux, double* uy, double* tIn, double* tOut, double* tTemp, double *tEq, int lx, int ly, int* cxNS, int* cyNS,int* cxT, int* cyT, double* tNS,double* tT, double Thot, int Tcold, double omegaNS, double omegaT, int* oppNS, int* stmNS, int* stmT){
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
			fTemp[spd*nnodes + nid] = fIn[spd*nnodes+nid]-omegaNS*(fIn[spd*nnodes+nid]-fEq[spd*nnodes+nid])+force[spd*nnodes+nid]; //collide fIn to fIn
		}
		
		for(int i = 0; i < 5; i++){
			double cu = 3*(cxT[i]*ux[nid] + cyT[i]*uy[nid]);
			tEq[i*nnodes+nid] = T[nid]*tT[i]*(1+cu);
			tTemp[i*nnodes+nid] = tIn[i*nnodes+nid]-omegaT*(tIn[i*nnodes+nid]-tEq[i*nnodes+nid]); //collide tIn to tIn
		}

		//micro boundary fluid
		/*
		if (bottom_node || top_node){
			double temp = fIn[1*nnodes+nid];
			fIn[1*nnodes+nid] = fIn[3*nnodes+nid];
			fIn[3*nnodes+nid] = temp;
			temp = fIn[2*nnodes+nid];
			fIn[2*nnodes+nid] = fIn[4*nnodes+nid];
			fIn[4*nnodes+nid] = temp;
			temp = fIn[5*nnodes+nid];
			fIn[5*nnodes+nid] = fIn[7*nnodes+nid];
			fIn[7*nnodes+nid] = temp;
			temp = fIn[6*nnodes+nid];
			fIn[6*nnodes+nid] = fIn[8*nnodes+nid];
			fIn[8*nnodes+nid] = temp;
		}
		*/

		if (bottom_node || top_node){
			for(int i = 0; i < 9; i++){
				fTemp[i*nnodes+nid] = fIn[oppNS[i]*nnodes + nid];
			}
		}
		//streaming
		for(int i = 0; i < 9; i++){
			fOut[i*nnodes+stmNS[i*nnodes+nid]] = fTemp[i*nnodes+nid]; // this is where even/odd ping pong occurs
		}
		for(int i = 0; i < 5; i++){
			tOut[i*nnodes+(stmNS[i*nnodes+nid])] = tTemp[i*nnodes+nid];
		}
	
	}
	

}

__global__ void micro_boundary_temp(double* tOut,int lx,int nnodes, double Thot, int Tcold){
	int nid = threadIdx.x+blockIdx.x*blockDim.x;
	if (nid < nnodes){
		if (nid < lx) { //bottom node
			tOut[2*nnodes+nid] = Thot-tOut[nid]-tOut[nnodes+nid]-tOut[3*nnodes+nid]-tOut[4*nnodes+nid];		
		}

		if(nid >= (nnodes-lx)){ //topnode
			tOut[4*nnodes+nid] = Tcold-tOut[0*nnodes+nid]-tOut[1*nnodes+nid]-tOut[2*nnodes+nid]-tOut[3*nnodes+nid];	
		}
	}
}

int main(int argc, char* argv[]) {
	//give to device as arguments in kernel
	int ly = 51;//kernel arg
	int aspect_ratio = 2;
	int lx = ly*aspect_ratio; //kernel arg
	int nnodes = lx*ly; //kernel arg
	double delta_x = 1.0/(ly-2);
	double Pr = 1.0;
	double Ra = 100000;
	double gr = 0.001;
	double Thot = 1.0; //kernel arg
	int Tcold = 0; //kernel arg
	//double T0 = (Thot+Tcold)/2;
	double delta_t = sqrt(gr*delta_x);
	double nu = (sqrt(Pr/Ra)*delta_t)/(delta_x*delta_x);
	double k = sqrt(1.0/(Pr*Ra))*delta_t/(delta_x*delta_x);
	double omegaNS = 1.0/(3*nu + 0.5);//kernel arg
	double omegaT = 1.0/(3*k + 0.5);//kernel arg
	
	//host variables
	int maxT = 10000;
	int Vis_ts = 100;
	int Vis_ind = 0;
	
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
	cudaMalloc((void**)&tT_d,5*sizeof(double));
	cudaMalloc((void**)&cxT_d,5*sizeof(int));
	cudaMalloc((void**)&cyT_d,5*sizeof(int));
	
	cudaMemcpy(tT_d, tT, 5*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(cxT_d, cxT, 5*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(cyT_d, cyT, 5*sizeof(int),cudaMemcpyHostToDevice);
	
	
	
	
	


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
	double* allData = new double[3*nnodes];

	double *rho_d, *T_d, *ux_d, *uy_d, *fEven_d, *fOdd_d, *fEq_d, *force_d, *tEven_d, *tOdd_d, *tEq_d, *fTemp_d, *tTemp_d;
	
	cudaMalloc((void**)&rho_d, nnodes*sizeof(double));
	cudaMalloc((void**)&T_d, nnodes*sizeof(double));
	cudaMalloc((void**)&ux_d, nnodes*sizeof(double));
	cudaMalloc((void**)&uy_d, nnodes*sizeof(double));
	cudaMalloc((void**)&fEven_d, 9*nnodes*sizeof(double));
	cudaMalloc((void**)&fOdd_d, 9*nnodes*sizeof(double));
	cudaMalloc((void**)&fEq_d, 9*nnodes*sizeof(double));
	cudaMalloc((void**)&force_d, 9*nnodes*sizeof(double));
	cudaMalloc((void**)&tEven_d, 5*nnodes*sizeof(double));
	cudaMalloc((void**)&tOdd_d, 5*nnodes*sizeof(double));
	cudaMalloc((void**)&tEq_d, 5*nnodes*sizeof(double));
	cudaMalloc((void**)&tTemp_d, 5*nnodes*sizeof(double));
	cudaMalloc((void**)&fTemp_d, 9*nnodes*sizeof(double));

	
	cudaMemcpy(fEven_d, fEven, 9*nnodes*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(tEven_d, tEven, 5*nnodes*sizeof(double), cudaMemcpyHostToDevice);
	
	dim3 BLOCKS(128,1,1);
	dim3 GRIDS((nnodes+127)/128);
	for(int cycle = 0; cycle < maxT; cycle++){
		if(cycle%2==0){
			timestep<<<GRIDS,BLOCKS>>>(fEven_d,fOdd_d,fTemp_d, fEq_d,force_d, rho_d, T_d, ux_d, uy_d, tEven_d, tOdd_d, tTemp_d, tEq_d, lx, ly, cxNS_d, cyNS_d,cxT_d,cyT_d, tNS_d,tT_d, Thot, Tcold, omegaNS, omegaT, oppNS_d, stmNS_d,stmT_d);
			micro_boundary_temp<<<GRIDS,BLOCKS>>>(tOdd_d,lx,nnodes,Thot,Tcold);
		}
		else{
			timestep<<<GRIDS,BLOCKS>>>(fOdd_d,fEven_d,fTemp_d, fEq_d,force_d, rho_d, T_d, ux_d, uy_d, tOdd_d, tEven_d, tTemp_d, tEq_d, lx, ly, cxNS_d, cyNS_d,cxT_d,cyT_d, tNS_d,tT_d, Thot, Tcold, omegaNS, omegaT, oppNS_d, stmNS_d,stmT_d);
			micro_boundary_temp<<<GRIDS,BLOCKS>>>(tEven_d,lx,nnodes,Thot,Tcold);

		}

		if (cycle%Vis_ts==0){
			cudaMemcpy(T,T_d,nnodes*sizeof(double),cudaMemcpyDeviceToHost);
			cudaMemcpy(ux,ux_d,nnodes*sizeof(double),cudaMemcpyDeviceToHost);
			cudaMemcpy(uy,uy_d,nnodes*sizeof(double),cudaMemcpyDeviceToHost);
			for(int i = 0; i < nnodes; i++){
				allData[i] = T[i];
				allData[nnodes+i] = ux[i];
				allData[2*nnodes+i] = uy[i];
			}
			write_vector_to_file(allData, 3*nnodes, Vis_ind);
			Vis_ind++;
		}
	
	}

	std::ofstream paramfile;
	paramfile.open("params.txt");
	paramfile << lx << "\n" << ly << "\n" << Vis_ind << "\n" << delta_x;

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
	delete[] allData;

}
