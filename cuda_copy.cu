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

void columnSum(double* sum, double* mat, int rows, int cols){
	for(int j = 0; j < cols; j++){
		double colsum = 0;
		for(int i = 0; i < rows; i++){
			colsum += mat[i*cols+j];
		}
		sum[j] = colsum;
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

__global__ void timestep(double* fIn, double* fOut,double* tIn, double* tOut,int lx, int ly, double Thot, int Tcold, double omegaNS, double omegaT, int* stm){
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	int nnodes = lx*ly;
	if (tid < nnodes){
		bool bottom_node = (tid<lx);
		bool top_node = (tid>=(nnodes-lx));
		double f0,f1,f2,f3,f4,f5,f6,f7,f8;
		double t0,t1,t2,t3,t4;
		double cu;
		
		f0=fIn[tid]; f1=fIn[nnodes+tid];
    		f2=fIn[2*nnodes+tid]; f3=fIn[3*nnodes+tid];
    		f4=fIn[4*nnodes+tid]; f5=fIn[5*nnodes+tid];
    		f6=fIn[6*nnodes+tid]; f7=fIn[7*nnodes+tid];
    		f8=fIn[8*nnodes+tid];

		t0=tIn[tid]; t1=tIn[nnodes+tid];
    		t2=tIn[2*nnodes+tid]; t3=tIn[3*nnodes+tid];
    		t4=tIn[4*nnodes+tid];

		double rho = f0+f1+f2+f3+f4+f5+f6+f7+f8;
		double ux = f1-f3+f5-f6-f7+f8; ux/=rho;
		double uy = f2-f4+f5+f6-f7-f8; uy/=rho;
		double T = t0+t1+t2+t3+t4;
		
		//compute equilibrium
		double fe0,fe1,fe2,fe3,fe4,fe5,fe6,fe7,fe8;
		double te0,te1,te2,te3,te4;
		double force0,force1,force2,force3,force4,force5,force6,force7,force8;
		
		fe0= (4.0/9)*rho*(1-(1.5)*(ux*ux+uy*uy));
		force0= 0;
		te0= T/3;
		
		cu= 3*(ux);
		fe1= (1.0/9)*rho*(1+cu+(0.5)*cu*cu - (1.5)*(ux*ux+uy*uy));
		force1=0;
		te1= T*(1.0/6)*(1+cu);
		
		cu= 3*(uy);
		fe2= (1.0/9)*rho*(1+cu+(0.5)*cu*cu - (1.5)*(ux*ux+uy*uy));
		force2= 3*(1.0/9)*rho*(T-((Thot+Tcold)/2))*(0.001)/(Thot-Tcold);
		te2= T*(1.0/6)*(1+cu);
		
		cu= 3*(-ux);
		fe3= (1.0/9)*rho*(1+cu+(0.5)*cu*cu - (1.5)*(ux*ux+uy*uy));
		force3=0;
		te3= T*(1.0/6)*(1+cu);
		
		cu= 3*(-uy);
		fe4= (1.0/9)*rho*(1+cu+(0.5)*cu*cu - (1.5)*(ux*ux+uy*uy));	
		force4= 3*(1.0/9)*rho*(T-((Thot+Tcold)/2))*(-0.001)/(Thot-Tcold);
		te4= T*(1.0/6)*(1+cu);
		
		cu= 3*(ux+uy);
		fe5= (1.0/36)*rho*(1+cu+(0.5)*cu*cu - (1.5)*(ux*ux+uy*uy));
		force5= 3*(1.0/36)*rho*(T-((Thot+Tcold)/2))*(0.001)/(Thot-Tcold);
		
		cu= 3*(-ux+uy);
		fe6= (1.0/36)*rho*(1+cu+(0.5)*cu*cu - (1.5)*(ux*ux+uy*uy));
		force6= 3*(1.0/36)*rho*(T-((Thot+Tcold)/2))*(0.001)/(Thot-Tcold);

		cu= 3*(-ux-uy);
		fe7= (1.0/36)*rho*(1+cu+(0.5)*cu*cu - (1.5)*(ux*ux+uy*uy));	
		force7= 3*(1.0/36)*rho*(T-((Thot+Tcold)/2))*(-0.001)/(Thot-Tcold);
		
		cu= 3*(ux-uy);
		fe8= (1.0/36)*rho*(1+cu+(0.5)*cu*cu - (1.5)*(ux*ux+uy*uy));
		force8= 3*(1.0/36)*rho*(T-((Thot+Tcold)/2))*(-0.001)/(Thot-Tcold);
		
		//collide
		f0 = f0-omegaNS*(f0-fe0)+force0;
		f1 = f1-omegaNS*(f1-fe1)+force1;
		f2 = f2-omegaNS*(f2-fe2)+force2;
		f3 = f3-omegaNS*(f3-fe3)+force3;
		f4 = f4-omegaNS*(f4-fe4)+force4;
		f5 = f5-omegaNS*(f5-fe5)+force5;
		f6 = f6-omegaNS*(f6-fe6)+force6;
		f7 = f7-omegaNS*(f7-fe7)+force7;
		f8 = f8-omegaNS*(f8-fe8)+force8;
		
		t0 = t0-omegaT*(t0-te0);
		t1 = t1-omegaT*(t1-te1);
		t2 = t2-omegaT*(t2-te2);
		t3 = t3-omegaT*(t3-te3);
		t4 = t4-omegaT*(t4-te4);
			
		
		//micro boundary fluid
		
		if (bottom_node || top_node){
			double temp = f1;
			f1 = f3;
			f3 = temp;
			temp = f2;
			f2 = f4;
			f4 = temp;
			temp = f5;
			f5 = f7;
			f7 = temp;
			temp = f6;
			f6 = f8;
			f8 = temp;
		}
		
		//streaming
		fOut[0*nnodes+stm[0*nnodes+tid]] = f0;
		fOut[1*nnodes+stm[1*nnodes+tid]] = f1;
		fOut[2*nnodes+stm[2*nnodes+tid]] = f2;
		fOut[3*nnodes+stm[3*nnodes+tid]] = f3;
		fOut[4*nnodes+stm[4*nnodes+tid]] = f4;
		fOut[5*nnodes+stm[5*nnodes+tid]] = f5;
		fOut[6*nnodes+stm[6*nnodes+tid]] = f6;
		fOut[7*nnodes+stm[7*nnodes+tid]] = f7;
		fOut[8*nnodes+stm[8*nnodes+tid]] = f8;
		
		tOut[0*nnodes+stm[0*nnodes+tid]] = t0;
		tOut[1*nnodes+stm[1*nnodes+tid]] = t1;
		tOut[2*nnodes+stm[2*nnodes+tid]] = t2;
		tOut[3*nnodes+stm[3*nnodes+tid]] = t3;
		tOut[4*nnodes+stm[4*nnodes+tid]] = t4;
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
	int ly = 7;//kernel arg
	int aspect_ratio = 2;
	int lx = ly*aspect_ratio; //kernel arg
	int nnodes = lx*ly; //kernel arg
	double delta_x = 1.0/(ly-2);
	double Pr = 1.0;
	double Ra = 100000;
	double gr = 0.001;
	double Thot = 1.0; //kernel arg
	int Tcold = 0; //kernel arg
	double delta_t = sqrt(gr*delta_x);
	double nu = (sqrt(Pr/Ra)*delta_t)/(delta_x*delta_x);
	double k = sqrt(1.0/(Pr*Ra))*delta_t/(delta_x*delta_x);
	double omegaNS = 1.0/(3*nu + 0.5);//kernel arg
	double omegaT = 1.0/(3*k + 0.5);//kernel arg
	
	//host variables
	int maxT = 5;
	int Vis_ts = 1;
	int Vis_ind = 0;
	
	//device needs to know these values
	double tNS[] = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};
	int cxNS[] = {0,1,0,-1,0,1,-1,-1,1};//needed to easily create stream target matrix
	int cyNS[] = {0,0,1,0,-1,1,1,-1,-1};//needed to easily create stream target matrix

	

	double tT[] = {1./3,1./6,1./6,1./6,1./6};//kernel arg
	//int cxT[] = {0,1,0,-1,0};//kernel arg
	//int cyT[] = {0,0,1,0,-1};//kernel arg

	
	
	

	//initialize fEven and tEven which will then be copied into fEven_d and tEven_d
	double* fMatrix = new double [9*nnodes]; //initialize fIn
	for(int s = 0; s < 9; s++){
		for(int n = 0; n < nnodes; n++){
			fMatrix[s*nnodes+n] = tNS[s];
		}
	}
	double* tMatrix = new double [5*nnodes]; //initialize tIn
	for(int s = 0; s < 5; s++){
		for(int n = 0; n < nnodes; n++){
			tMatrix[s*nnodes+n] = tT[s]*Tcold;
		}
	}
	for(int s = 0; s < 5; s++){ //go along bottom nodes, set equal to Thot
		for(int b = 0; b < lx; b++){ 
			tMatrix[s*nnodes+b] = tT[s]*Thot;
		}
	}
	
	for(int s = 0; s < 5; s++){ //create asymmetry
		tMatrix[s*nnodes+int(1.5*lx)] = tT[s]*(1.1*Thot);
	}
	double *fEven_d, *fOdd_d, *tEven_d, *tOdd_d;
	cudaMalloc((void**)&fEven_d, 9*nnodes*sizeof(double));
	cudaMalloc((void**)&fOdd_d, 9*nnodes*sizeof(double));
	cudaMalloc((void**)&tEven_d, 5*nnodes*sizeof(double));
	cudaMalloc((void**)&tOdd_d, 5*nnodes*sizeof(double));
	cudaMemcpy(fEven_d, fMatrix, 9*nnodes*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(tEven_d, tMatrix, 5*nnodes*sizeof(double), cudaMemcpyHostToDevice);
	
	
	
	
	//create stream target matrix
	int* stm = new int [9*nnodes];
	int* ind = new int [nnodes];
	int* tInd = new int [nnodes];
	for (int i = 0; i < nnodes; i++){
		ind[i] = i;
	}
	for (int i = 0; i < 9; i++){
		circshift(tInd, ind, ly, lx, -cyNS[i], -cxNS[i]);
		for(int n = 0; n < nnodes; n++){
			stm[i*nnodes+n] = tInd[n];
		}
	}
	delete[] ind;
	delete[] tInd;
	int *stm_d;
	cudaMalloc((void**)&stm_d, 9*nnodes*sizeof(int));
	cudaMemcpy(stm_d, stm, 9*nnodes*sizeof(int), cudaMemcpyHostToDevice);

	
	//declare matrices for data viz on host
	double* T = new double[nnodes];
	double* rho = new double[nnodes];
	double* ux = new double[nnodes];
	double* uy = new double[nnodes];
	double* allData = new double[3*nnodes];

	
	

	for(int cycle = 0; cycle < maxT; cycle++){
		if(cycle%2==0){
			timestep<<<(nnodes+127)/128,128>>>(fEven_d, fOdd_d, tEven_d, tOdd_d, lx, ly, Thot, Tcold, omegaNS, omegaT, stm_d);
			micro_boundary_temp<<<(nnodes+127)/128,128>>>(tOdd_d,lx,nnodes,Thot,Tcold);
		}
		else{
			timestep<<<(nnodes+127)/128,128>>>(fOdd_d, fEven_d, tOdd_d, tEven_d, lx, ly, Thot, Tcold, omegaNS, omegaT, stm_d);
			micro_boundary_temp<<<(nnodes+127)/128,128>>>(tEven_d,lx,nnodes,Thot,Tcold);

		}

		if (cycle%Vis_ts==0){
			//std::cout<<"Executing time step " << cycle << ".\n";	
			cudaMemcpy(fMatrix,fOdd_d,9*nnodes*sizeof(double),cudaMemcpyDeviceToHost);
			cudaMemcpy(tMatrix,tOdd_d,9*nnodes*sizeof(double),cudaMemcpyDeviceToHost);
			columnSum(T,tMatrix,5,nnodes); //populates T matrix
			columnSum(rho,fMatrix,9,nnodes); //populates rho matrix
			printMat(fMatrix,9,nnodes);
			/*for(int n = 0; n < nnodes; n++){
				ux[n] = fMatrix[1*nnodes+n]-fMatrix[3*nnodes+n]+fMatrix[5*nnodes+n]-fMatrix[6*nnodes+n]-fMatrix[7*nnodes+n]+fMatrix[8*nnodes+n];
				ux[n]/=rho[n];
				uy[n] = fMatrix[2*nnodes+n]-fMatrix[4*nnodes+n]+fMatrix[5*nnodes+n]+fMatrix[6*nnodes+n]-fMatrix[7*nnodes+n]-fMatrix[8*nnodes+n];
				uy[n]/=rho[n];			
			}
			
			for(int i = 0; i < nnodes; i++){
				allData[i] = T[i];
				allData[nnodes+i] = ux[i];
				allData[2*nnodes+i] = uy[i];
			}
			write_vector_to_file(allData, 3*nnodes, Vis_ind);
			Vis_ind++;
			*/
		}
	
	}

	std::ofstream paramfile;
	paramfile.open("params.txt");
	paramfile << lx << "\n" << ly << "\n" << Vis_ind << "\n" << delta_x;

	

	cudaFree(fEven_d);
	cudaFree(fOdd_d);
	cudaFree(tEven_d);
	cudaFree(tOdd_d);
	cudaFree(stm_d);
	
	delete[] fMatrix;
	delete[] tMatrix;
	delete[] stm;	
	delete[] T;
	delete[] ux;
	delete[] uy;
	delete[] rho;
	delete[] allData;

}
