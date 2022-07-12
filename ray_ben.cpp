#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <math.h>


void circshift(int* out, int* in, int numrows, int numcols, int rowshift, int colshift){
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

void multiply(double* result, int* mat1, int row1, int col1, double* mat2, int row2, int col2){
	for(int i = 0; i < row1; i++) {
		for(int j = 0; j < col2; j++){
			int sum = 0;
			for (int k = 0; k < col1; k++){
				sum = sum + mat1[i*col1 + k]*mat2[k*col2 + j];
			}
			result[i*col2 + j] = sum;
		}
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
	
	double tNS[] = {4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36};
	int cxNS[] = {0,1,0,-1,0,1,-1,-1,1};
	int cyNS[] = {0,0,1,0,-1,1,1,-1,-1};
	int oppNS[] = {0,3,4,1,2,7,8,5,6};
	
	double tT[] = {1/3,1/6,1/6,1/6,1/6};
	int cxT[] = {0,1,0,-1,0};
	int cyT[] = {0,0,1,0,-1};
	int oppT[] = {0,3,4,1,2};
	
	/*initialize grid coordinates
	int* xcoords = new int [nnodes];
	int* ycoords = new int [nnodes];
	
	for(int r = 0; r < ly; r++){
		for(int c = 0; c < lx; c++){
			xcoords[r*lx+c] = c;
			ycoords[r*lx+c] = r;
		}
		
	}
	*/
	
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
	
	delete [] ind;
	delete[] tInd;
	double* rho = new double [nnodes];
	double* T = new double [nnodes];
	double* ux = new double [nnodes];
	double* uy = new double [nnodes];
	double* fEq = new double [9*nnodes];
	double* force = new double [9*nnodes];
	double* fOut = new double [9*nnodes];
	double cu;
	double* tEq = new double [5*nnodes];
	double* tOut = new double [5*nnodes];
	//main loop
	for(int cycle = 0; cycle < maxT; cycle++){
		
		columnSum(rho, fIn, 9, nnodes);
		columnSum(T, tIn, 5, nnodes);
		multiply(ux, cxNS, 1, 9, fIn, 9, nnodes);
		multiply(uy, cyNS, 1, 9, fIn, 9, nnodes);
		for(int i = 0; i < nnodes; i++){
			ux[i] = ux[i]/rho[i];
			uy[i] = uy[i]/rho[i];
		}
		
		//collision step fluid
		
		for(int spd = 0; spd < 9; spd++){
			for(int n = 0; n < nnodes; n++){
				cu = 3*(cxNS[spd]*ux[n] + cyNS[spd]*uy[n]);
				fEq[spd*nnodes+n] = tNS[spd]*rho[n]*(1+cu+(0.5)*cu*cu - (1.5)*(ux[n]*ux[n]+uy[n]*uy[n]));
				force[spd*nnodes+n] = 3*tNS[spd]*rho[n]*(T[n]-T0)*(cxNS[spd]*buoyancy[0]+cyNS[spd]*buoyancy[1])/(Thot-Tcold);
				fOut[spd*nnodes+n] = fIn[spd*nnodes+n] - omegaNS*(fIn[spd*nnodes+n]-fEq[spd*nnodes+n])+force[spd*nnodes+n];
			}
		}
		//collision step temperature
		
		for(int i = 0; i < 5; i++){
			for(int n = 0; n < nnodes; n++){
				cu = 3*(cxT[i]*ux[n] + cyT[i]*uy[n]);
				tEq[i*nnodes+n] = T[n]*tT[i]*(1+cu);
				tOut[i*nnodes+n] = tIn[i*nnodes+n] - omegaT*(tIn[i*nnodes+n]-tEq[i*nnodes+n]);
			}
		}
		// microscopic boundary conditions for fluid
		for (int i = 0; i < 9; i++){
			for(int n = 0; n < lx; n++){
				fOut[i*nnodes + top_nodes[n]] = fIn[oppNS[i]*nnodes + top_nodes[n]];
				fOut[i*nnodes + n] = fIn[oppNS[i]*nnodes + n]; //assuming bottom nodes are nodes 0,1,2,...,lx-1
			}
		}
		//streaming
		for(int i = 0; i < 9; i++){
			for(int n = 0; n < nnodes; n++){
				fIn[i*nnodes+stmNS[i*nnodes+n]] = fOut[i*nnodes+n];//stream fluid
				if (i < 5){
					tIn[i*nnodes+stmT[i*nnodes+n]] = tOut[i*nnodes+n];//stream temperature
				}
			}
		}

		//microscopic boundary conditions for temperature
		for(int n = 0; n < lx; n++){
			tIn[2*nnodes+n] = Thot-tIn[0*nnodes+n]-tIn[1*nnodes+n]-tIn[3*nnodes+n]-tIn[4*nnodes+n];  // bottom nodes
			tIn[4*nnodes+top_nodes[n]] = Tcold-tIn[0*nnodes+top_nodes[n]]-tIn[1*nnodes+top_nodes[n]]-tIn[2*nnodes+top_nodes[n]]-tIn[3*nnodes+top_nodes[n]]; // top nodes
		}
		
		std::cout<<cycle<<"\n";
	
	}
	delete[] top_nodes;
	delete[] bottom_nodes;
	delete[] fIn;
	delete[] tIn;
	delete[] stmNS;
	delete[] stmT;
	delete[] rho;
	delete[] T;
	delete[] ux;
	delete[] uy;
	delete[] fEq;
	delete[] force;
	delete[] fOut;
	delete[] tEq;
	delete[] tOut;
	return 0;
}

 

