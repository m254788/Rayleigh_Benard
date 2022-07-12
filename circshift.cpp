#include <iostream>

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

void printMat(int* mat, int rows, int cols){
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			std::cout << mat[i*cols+j]<<" ";
		}
		std::cout<<"\n";
	}
}

int main(int argc, char* argv[]){
	int* matrix = new int[10];
	for (int i = 0; i < 10; i++){
		matrix[i] = i;
	}
	
	int* newmatrix = new int[10];
	printMat(matrix, 2, 5);
	circshift(newmatrix, matrix, 2, 5, -1, -1);
	printMat(newmatrix,2,5);
}
