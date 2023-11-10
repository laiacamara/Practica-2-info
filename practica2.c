#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 512

float Mat[N][N];
float MatDD[N][N];

float V1[N];
float V2[N];
float V3[N];
float V4[N];

void initData(){
int i, j;
srand(4422543);
for( i = 0; i < N; i++ )
	for( j = 0; j < N; j++ ){
		Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
		if ( (abs(i - j) <= 3) && (i != j))
 			MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
 		else if ( i == j )
 			MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
 		else MatDD[i][j] = 0.0;
 	}
for( i = 0; i < N; i++ ){
	V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 	V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
	V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
	}
}

/*1*/

void PrintVect( float vect[N], int from, int numel) {
	for (int i = from; i < from + numel; i++) {
                printf("%.6f  ",vect[i]);
                }
                printf("\n");
}

/*2*/

void PrintRow( float Mat[N][N], int row, int from, int numel ) {
	for (int i = from; i < from + 10 && i < from + row + numel; i++) {
		printf("%.6f  ",Mat[row][i]);
		}
		printf("\n");
}

/*3*/

void MultEscalar(float vect[N], float vectres[N], float alfa) {
	for (int i = 0; i < N; i++) {
		vectres[i] = vect[i] * alfa;
                }
}

/*4*/

float Scalar(float vect1[N], float vect2[N]) {
	float Resultat = 0.0;
	for (int i = 0; i < N; i++) {
        	Resultat += (vect1[i] * vect2[i]);
	}
    return Resultat;
}

/*5*/

float Magnitude(float vect1[N]) {
  float magnitude = 0.0;
 for (int i = 0; i < N; i++) {
    magnitude += pow(vect1[i], 2);
  }
  magnitude = sqrt(magnitude);
  return magnitude;
}
/*6*/

int Ortogonal(float vect1[N], float vect2[N]) {
	float producteEscalar = 0.0;
	for (int i = 0; i < N; i++) {
        	producteEscalar += (vect1[i] * vect2[i]);
	}
	if (producteEscalar == 0) {
        	return 1;
	} else {
		return 0;
	}
}

/*7*/

void Projection(float vect1[N], float vect2[N], float vectres[N]) {
	float scalar = Scalar(vect1, vect2) / Magnitude(vect2);
	MultEscalar(vect2, vectres, scalar);
}

/*8*/

float Infininorm(float M[N][N]) {
	float maxNorm = 0;
	for (int i = 0; i < N; i++) {
        	float rowSum = 0;
        	for (int j = 0; j < N; j++) {
            		rowSum += fabs(M[i][j]);
        	}
        if (rowSum > maxNorm) {
            maxNorm = rowSum;
	}
	}
	return maxNorm;
}

/*9*/

float Onenorm(float M[N][N]) {
    float norm = 0.0;

    for (int j = 0; j < N; j++) {
        float sum = 0.0;
        for (int i = 0; i < N; i++) {
            sum += fabs(M[i][j]);
        }
        if (sum > norm) {
            norm = sum;
        }
    }

    return norm;
}


/*10*/

float NormFrobenius(float M[N][N]) {
	float frobeniusNorm = 0.0;
	for (int i = 0; i < N; i++) {
        	for (int j = 0; j < N; j++) {
			frobeniusNorm += M[i][j] * M[i][j];
		}
	}

        frobeniusNorm = sqrt(fabs(frobeniusNorm));
	return frobeniusNorm;

}
/*11*/

int DiagonalDom(float M[N][N]) {
    for (int i = 0; i < N; i++) {
        float diag_value = fabs(M[i][i]);
        float row_sum = 0.0;

        for (int j = 0; j < N; j++) {
            row_sum += fabs(M[i][j]);
        }

        row_sum -= diag_value;

        if (diag_value <= 2 * row_sum) {
            continue;
        } else {
            return 0;
        }
    }

    return 1;
}


/*12*/


int Jacobi( float M[N][N] , float vect[N], float vectres[N], unsigned int iter ){
	if (!DiagonalDom(M)) {
		printf("La matriu M no és diagonal dominant, no es pot aplicar Jacobi\n");
	return 0;
        }
	float temp[N];
	unsigned int k;
	for (k=0; k< iter;k++){
		for (int i = 0; i < N; i++) {
        	temp[i]=vect[i];

        for (int j=0; j<N; j++){
        	if (j!=i){
        		temp[i]-=M[i][j]*vectres[j];
        	}
        }
        temp[i] /=M[i][i];
        }
      	for (int i=0; i<N; i++){
        	vectres[i]=temp[i];
      	}
    	}
    return 1;
}


/*main*/

int main() {
    initData();

/*a*/
	printf("V1 del 0 al 9 i del 256 al 265 \n");
	PrintVect( V1, 0, 10 );
        printf("V2 del 0 al 9 i del 256 al 265 \n");
        PrintVect( V2, 0, 10 );
        printf("V3 del 0 al 9 i del 256 al 265 \n");
        PrintVect( V3, 0, 10 );
	printf("\n");
/*b*/
        printf("Mat fila 0 i fila 100 del 0 al 9 \n");
	PrintRow( Mat, 0, 0, 9);
        PrintRow( Mat, 100, 0, 9);
	printf("\n");

/*c*/
        printf("MatDD fila 0 del 0 al 9 i fila 100 del 95 al 104 \n");
        PrintRow( MatDD, 0, 0, 9);
        PrintRow( MatDD, 100, 90,99);
        printf("\n");
/*d*/
        float normaInfini = Infininorm(Mat);
        printf("Infininorma de Mat: %f\n", normaInfini);
	printf("Norma ú de Mat = %f\n", Onenorm(Mat));
        float M[N][N];
        float normaFrobenius = NormFrobenius(Mat);
        printf("Norma de Frobenius de Mat: %f\n", normaFrobenius);
	if (DiagonalDom(Mat) == 1) {
		printf("La matriu Mat no és diagonal dominant\n");
	} else {
		printf("La matriu Mat és diagonal dominant\n");
	}
        printf("\n");


	float normaInfini2 = Infininorm(MatDD);
        printf("Infiniorma de MatDD: %f\n", normaInfini2);
	printf("Norma ú de MatDD = %f\n", Onenorm(MatDD));
        float normaFrobenius2 = NormFrobenius(MatDD);
        printf("Norma de Frobenius de MattDD: %f\n", normaFrobenius2);
                if (DiagonalDom(MatDD) == 1) {
                printf("La matriu Mat no és diagonal dominant\n");
        } else {
                printf("La matriu Mat és diagonal dominant\n");
        }
	printf("\n");

/*e*/
        printf("\n");
        float Resultat1 = Scalar(V1,V2);
        printf("Escalar V1,V2 = %.6f\n", Resultat1);
        float Resultat2 = Scalar(V3, V1);
        printf("Escalar V1,V3 = %.6f\n", Resultat2);
        float Resultat3 = Scalar(V3, V2);
        printf("Escalar V2,V3 = %.6f\n", Resultat3);
        printf("\n");
/*f*/

        float mag1 = Magnitude(V1);
        float mag2 = Magnitude(V2);
        float mag3 = Magnitude(V3);
        printf("Magnitud V1 = %f\n", mag1);
        printf("Magnitud V2 = %f\n", mag2);
        printf("Magnitud V3 = %f\n", mag3);
        printf("\n");

/*g*/

        if (Ortogonal(V1, V2)) {
                printf("V1 i V2 són ortogonals.\n");
        } else {
                printf("V1 i V2 no són ortogonals.\n");
        }

        if (Ortogonal(V1, V3)) {
                printf("V1 i V3 són ortogonals.\n");
        } else {
                printf("V1 i V3 no són ortogonals.\n");
        }

        if (Ortogonal(V2, V3)) {
                printf("V2 i V3 són ortogonals.\n");
        } else {
                printf("V2 i V3 no són ortogonals.\n");
        }
        printf("\n");

/*h*/
        printf("Els elements 0 al 9 i 256 al 265 del resultat de multiplicar V3x2.0 són: \n");
	float alfa = 2.0;
	float V3Resultat[N];
        MultEscalar( V3, V3Resultat, alfa);
	for (int i = 0; i < 10; i++) {
		printf("%.6f ", V3Resultat[i]);
	}
	for (int i = 256; i < 266; i++) {
    		printf("%.6f ", V3Resultat[i]);
	}
	printf("\n");

/*i*/
	printf("\n");
	printf("Els elements 0 a 9 del resultat de la projecció de V2 sobre V3 són: ");
	Projection(V2, V3, V4);
	PrintVect(V4, 0, 10);
	printf("Els elements 0 a 9 del resultat de la projecció de V1 sobre V2 són:\n");
	Projection(V1, V2, V4);
	PrintVect(V4, 0, 10);
	printf("\n");

/*j*/

	float Mat[N][N], Vect[N], Result[N];
	float solJac1[N];
	float solJac2[N];
	Jacobi(MatDD, V3, solJac1, 1);
	Jacobi(MatDD, V3, solJac2, 1000);
	printf("Els elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
	PrintVect(solJac1, 0, 10);
	printf("Els elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
	PrintVect(solJac2, 0, 10);
	if (DiagonalDom(Mat) == 0) {
		printf("La matriu M no és diagonal dominant, no es pot aplicar Jacobi\n");
	} else {
		printf("La matriu M és diagonal dominant\n");
	}

	return 0;

}

