import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

import static java.lang.Math.*;

public class Main {
    public static void main(String[] args) throws IOException {

        int n;
        double[][] matrixCreated;
        double[] matrixB;

        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));

        System.out.println("Enter the number of variables in the equation:");
        n = Integer.parseInt(reader.readLine());
        matrixCreated= new double[n][n];
        matrixB=new double[n];
        System.out.println("Enter the matrix:");

        for (int i = 0; i < n; i++)
        {
            StringTokenizer matrixToken = new StringTokenizer(reader.readLine());

            while (matrixToken.hasMoreTokens())
                for (int j = 0; j < n && matrixToken.hasMoreTokens(); j++)
                    matrixCreated[i][j] = Integer.parseInt(matrixToken.nextToken());
        }

        System.out.println("enter the right side of the equations");
        StringTokenizer results=new StringTokenizer(reader.readLine());
        while(results.hasMoreTokens()){
            for(int k=0;k<n;k++)
                matrixB[k]=Integer.parseInt(results.nextToken());
        }

        double[] x = lsolve(matrixCreated, matrixB);//solve the matrix given


        // print results
        for (int i = 0; i < n; i++) {
            System.out.println(x[i]);
        }

    }

    // Gaussian elimination with partial pivoting
    public static double[] lsolve(double[][] mainMatrix, double[] rightSide) {
        int n = rightSide.length;

        for (int p = 0; p < n; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (abs(mainMatrix[i][p]) > abs(mainMatrix[max][p])) {//find largest absolute value for pivot
                    max = i;
                }
            }
            double[] temp = mainMatrix[p]; mainMatrix[p] = mainMatrix[max]; mainMatrix[max] = temp;
            double   t;
            t = rightSide[p];
            rightSide[p] = rightSide[max]; rightSide[max] = t;

            // pivot within A
            for (int i = p + 1; i < n; i++) {
                double alpha = mainMatrix[i][p] / mainMatrix[p][p];
                rightSide[i] -= alpha * rightSide[p];
                for (int j = p; j < n; j++) {
                    mainMatrix[i][j] -= alpha * mainMatrix[p][j];
                }
            }
        }

        // back substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += mainMatrix[i][j] * x[j];
            }
            x[i] = (rightSide[i] - sum) / mainMatrix[i][i];
        }
        return x;
    }

}