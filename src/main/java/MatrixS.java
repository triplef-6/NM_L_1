import static java.lang.Long.MIN_VALUE;

public class MatrixS {
    public static double[][] gauss(double[][] A, int n) {
        double[][] AE = addE(A, n);

        // приводим к верхнетреугольному виду
        for (int k = 0; k < n - 1; k++) {
            double max_a = MIN_VALUE;
            int max_i = k;
            for (int i = k; i < n; i++) {
                if (AE[i][k] > max_a) {
                    max_a = AE[i][k];
                    max_i = i;
                }
            }
            if (max_i != k) {
                swapLines(AE, n * 2, k, max_i);
            }

            for (int i = k + 1; i < n; i++) {
                double t1 = AE[k][k];
                double t2 = AE[i][k];
                for (int j = k; j < n * 2; j++) {
                    AE[k][j] /= t1;
                    AE[i][j] -= AE[k][j] * t2;
                }
            }
        }
        double t = AE[n - 1][n - 1];
        for (int j = n - 1; j < n * 2; j++) {
            AE[n - 1][j] /= t;
        }

        // приводим к диаганальному виду
        for (int k = n - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                double t2 = AE[i][k];
                for (int j = k; j < n * 2; j++) {
                    AE[i][j] = AE[i][j] - AE[k][j] * t2;
                }
            }
        }

        return removeE(AE, n);
    }

    private static void swapLines(double[][] A, int n, int I, int J) {
        for (int k = 0; k < n; k++) {
            double tmp = A[I][k];
            A[I][k] = A[J][k];
            A[J][k] = tmp;
        }
    }

    private static double[][] addE(double[][] A, int n) {
        double[][] AE = new double[n][n * 2];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n * 2; j++) {
                if (j < n) {
                    AE[i][j] = A[i][j];
                } else if (i == j - n) {
                    AE[i][j] = 1;
                } else {
                    AE[i][j] = 0;
                }
            }
        }
        return AE;
    }

    private static double[][] removeE(double[][] AE, int n) {
        double[][] A_ = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(AE[i], n, A_[i], 0, n);
        }
        return A_;
    }
}
