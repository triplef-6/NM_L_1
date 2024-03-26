import java.util.stream.IntStream;

public class Table {
    public static void printTableBeta(int n, int maxSBeta) {
        System.out.println("""
                        +-----+-----+----------+-----------+---------+----------+-------+----------+
                        !  α  !  β  !  ||A||∞  !  ||A_||∞  !  ν∞(A)  !  ||z||∞  !   ζ   !  ||R||∞  !
                        +-----+-----+----------+-----------+---------+----------+-------+----------+""");
        double alpha = 1;
        double beta = 10;
        double maxBeta = Math.pow(10, maxSBeta);
        int k = 2;
        do {
            // ||A||∞  ||A_||∞  ν∞(A)
            double[][] A = new double[n][n];
            double[][] A_gen = new double[n][n];
            Gen gen = new Gen();
            GenMatrixParam genMatrixParam = gen.mygen(A, A_gen, n, alpha, beta, 0, 0, 2, 1, false);

            // ||z||∞ ||R||∞
            double[][] AN = new double[n][n];
            IntStream.range(0, n).forEach(i -> System.arraycopy(A[i], 0, AN[i], 0, n));
            double[][] A_ = gauss(AN, n);

            double z_norm = z_norm(A_, A_gen, n);
            double r_norm = r_norm(A, A_, n);

            System.out.printf("|%5.0e|%5.0e|%10.4e|%11.5e|%9.3e|%10.4e|%7.1e|%9.4e|", alpha, beta,
                    genMatrixParam.getA_norm(), genMatrixParam.getA_inv_norm(), genMatrixParam.getObusl(),
                    z_norm, z_norm / genMatrixParam.getA_inv_norm(), r_norm
            );
            System.out.println("\n+-----+-----+----------+-----------+---------+----------+-------+----------+");

            beta = Math.pow(10, k);
            k++;
        } while (beta < maxBeta); // 10 ** 12
    }

    public static void printTableAlpha(int n, int maxSAlpha) {
        System.out.println("""
                        +-----+-----+----------+-----------+---------+----------+-------+----------+
                        !  α  !  β  !  ||A||∞  !  ||A_||∞  !  ν∞(A)  !  ||z||∞  !   ζ   !  ||R||∞  !
                        +-----+-----+----------+-----------+---------+----------+-------+----------+""");
        double alpha = 0.1;
        double maxAlpha = Math.pow(10, -maxSAlpha);
        double beta = 1;
        int k = 2;
        do {
            // ||A||∞  ||A_||∞  ν∞(A)
            double[][] A = new double[n][n];
            double[][] A_gen = new double[n][n];
            Gen gen = new Gen();
            GenMatrixParam genMatrixParam = gen.mygen(A, A_gen, n, alpha, beta, 0, 0, 2, 1, false);

            // ||z||∞ ||R||∞
            double[][] AN = new double[n][n];
            IntStream.range(0, n).forEach(i -> System.arraycopy(A[i], 0, AN[i], 0, n));
            double[][] A_ = gauss(AN, n);

            double z_norm = z_norm(A_, A_gen, n);
            double r_norm = r_norm(A, A_, n);

            System.out.printf("|%5.0e|%5.0e|%10.4e|%11.5e|%9.3e|%10.4e|%7.1e|%9.4e|", alpha, beta,
                    genMatrixParam.getA_norm(), genMatrixParam.getA_inv_norm(), genMatrixParam.getObusl(),
                    z_norm, z_norm / genMatrixParam.getA_inv_norm(), r_norm
            );
            System.out.println("\n+-----+-----+----------+-----------+---------+----------+-------+----------+");

            alpha = Math.pow(10, -k);
            k++;
        } while (alpha > maxAlpha); //
    }

    public static double[][] gauss(double[][] A, int n) {
        Gen g = new Gen();
        double[][] E = genE(n);

        // приводим к верхнетреугольному виду
        for (int k = 0; k < n - 1; k++) {
            double max_a = 0;
            int max_i = k;
            for (int i = k; i < n; i++) {
                if (Math.abs(A[i][k]) > max_a) {
                    max_a = A[i][k];
                    max_i = i;
                }
            }

            if (max_i != k) {
                swapLines(A, n, k, max_i);
                swapLines(E, n, k, max_i);
            }

            for (int i = k + 1; i < n; i++) {
                double t1 = A[k][k];
                double t2 = A[i][k];
                for (int j = k; j < n; j++) {
                    A[k][j] /= t1;
                    E[k][j] /= t1;

                    A[i][j] -= A[k][j] * t2;
                    E[i][j] -= E[k][j] * t2;
                }
            }
        }

        double t = A[n - 1][n - 1];
        A[n - 1][n - 1] /= t;
        for (int j = 0; j < n; j++) {
            E[n - 1][j] /= t;
        }

        // приводим к диаганальному виду
        for (int k = n - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                double t2 = A[i][k];

                for (int j = 0; j < k; j++) {
                    E[i][j] = E[i][j] - E[k][j] * t2;
                }
                for (int j = k; j < n; j++) {
                    A[i][j] = A[i][j] - A[k][j] * t2;
                    E[i][j] = E[i][j] - E[k][j] * t2;
                }
            }
        }

        return E;
    }

    private static double z_norm(double[][] A_, double[][] A_gen, int n) {
        double[][] Z = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Z[i][j] = A_[i][j] - A_gen[i][j];
            }
        }
        Gen gen = new Gen();
        return gen.matr_inf_norm(Z, n);
    }

    private static double r_norm(double[][] A, double[][] A_, int n) {
        Gen gen = new Gen();
        double[][] A_A = gen.matr_mul(A, A_, n);

        double[][] R = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    R[i][j] = A_A[i][j] - 1.;
                } else {
                    R[i][j] = A_A[i][j] - 0.;
                }
            }
        }
        return gen.matr_inf_norm(R, n);
    }

    private static void swapLines(double[][] A, int n, int I, int J) {
        for (int k = 0; k < n; k++) {
            double tmp = A[I][k];
            A[I][k] = A[J][k];
            A[J][k] = tmp;
        }
    }

    private static double[][] genE(int n) {
        double[][] E = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    E[i][j] = 1;
                } else {
                    E[i][j] = 0;
                }
            }
        }
        return E;
    }
}

