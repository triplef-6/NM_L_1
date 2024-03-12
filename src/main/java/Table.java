import java.util.stream.IntStream;

public class Table {
    public static void printTable(int n/*, double alpha, double beta*/) {
        System.out.println(
                "+-----+-----+--------+---------+-------+--------+-----+--------+\n" +
                "!  α  !  β  ! ||A||∞ ! ||A_||∞ ! ν∞(A) ! ||z||∞ !  ζ  ! ||R||∞ !\n" +
                "+-----+-----+--------+---------+-------+--------+-----+--------+");
        double alpha = 1;
        double beta = 10;
        int k = 2;
        do {
            // ||A||∞  ||A_||∞  ν∞(A)
            double[][] A = new double[n][n];
            double[][] A_gen = new double[n][n];
            Gen gen = new Gen();
            GenMatrixParam genMatrixParam = gen.mygen(A, A_gen, n, alpha, beta, 0, 0, 2, 1);

            // ||z||∞ ||R||∞
            double[][] AN = new double[n][n];
            IntStream.range(0, n).forEach(i -> System.arraycopy(A[i], 0, AN[i], 0, n));
            double[][] A_ = MatrixS.gauss(AN, n);

            double z_norm = z_norm(A_, A_gen, n);
            double r_norm = r_norm(A, A_, n);

            System.out.printf("|%5.0e|%5.0e|%8.2e|%9.3e|%7.1e|%8.2e|%5.0e|%7.2e|", alpha, beta,
                    genMatrixParam.getA_norm(), genMatrixParam.getA_inv_norm(), genMatrixParam.getObusl(),
                    z_norm, z_norm / genMatrixParam.getA_inv_norm(), r_norm
            );
            System.out.println("\n+-----+-----+--------+---------+-------+--------+-----+--------+");
            beta = Math.pow(10, k);
            k++;
        } while (beta < 1.0E20); // 10 ** 12






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
        double[][] A_A = new double[n][n];
        gen.matr_mul(A, A_, A_A, n);

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






}
