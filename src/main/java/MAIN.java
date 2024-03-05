public class MAIN {
    public static void main(String[] args) {
        int n = 5;
        double[][] A = new double[n][n];
        double[][] A_ = new double[n][n];

        Gen gen = new Gen();

        gen.mygen(A, A_, n, 1, 5, 0, 0, 2, 1);

        gen.print_matr(A_, n);
        gen.print_matr(MatrixS.gauss(A, n), n);
        double[][] A_A = new double[n][n];
        gen.matr_mul(A, A_, A_A, n);
        gen.print_matr(A_A, n);





    }
}
