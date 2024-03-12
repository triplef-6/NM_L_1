public class MAIN {
    public static void main(String[] args) {
        int n = 10;
        double[][] A = new double[n][n];
        double[][] A_ = new double[n][n];

        Gen gen = new Gen();
        gen.mygen(A, A_, n, 1, 2, 0, 0, 2, 1);
//
//        gen.print_matr(A_, n);
//        gen.print_matr(MatrixS.gauss(A, n), n);
        Table.printTable(n);





    }
}
