public class GenMatrixParam {
    private double A_norm; // ||  A  || ||A||∞
    private double A_inv_norm; // ||A_inv|| ||A_||∞
    private double obusl; // obusl ν∞(A)

    public GenMatrixParam() {
    }

    public double getA_norm() {
        return A_norm;
    }

    public void setA_norm(double a_norm) {
        A_norm = a_norm;
    }

    public double getA_inv_norm() {
        return A_inv_norm;
    }

    public void setA_inv_norm(double a_inv_norm) {
        A_inv_norm = a_inv_norm;
    }

    public double getObusl() {
        return obusl;
    }

    public void setObusl(double obusl) {
        this.obusl = obusl;
    }
}
