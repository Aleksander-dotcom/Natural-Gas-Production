package zadania.EGZ;

public class Well {
    private int H, re;
    private double rw, k, D, Sm, h, beta;

    public Well(int H, int re, double h, double Sm, double beta, double rw, double k, double D) {
        this.H = H;
        this.re = re;
        this.h = h;
        this.Sm = Sm;
        this.beta = beta;
        this.rw = rw;
        this.k = k;
        this.D = D;
    }

    public int getH() {return H;}
    public int getRe() {return re;}
    public double geth() {return h;}
    public double getSm() {return Sm;}
    public double getBeta() {return beta;}
    public double getRw() {return rw;}
    public double getK() {return k;}
    public double getD() {return D;}
}