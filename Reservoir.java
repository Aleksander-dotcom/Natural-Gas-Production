package zadania.EGZ;

public class Reservoir {

    private double Ai, Pi, hsr, S, dS, Tzl, Tgl, fi, Sw;

    public Reservoir(double Ai, double Pi, double hsr, double S, double dS, double Tzl, double Tgl, double fi, double Sw) {
        this.Ai = Ai;
        this.hsr = hsr;
        this.S = S;
        this.dS = dS;
        this.Tzl = Tzl;
        this.Tgl = Tgl;
        this.fi = fi;
        this.Sw = Sw;
        this.Pi = Pi;
    }

    public double getAi() {return Ai;}
    public double getHsr() {return hsr;}
    public double getS() {return S;}
    public double getdS() {return dS;}
    public double getTzl() {return Tzl;}
    public double getTgl() {return Tgl;}
    public double getFi() {return fi;}
    public double getSw() {return Sw;}
    public double getPi() {return Pi;}
}
