package zadania.EGZ;

public class Gas {
    public final double g = 9.80665;
    public final double Tn = 273.15;
    public final double Ru = 8314.13;
    public final double Pn = 101_325;
    private double M,Tc,Pc,w;

    public Gas(double[] x) {
        double[] molarMass = new double[]{16.043, 30.070, 44.097, 28.014, 44.010};
        double[] criticalTemperature = new double[]{190.564, 305.320, 369.830, 126.200, 304.210};
        double[] criticalPressure = new double[]{4590000, 4850000, 4210000, 3390000, 7390000};
        double[] acentricFactor = new double[]{0.011, 0.098, 0.149, 0.037, 0.224};

        Tc = Mean(criticalTemperature,x);
        M = Mean(molarMass,x);
        Pc = Mean(criticalPressure,x);
        w = Mean(acentricFactor,x);
    }

    public double Mean (double[] parametr, double[] x){
        double mean = 0;
        for (int i = 0; i<x.length; i++){
            double current = parametr[i] * x[i];
            mean += current;
        }
        return mean;
    }

    public double getM() {return M;}
    public double getTc() {return Tc;}
    public double getPc() {return Pc;}
    public double getW() {return w;}
    public double getG() {return g;}
    public double getTn() {return Tn;}
    public double getRu() {return Ru;}
    public double getPn() {return Pn;}
}
