package zadania.EGZ;

import java.math.BigDecimal;
import java.math.RoundingMode;

public class Functions {

    Reservoir reservoir = new Reservoir(22.5E6, 10_800_000, 44.3, 0.01, 0.0025, 300, 273.15, 0.2, 0.16);

    Gas gas = new Gas(new double[]{0.7, 0.03, 0.02, 0.25, 0.0});

    Well well1 = new Well(1460, 500, 13.2, 1.2, 1E8, 0.15, 56E-15, 0.0603);
    Well well2 = new Well(1310, 550, 13.4, 1.3, 1E8, 0.12, 64E-15, 0.0730);
    Well well3 = new Well(1360, 550, 13.4, 1.4, 1E8, 0.12, 86E-15, 0.0762);
    Well well4 = new Well(1390, 550, 13.8, 1.3, 1E8, 0.12, 72E-15, 0.0730);
    Well well5 = new Well(1430, 600, 13.3, 1.6, 1E8, 0.12, 72E-15, 0.0762);
    Well well6 = new Well(1360, 700, 14.1, 1.1, 1E8, 0.15, 56E-15, 0.0730);
    Well well7 = new Well(1500, 600, 14.1, 1.0, 1E8, 0.15, 72E-15, 0.0762);
    Well well8 = new Well(1490, 450, 14.2, 0, 1E8, 0.175, 54E-15, 0.0730);
    Well well9 = new Well(1320, 650, 14.0, -1.0, 1E8, 0.175, 56E-15, 0.0762);
    Well well10 = new Well(1340, 570, 13.7, 1.0, 1E8, 0.175, 86E-15, 0.0730);
    Well well11 = new Well(1420, 520, 13.9, -1.0, 1E8, 0.175, 56E-15, 0.0762);
    Well well12 = new Well(1470, 510, 11, 0, 1E8, 0.15, 42E-15, 0.0762);
    Well well13 = new Well(1490, 560, 13.8, 0, 1E8, 0.12, 72E-15, 0.0762);
    Well well14 = new Well(1370, 550, 13.6, 1.5, 1E8, 0.15, 56E-15, 0.0603);
    Well well15 = new Well(1390, 500, 13.4, 1.0, 1E8, 0.12, 86E-15, 0.0730);
    Well well16 = new Well(1410, 610, 13.5, 1.0, 1E8, 0.15, 64E-15, 0.0730);
    Well well17 = new Well(1490, 630, 13.6, 1, 1E8, 0.12, 86E-15, 0.0603);
    Well well18 = new Well(1370, 580, 13.9, 0, 1E8, 0.12, 54E-15, 0.0730);

    Well[] wells = new Well[]{well1, well2, well3, well4, well5, well6, well7, well8, well9, well10, well11, well12, well13, well14, well15, well16, well17, well18};


    public double round(double value, int precision) {
        BigDecimal bd = new BigDecimal(value).setScale(precision, RoundingMode.HALF_UP);

        return bd.doubleValue();
    }

    public String format(double value, int precision, int w) {
        if (w == 0){
            return String.format("%6.2E", value);
        }
        return String.format("%,." + precision + "f", value);
    }

    public double pseudoz(double p) {
        double Tpr, m, alfa, aa, bb, AA, BB, ZZ, b2, b1, b0, pp, qq, dd;
        Tpr = reservoir.getTzl() / gas.getTc();
        m = 0.48508 + 1.55171 * gas.getW() - 0.15613 * Math.pow(gas.getW(), 2);
        alfa = Math.pow(1 + m * (1 - Math.sqrt(Tpr)), 2);
        aa = 0.42747 * Math.pow(gas.getRu(), 2) * Math.pow(gas.getTc(), 2) / gas.getPc();
        bb = 0.08664 * gas.getRu() * gas.getTc() / gas.getPc();
        AA = aa * alfa * p / Math.pow(gas.getRu(), 2) / Math.pow(reservoir.getTzl(), 2);
        BB = bb * p / gas.getRu() / reservoir.getTzl();
        b2 = -1;
        b1 = AA - BB - Math.pow(BB, 2);
        b0 = -AA * BB;
        pp = b1 - (Math.pow(b2, 2)) / 3;
        qq = b0 - (b1 * b2) / 3 + (2 * Math.pow(b2, 3)) / 27;
        dd = Math.pow(pp / 3, 3) + Math.pow(qq / 2, 2);
        if (dd > 0) {
            ZZ = Math.pow((-qq / 2 + Math.sqrt(dd)), 1.0 / 3.0) - b2 / 3 + Math.pow((-qq / 2 - Math.sqrt(dd)), 1.0 / 3.0);
        } else {
            return -1;
        }
        return ZZ;
    }

    public double density(double p) {
        return p / (pseudoz(p) * (gas.getRu() / gas.getM()) * reservoir.getTzl());
    }

    public double density(double p, double T) {
        return p / (pseudoz(p) * (gas.getRu() / gas.getM()) * T);
    }

    public double viscosity(double p) {
        double x = 3.5 + 986 / (1.8 * reservoir.getTzl()) + 0.01 * gas.getM();
        double y = 2.4 - 0.2 * x;
        double u = (9.4 + 0.02 * gas.getM()) * (Math.pow(1.8 * reservoir.getTzl(), 1.5)) / (209 + 1.8 * reservoir.getTzl() + 19 * gas.getM());
        return u * (Math.exp(x * (Math.pow(density(p) / 1000, y)))) * 1E-7;

    }

    public double[] flow(double Pzl, double Pdd) {
        double[] result = new double[wells.length];
        for (int i = 0; i < wells.length; i++) {
            double a = (viscosity(Pzl) * gas.getPn() * pseudoz(Pzl) * reservoir.getTzl() * ((Math.log10(wells[i].getRe() / wells[i].getRw())) - 3 / 4.0 + wells[i].getSm())) /
                    (Math.PI * wells[i].getK() * wells[i].geth() * gas.getTn());
            double b = ((viscosity(Pzl) * gas.getPn() * pseudoz(Pzl) * reservoir.getTzl()) / (Math.PI * wells[i].getK() * wells[i].geth() * gas.getTn())) *
                    (wells[i].getBeta() * density(gas.getPn(), gas.getTn()) * wells[i].getK()) / (2 * Math.PI * wells[i].geth() * wells[i].getRw() * viscosity(Pdd));
            double q = (-a + Math.sqrt(Math.pow(a, 2) + 4 * b * (Math.pow(Pzl, 2) - Math.pow(Pdd, 2)))) / (2 * b);
            result[i] = round(q, 3);
        }
        return result;
    }

    public double[] conditions(double Pzl, double Pdd, int t) {
        int i = 0;
        double[] result = new double[wells.length];
        while (true) {
            if (t == 0 && i < 18) {
                result[i] = 0;
                i++;
            } else if (t == 1 && i < 3) {
                result[i] = flow(Pzl, Pdd)[i];
                i++;
            } else if (t == 2 && i < 6) {
                result[i] = flow(Pzl, Pdd)[i];
                i++;
            } else if (t == 3 && i < 9) {
                result[i] = flow(Pzl, Pdd)[i];
                i++;
            } else if (t == 4 && i < 12) {
                result[i] = flow(Pzl, Pdd)[i];
                i++;
            } else if (t == 5 && i < 15) {
                result[i] = flow(Pzl, Pdd)[i];
                i++;
            } else if (t >= 6 && i < 18) {
                result[i] = flow(Pzl, Pdd)[i];
                i++;
            } else {
                break;
            }
        }
        return result;
    }

    public double[] gasProd(double[] q) {
        double[] result = new double[wells.length];
        for (int i = 0; i < wells.length; i++) {
            double Gp = q[i] * 365 * 86400;
            result[i] = Gp;
        }
        return result;
    }

    public double cumGasProduction(double[] gasProd, double cumGp) {
        double cumGpAnnual = 0;
        for (double singleGp : gasProd) {
            cumGpAnnual += singleGp;
        }
        cumGp += cumGpAnnual;
        return cumGp;
    }

    public double[] lambda(double Pzl, double Pdd) {
        double[] result = new double[wells.length];
        double E = 0.01;
        for (int i = 0; i < wells.length; i++) {
            double q = flow(Pzl, Pdd)[i];
            double v = q / (Math.PI * Math.pow(wells[i].getD(), 2) / 4);
            double Re = wells[i].getD() * v * density(Pdd) * wells[i].getH() / viscosity(Pdd);
            double f = Math.pow((1 / (-2 * Math.log10((E / 1000 / wells[i].getD() / 3.71 + 15 / Re)))), 2);
            result[i] = f;
        }
        return result;
    }

    public double[] wellheadPressure(double[] q, double Pzl, double Pdd) {
        double[] result = new double[wells.length];
        for (int i = 0; i < wells.length; i++) {
            double exp = Math.exp((2 * gas.getG() * wells[i].getH()) / (pseudoz(Pdd) * (gas.getRu() / gas.getM()) * (reservoir.getTzl() + reservoir.getTgl()) / 2));
            double B = 8 * lambda(Pzl, Pdd)[i] * Math.pow(gas.getPn(), 2) * Math.pow(pseudoz(Pdd), 2) * Math.pow((reservoir.getTzl() + reservoir.getTgl()) / 2, 2) * (exp - 1) /
                    (Math.pow(Math.PI, 2) * Math.pow(wells[i].getD(), 5) * Math.pow(gas.getTn(), 2) * gas.getG());
            double Pgd = Math.sqrt((Math.pow(Pdd, 2) - B * Math.pow(q[i], 2)) / exp);
            result[i] = Pgd;
        }
        return result;
    }

    public double wellboreFlowPressure(double Pzl, int t) {
        double Pdd;
        if (t <= 1) {
            Pdd = Pzl * (1 - reservoir.getS());
        } else {
            Pdd = Pzl * (1 - reservoir.getS() - reservoir.getdS() * (t - 1));
        }
        return Pdd;
    }

    public double Gi() {
        double Bgi = pseudoz(reservoir.getPi()) * gas.getPn() * reservoir.getTzl() / reservoir.getPi() / gas.getTn();
        return reservoir.getAi() * reservoir.getHsr() * reservoir.getFi() * (1 - reservoir.getSw()) / Bgi;
    }

    public double newPressure(double cumGp) {
        double zz = 1, dz, pNew;
        do {
            pNew = zz * ((1 - cumGp / Gi()) * (reservoir.getPi() / pseudoz(reservoir.getPi())));
            double Znew = pseudoz(pNew);
            dz = Math.abs(Znew - zz);
            zz = Znew;
        }
        while (dz > 0.001);
        return pNew;
    }

    public void initialization() {
        double Pzl = reservoir.getPi();
        double cumGp = 0;
        double zPzl, Pdd, zPdd, densityPzl, densityPdd, viscosityPzl, viscosityPdd;
        double[] flows, gasProd, Pgd;
        System.out.println("Time [yr]" + "     " + "Pressure [Pa]" + "      " + "Z factor (P)" + "       " + "Flowing bottom-hole pressure [Pa]"
                            + "       " + "Z factor (Pwf)" + "       " + "Density (P) [kg/m3]" + "       " + "Density (Pwf) [kg/m3]"
                            + "       " + "Viscosity (P) [Pas]" + "       " + "Viscosity (Pwf) [Pas]" + "       " + "Flow (1 well) [m3/s]"
                            + "       " + "Flow (2 well) [m3/s]" + "       " + "Flow (3 well) [m3/s]" + "       " + "Flow (4 well) [m3/s]"
                            + "       " + "Flow (5 well) [m3/s]" + "       " + "Flow (6 well) [m3/s]" + "       " + "Flow (7 well) [m3/s]"
                            + "       " + "Flow (8 well) [m3/s]" + "       " + "Flow (9 well) [m3/s]" + "       " + "Flow (10 well) [m3/s]"
                            + "       " + "Flow (11 well) [m3/s]" + "       " + "Flow (12 well) [m3/s]" + "       " + "Flow (13 well) [m3/s]"
                            + "       " + "Flow (14 well) [m3/s]" + "       " + "Flow (15 well) [m3/s]" + "       " + "Flow (16 well) [m3/s]"
                            + "       " + "Flow (17 well) [m3/s]" + "       " + "Flow (18 well) [m3/s]" + "       " + "Gp (1 well) [m3]"
                            + "       " + "Gp (2 well) [m3]" + "       " + "Gp (3 well) [m3]" + "       " + "Gp (4 well) [m3]"
                            + "       " + "Gp (5 well) [m3]" + "       " + "Gp (6 well) [m3]" + "       " + "Gp (7 well) [m3]"
                            + "       " + "Gp (8 well) [m3]" + "       " + "Gp (9 well) [m3]" + "       " + "Gp (10 well) [m3]"
                            + "       " + "Gp (11 well) [m3]" + "       " + "Gp (12 well) [m3]" + "       " + "Gp (13 well) [m3]"
                            + "       " + "Gp (14 well) [m3]" + "       " + "Gp (15 well) [m3]" + "       " + "Gp (16 well) [m3]"
                            + "       " + "Gp (17 well) [m3]" + "       " + "Gp (18 well) [m3]" + "       " + "CUM Gp [m3]"
                            + "          " + "Pwh (1 well) [Pa]" + "          " + "Pwh (2 well) [Pa]" + "          " + "Pwh (3 well) [Pa]"
                            + "          " + "Pwh (4 well) [Pa]" + "          " + "Pwh (5 well) [Pa]" + "          " + "Pwh (6 well) [Pa]"
                            + "          " + "Pwh (7 well) [Pa]" + "          " + "Pwh (8 well) [Pa]" + "          " + "Pwh (9 well) [Pa]"
                            + "          " + "Pwh (10 well) [Pa]" + "          " + "Pwh (11 well) [Pa]" + "          " + "Pwh (12 well) [Pa]"
                            + "          " + "Pwh (13 well) [Pa]" + "          " + "Pwh (14 well) [Pa]" + "          " + "Pwh (15 well) [Pa]"
                            + "          " + "Pwh (16 well) [Pa]" + "          " + "Pwh (17 well) [Pa]" + "          " + "Pwh (18 well) [Pa]");
        for (int t = 0; t < 50 + 1; t++) {
            zPzl = round(pseudoz(Pzl), 3);
            Pdd = round(wellboreFlowPressure(Pzl, t), 3);
            zPdd = round(pseudoz(Pdd), 3);
            densityPzl = round(density(Pzl), 3);
            densityPdd = round(density(Pdd), 3);
            viscosityPzl = round(viscosity(Pzl), 8);
            viscosityPdd = round(viscosity(Pdd), 8);
            flows = conditions(Pzl, Pdd, t);
            gasProd = gasProd(flows);
            Pgd = wellheadPressure(flows, Pzl, Pdd);
            cumGp = round(cumGasProduction(gasProd, cumGp), 3);

            System.out.println("   " + t
                    + "         " + format(Pzl, 3,1)
                    + "         " + format(zPzl, 3,1)
                    + "                    " + format(Pdd, 3,1)
                    + "                     " + format(zPdd,3,1)
                    + "                  " + format(densityPzl,3,1)
                    + "                     " + format(densityPdd,3,1)
                    + "                    " + format(viscosityPzl,8,0)
                    + "                   " + format(viscosityPdd,8,0)

                    + "                      " + format(flows[0],3,1)       + "                      " + format(flows[1],3,1)
                    + "                      " + format(flows[2],3,1)       + "                      " + format(flows[3],3,1)
                    + "                      " + format(flows[4],3,1)       + "                      " + format(flows[5],3,1)
                    + "                      " + format(flows[6],3,1)       + "                      " + format(flows[7],3,1)
                    + "                      " + format(flows[8],3,1)       + "                      " + format(flows[9],3,1)
                    + "                      " + format(flows[10],3,1)      + "                        " + format(flows[11],3,1)
                    + "                       " + format(flows[12],3,1)     + "                       " + format(flows[13],3,1)
                    + "                       " + format(flows[14],3,1)     + "                       " + format(flows[15],3,1)
                    + "                       " + format(flows[16],3,1)     + "                       " + format(flows[17],3,1)

                    + "                   " + format(gasProd[0],3,0)        + "               " + format(gasProd[1],3,0)
                    + "               " + format(gasProd[2],3,0)            + "               " + format(gasProd[3],3,0)
                    + "               " + format(gasProd[4],3,0)            + "               " + format(gasProd[5],3,0)
                    + "               " + format(gasProd[6],3,0)            + "               " + format(gasProd[7],3,0)
                    + "               " + format(gasProd[8],3,0)            + "               " + format(gasProd[9],3,0)
                    + "               " + format(gasProd[10],3,0)           + "                 " + format(gasProd[11],3,0)
                    + "                " + format(gasProd[12],3,0)          + "                " + format(gasProd[13],3,0)
                    + "                " + format(gasProd[14],3,0)          + "                " + format(gasProd[15],3,0)
                    + "                " + format(gasProd[16],3,0)          + "                " + format(gasProd[17],3,0)
                    + "             " + format(cumGp,3,0)
                    + "                " + format(Pgd[0],3,0)               + "                   " + format(Pgd[1],3,0)
                    + "                   " + format(Pgd[2],3,0)            + "                   " + format(Pgd[3],3,0)
                    + "                   " + format(Pgd[4],3,0)            + "                   " + format(Pgd[5],3,0)
                    + "                   " + format(Pgd[6],3,0)            + "                   " + format(Pgd[7],3,0)
                    + "                   " + format(Pgd[8],3,0)            + "                   " + format(Pgd[9],3,0)
                    + "                   " + format(Pgd[10],3,0)           + "                     " + format(Pgd[11],3,0)
                    + "                    " + format(Pgd[12],3,0)          + "                    " + format(Pgd[13],3,0)
                    + "                    " + format(Pgd[14],3,0)          + "                    " + format(Pgd[15],3,0)
                    + "                    " + format(Pgd[16],3,0)          + "                    " + format(Pgd[17],3,0));

            Pzl = round(newPressure(cumGp), 3);
        }
    }
}