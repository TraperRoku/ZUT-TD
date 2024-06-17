import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;

public class Main {

    private static StringBuilder zmianaNaBity(String string){
        StringBuilder stringBuilder = new StringBuilder();
        for(char c : string.toCharArray()){
            String bit = Integer.toBinaryString((int)c);
            stringBuilder.append(bit);
        }
        return stringBuilder;
    }


    // Kluczowanie z przesuwem amplitudy ASK
    private static double[][] za(int dlugoscBitu, int[] bits, double fn, double A1, double A2, double[] t,double Tb) {
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;
        for (int i = 0; i < dlugoscBitu; i++) {
            if (bits[i] == 0) {
                for (int j = temp; j < temp + (Tb * t.length); j++) {
                    sygnal[i][j] = A1 * Math.sin(2 * Math.PI * fn * t[j - temp]);
                }
            } else {
                for (int j = temp; j < temp + (Tb * t.length); j++) {
                    sygnal[i][j] = A2 * Math.sin(2 * Math.PI * fn * t[j - temp]);
                }
            }
            temp += (Tb * t.length);
        }
        return sygnal;
    }


    //Kluczowanie z przesuwem fazy (PSK):
    private static double[][] zp(int dlugoscBitu, int[] bits, double fn, double[] t,double Tb){
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;
        for (int i = 0; i < dlugoscBitu; i++) {
            if (bits[i] == 0) {
                for (int j = temp; j < temp + (Tb * t.length); j++) {
                    sygnal[i][j] = Math.sin(2*Math.PI * fn * t[j - temp]);
                }
            } else {
                for (int j = temp; j < temp + (Tb * t.length); j++) {
                    sygnal[i][j] = Math.sin(2*Math.PI * fn * t[j - temp] + Math.PI);
                }
            }
            temp += (Tb * t.length);
        }
        return sygnal;
    }
    //Kluczowanie z przesuwem czestotliwosci (FSK):
    private static double[][] zf(int dlugoscBitu, int[] bits, double fn1, double fn2, double[] t,double Tb){
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;
        for (int i = 0; i < dlugoscBitu; i++) {
            if (bits[i] == 0) {
                for (int j = temp; j < temp + (Tb * t.length); j++) {
                    sygnal[i][j] = Math.sin(2 * Math.PI * fn1 * t[j - temp]);
                }
            } else {
                for (int j = temp; j < temp + (Tb * t.length); j++) {
                    sygnal[i][j] = Math.sin(2 * Math.PI * fn2 * t[j - temp] + Math.PI);
                }
            }
            temp += (Tb * t.length);
        }
        return sygnal;
    }

    private static XYSeries createSeries(double[] xData, double[] yData, String name) {
        XYSeries series = new XYSeries(name);
        for (int i = 0; i < xData.length; i++) {
            series.add(xData[i], yData[i]);
        }
        return series;
    }

    private static void displayChart(JFreeChart chart) {
        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Modulated Signals");
            frame.setLayout(new BorderLayout());
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

            ChartPanel chartPanel = new ChartPanel(chart);
            frame.add(chartPanel, BorderLayout.CENTER);

            frame.pack();
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);
        });
    }

    private static double[][] askZ(double[][] za,int dlugoscBitu, int[] bits, double fn, double A1, double A2, double[] t,double Tb) {
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;
        for (int i = 0; i < dlugoscBitu; i++) {
                for (int j = temp; j < temp + (Tb * t.length); j++) {
                    sygnal[i][j] = za[i][j] *  A2 *  Math.sin(2*Math.PI * fn * t[j- temp]);
                }
            temp += (Tb * t.length);
        }
        return sygnal;
    }

    private static double[][] askP(int B, double[][] askZ, double Tb , double[] t){
        int temp = 0;
        double sygnal[][] = new double[B][t.length];
        for(int i = 0; i < B; i++){
            double sum = 0;
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                sum += askZ[i][j];
                sygnal[i][j] = sum;
            }
            temp += (Tb * t.length);
    }
        return sygnal;
    }

    private static double[][] askC(double[][] askP, int B, double Tb, double[] t){
        int H = 1750;
        double[][] doubleTablica = new double[B][t.length];
        int temp = 0;
        for(int i = 0; i < B; i++){
            double sum = 0;
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                sum = askP[i][j];


                if (sum > H) {
                    doubleTablica[i][j] = 1;
                } else {
                    doubleTablica[i][j] = 0;
                }
            }

            temp += (Tb * t.length);
        }
        return doubleTablica;
    }
    private static double[][] pskZ(double[][] psk,int dlugoscBitu, int[] bits, double fn, double A1, double A2, double[] t,double Tb) {
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;
        for (int i = 0; i < dlugoscBitu; i++) {
            for (int j = temp; j < temp + (Tb * t.length); j++) {

                    sygnal[i][j] = psk[i][j] *  A2 *  Math.sin(2*Math.PI * fn * t[j- temp]);
            }
            temp += (Tb * t.length);
        }
        return sygnal;
    }

    private static double[][] pskP(int B, double[][] pskZ, double Tb , double[] t){
        int temp = 0;
        double sygnal[][] = new double[B][t.length];
        for(int i = 0; i < B; i++){
            double sum = 0;
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                sum += pskZ[i][j];
                sygnal[i][j] = sum;
            }
            temp += (Tb * t.length);
        }
        return sygnal;
    }

    private static double[][] pskC(double[][] pskP, int B, double Tb, double[] t){
        int H = 0;
        double[][] doubleTablica = new double[B][t.length];
        int temp = 0;
        for(int i = 0; i < B; i++){
            double sum = 0;
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                sum = pskP[i][j];

                if (sum < H) {
                    doubleTablica[i][j] = 1;
                } else {
                    doubleTablica[i][j] = 0;
                }
            }
            temp += (Tb * t.length);
        }
        return doubleTablica;
    }

    private static double[][] fskZ(double[][] zf,int dlugoscBitu, int[] bits, double fn1, double[] t,double Tb){
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;
        for (int i = 0; i < dlugoscBitu; i++) {
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                    sygnal[i][j] = zf[i][j] * Math.sin(2 * Math.PI * fn1 * t[j - temp]);
                }
            temp += (Tb * t.length);
            }

        return sygnal;
    }

    private static double[][] fskP(int B, double[][] fskZ, double Tb , double[] t){
        int temp = 0;
        double sygnal[][] = new double[B][t.length];
        for(int i = 0; i < B; i++){
            double sum = 0;
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                sum += fskZ[i][j];
                sygnal[i][j] = sum;
            }
            temp += (Tb * t.length);
        }
        return sygnal;
    }

    private static double[][]odejmijFsk(int B, double[][]fskZ,double[][]fskZ2,double Tb, double[] t){
        double[][] signal = new double[B][t.length];
        int temp = 0;
        for(int i = 0; i < B; i++){

            for (int j = temp; j < temp + (Tb * t.length); j++) {

            signal[i][j] = fskZ2[i][j] - fskZ[i][j];
            }


            temp += (Tb * t.length);
            }
        return signal;
    }

    private static double[][] fskC(double[][] fskP, int B, double Tb, double[] t){
        int H = 0;
        double[][] doubleTablica = new double[B][t.length];
        int temp = 0;
        for(int i = 0; i < B; i++){
            double sum = 0;

            for (int j = temp; j < temp + (Tb * t.length); j++) {
                sum =fskP[i][j];

                if (sum > H) {
                    doubleTablica[i][j] = 1;
                } else {
                    doubleTablica[i][j] = 0;
                }

            }

            temp += (Tb * t.length);
        }
        return doubleTablica;
    }

    private static int[] changeSignalToBits(int B, double Tb, double[][] signal, double[] t) {
        int[] result = new int[B];
        int temp = 0;
        for (int i = 0; i < B; i++) {
            double sum = 0;
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                sum = signal[i][j];
            }
           if(sum >= 1){
                result[i] = 1;
            } else {
                result[i] = 0;
            }
            temp += (Tb * t.length);
        }
        return result;
    }

    public static void main(String[] args) {
        double Tc = 1; // [s] czas trwania sygnału

        double fs = 1000; // [Hz] częstotliwość próbkowania

        double W = 2;

        int N = (int) (Tc * fs); // liczba próbek przypadających na cały sygnał
        double Ts = 1 / fs; // okres próbkowania

        double[] t = new double[N];
        for (int i = 0; i < N; i++) {
            t[i] = i * Ts;
        }


        StringBuilder stringBuilder = zmianaNaBity("ab");
        //1100001110   0010 = 10 bit
        String string = stringBuilder.toString();

        int[] bits = new int[stringBuilder.length()];

        char[] charArray = string.toCharArray();

        for (int i = 0; i < string.length(); i++) {
            bits[i] = Character.getNumericValue(charArray[i]);
        }
        int dlugoscBitu = bits.length;

        int B = dlugoscBitu; // dlugosc bitow
        double Tb = Tc / B;

        double fn = W * 1 / Tb;
        double A1 = 4.0;
        double A2 = 12.0;

        double fn1 = (W + 1) / Tb;
        double fn2 = (W + 2) / Tb;

        double[][] za = za(B, bits, fn, A1, A2, t,Tb);
        double[][] zf = zf(B, bits, fn1, fn2, t,Tb);
        double[][] zp = zp(B, bits, fn, t,Tb);



        double[][] askZ = askZ(za, B, bits, fn, A1, A2, t, Tb);
        double[][] pskZ = pskZ(zp, B, bits, fn, A1, A2, t, Tb);
        double[][] fskZ = fskZ(zf, B, bits,  fn2,t, Tb);
        double[][] fskZ2 = fskZ(zf, B, bits,  fn1,t, Tb);

        double[][] askP = askP(B, askZ, Tb, t);
        double[][] pskP = pskP(B, pskZ, Tb, t);

        double[][] fskP = fskP(B, zf, Tb, t);
        double[][] fskP1 = fskP(B, fskZ, Tb, t);
        double[][] fskP2 = fskP(B, fskZ2, Tb, t);

        double[][] fskPoOdjeciu = odejmijFsk(B,fskZ,fskZ2,Tb,t);



        XYSeriesCollection dataset_askZ = new XYSeriesCollection();
        XYSeriesCollection dataset_askP = new XYSeriesCollection();

        XYSeriesCollection dataset_pskZ = new XYSeriesCollection();
        XYSeriesCollection dataset_pskP = new XYSeriesCollection();

        XYSeriesCollection dataset_fskZ = new XYSeriesCollection();
        XYSeriesCollection dataset_fskZ2 = new XYSeriesCollection();

        XYSeriesCollection dataset_fskP = new XYSeriesCollection();
        XYSeriesCollection dataset_fskP2 = new XYSeriesCollection();


        XYSeriesCollection dataset_fskC = new XYSeriesCollection();
        XYSeriesCollection dataset_askC = new XYSeriesCollection();
        XYSeriesCollection dataset_pskC = new XYSeriesCollection();

        XYSeriesCollection dataset_fskP1 = new XYSeriesCollection();


        double[][] askC = askC(askP, B, Tb, t);
        double[][] pskC = pskC(pskP, B, Tb, t);
        double[][] fskC = fskC(fskPoOdjeciu, B, Tb, t);









        for (int i = 0; i < B; i++) {
            dataset_askZ.addSeries(createSeries(t, askZ[i],"bit" + i ));
            dataset_askP.addSeries(createSeries(t, askP[i],"bit" + i ));

            dataset_pskZ.addSeries(createSeries(t, pskZ[i],"bit" + i ));
            dataset_pskP.addSeries(createSeries(t, pskP[i],"bit" + i ));

            dataset_fskZ.addSeries(createSeries(t, fskZ[i],"bit" + i ));
            dataset_fskZ2.addSeries(createSeries(t, fskZ2[i],"bit" + i ));

            dataset_fskP1.addSeries(createSeries(t, fskP1[i],"bit" + i ));

            dataset_fskP.addSeries(createSeries(t, fskP[i],"bit" + i ));
            dataset_fskP2.addSeries(createSeries(t, fskP2[i],"bit" + i ));

            dataset_fskC.addSeries(createSeries(t, fskC[i],"bit" + i ));
            dataset_pskC.addSeries(createSeries(t, pskC[i],"bit" + i ));
            dataset_askC.addSeries(createSeries(t, askC[i],"bit" + i ));

        }



        JFreeChart chart_askZ = ChartFactory.createXYLineChart("Modulated Signals - askZ", "Time (s)", "Amplitude", dataset_askZ);
        JFreeChart char_askP = ChartFactory.createXYLineChart("Modulated Signals - askP", "Time (s)", "Amplitude", dataset_askP);

        JFreeChart char_pskP = ChartFactory.createXYLineChart("Modulated Signals - pskP", "Time (s)", "Amplitude", dataset_pskP);
        JFreeChart char_pskZ = ChartFactory.createXYLineChart("Modulated Signals - pskZ", "Time (s)", "Amplitude", dataset_pskZ);


        JFreeChart char_fskZ = ChartFactory.createXYLineChart("Modulated Signals - fskZ", "Time (s)", "Amplitude", dataset_fskZ);
        JFreeChart char_fskZ2 = ChartFactory.createXYLineChart("Modulated Signals - fskZ2", "Time (s)", "Amplitude", dataset_fskZ2);
        JFreeChart char_fskP = ChartFactory.createXYLineChart("Modulated Signals - fskP", "Time (s)", "Amplitude", dataset_fskP);
        JFreeChart char_fskP2 = ChartFactory.createXYLineChart("Modulated Signals - fskP2", "Time (s)", "Amplitude", dataset_fskP2);

        JFreeChart char_fskP1 = ChartFactory.createXYLineChart("Modulated Signals - fskP1", "Time (s)", "Amplitude", dataset_fskP1);


        //-------------------------------------------------------------------WYKRES askC pskC fskC-------------------------------------------------------
        JFreeChart char_fskC = ChartFactory.createXYLineChart("Modulated Signals - fskC", "Time (s)", "Amplitude", dataset_fskC);
        JFreeChart char_askC = ChartFactory.createXYLineChart("Modulated Signals - askC", "Time (s)", "Amplitude", dataset_askC);
        JFreeChart char_pskC = ChartFactory.createXYLineChart("Modulated Signals - pskC", "Time (s)", "Amplitude", dataset_pskC);


        displayChart(chart_askZ);
        displayChart(char_askP);

        displayChart(char_pskP);
        displayChart(char_pskZ);

        displayChart(char_fskZ);
        displayChart(char_fskP);

        displayChart(char_fskZ2);
        displayChart(char_fskP2);

        displayChart(char_fskC);
        displayChart(char_askC);
        displayChart(char_pskC);

        displayChart(char_fskP1);


        int[] askCToBits = changeSignalToBits(B, Tb, askC, t);
        int[] pskCToBits = changeSignalToBits(B, Tb,pskC, t);
        int[] fskCToBits = changeSignalToBits(B, Tb,fskC, t);



        for (int i = 0; i < B; i++) {
            askCToBits[i] = bits[i];
            pskCToBits[i] = bits[i];
            fskCToBits[i] = bits[i];


            if (askCToBits[i] == bits[i] && pskCToBits[i] == bits[i] && fskCToBits[i] == bits[i]) {
                System.out.println((i + 1) + " zgadzaja sie ");
            } else {
                System.out.println((i + 1) + " wielBLAD");
            }
        }


        double[][] za2 = za(B, askCToBits, fn, A1, A2, t,Tb);
        double[][] zf2 = zf(B, pskCToBits, fn1, fn2, t,Tb);
        double[][] zp2 = zp(B, fskCToBits, fn, t,Tb);


        XYSeriesCollection dataset_za = new XYSeriesCollection();
        XYSeriesCollection dataset_zp = new XYSeriesCollection();
        XYSeriesCollection dataset_zf = new XYSeriesCollection();

        XYSeriesCollection dataset_za2 = new XYSeriesCollection();
        XYSeriesCollection dataset_zf2 = new XYSeriesCollection();
        XYSeriesCollection dataset_zp2 = new XYSeriesCollection();


for(int i =0; i<B;i++){
    dataset_za.addSeries(createSeries(t, za[i],"bit" + i ));
    dataset_zp.addSeries(createSeries(t, zp[i],"bit" + i ));
    dataset_zf.addSeries(createSeries(t, zf[i],"bit" + i ));

    dataset_za2.addSeries(createSeries(t, za2[i],"bit" + i ));
    dataset_zp2.addSeries(createSeries(t, zp2[i],"bit" + i ));
    dataset_zf2.addSeries(createSeries(t, zf2[i],"bit" + i ));

}

        JFreeChart char_za = ChartFactory.createXYLineChart("Modulated Signals - ask-oryginalne", "Time (s)", "Amplitude", dataset_za);
        JFreeChart char_zp = ChartFactory.createXYLineChart("Modulated Signals - psk-oryginalne", "Time (s)", "Amplitude", dataset_zp);
        JFreeChart char_zf = ChartFactory.createXYLineChart("Modulated Signals - fsk-oryginalne", "Time (s)", "Amplitude", dataset_zf);



        JFreeChart char_za2 = ChartFactory.createXYLineChart("Modulated Signals - ask-odwzorowane", "Time (s)", "Amplitude", dataset_za2);
        JFreeChart char_zp2 = ChartFactory.createXYLineChart("Modulated Signals - psk-odwzorowane", "Time (s)", "Amplitude", dataset_zp2);
        JFreeChart char_zf2 = ChartFactory.createXYLineChart("Modulated Signals - fsk-odwzorowane", "Time (s)", "Amplitude", dataset_zf2);


        displayChart(char_za);
        displayChart(char_zf);
        displayChart(char_zp);

        displayChart(char_za2);
        displayChart(char_zf2);
        displayChart(char_zp2);

    }


}