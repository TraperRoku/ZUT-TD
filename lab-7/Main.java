import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Objects;

public class Main {



    private static StringBuilder zmianaNaBity(String string){
        StringBuilder stringBuilder = new StringBuilder();
        for(char c : string.toCharArray()){
            String bit = Integer.toBinaryString((int)c);
            stringBuilder.append(bit);
        }
        return stringBuilder;
    }


    private static int[] kodowanieHamminga74(int [] bity) {
        int x3 = bity[0];
        int x5 = bity[1];
        int x6 = bity[2];
        int x7 = bity[3];


        int x1 = ((x3 + x5) % 2 + x7) % 2;
        int x2 = ((x3 + x6) % 2 + x7) % 2;
        int x4 = ((x5 + x6) % 2 + x7) % 2;

        return new int[]{x1, x2, x3, x4, x5, x6, x7};
    }




    private static int[] dekoderHamminga74(int[] bity){
        int x1 = bity[0];
        int x2 = bity[1];
        int x3 = bity[2];
        int x4 = bity[3];
        int x5 = bity[4];
        int x6 = bity[5];
        int x7 = bity[6];

        int[] tabDoZwrotu = new int[4];

        int x1prim = ((x3+x5)%2+x7)%2;
        int x2prim = ((x3+x6)%2+x7)%2;
        int x4prim = ((x5+x6)%2+x7)%2;

        int x1zPodlogaUp = (x1 + x1prim)%2;
        int x2zPodlogaUp = (x2 + x2prim)%2;
        int x4zPodlogaUp = (x4 + x4prim)%2;

        int syndrom= x1zPodlogaUp + x2zPodlogaUp*2 + x4zPodlogaUp *4 ;

        if(syndrom != 0){
            if(bity[syndrom-1] == 1){
                bity[syndrom-1] = 0;
            }else {
                bity[syndrom - 1] = 1;
            }
            //ilosc bledow bład / iloscbitow

        }

        tabDoZwrotu[0] = bity[2];
        tabDoZwrotu[1] = bity[4];
        tabDoZwrotu[2] = bity[5];
        tabDoZwrotu[3] = bity[6];

        return tabDoZwrotu;
    }


    // Kluczowanie z przesuwem amplitudy ASK
    private static double[][] za(int dlugoscBitu, ArrayList<Integer> bits, double fn, double A1, double A2, double[] t,double Tb) {
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;

        if(A1 > A2){
            double temp2 = A1;
            A1 = A2;
            A2 = temp2;
        }
        for (int i = 0; i < dlugoscBitu; i++) {
            if (bits.get(i) == 0) {
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
    private static double[][] zp(int dlugoscBitu, ArrayList<Integer> bits, double fn, double[] t,double Tb){
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;
        for (int i = 0; i < dlugoscBitu; i++) {
            if (bits.get(i) == 0) {
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
    private static double[][] zf(int dlugoscBitu, ArrayList<Integer>  bits, double fn1, double fn2, double[] t,double Tb){
        double sygnal[][] = new double[dlugoscBitu][t.length];
        int temp = 0;
        for (int i = 0; i < dlugoscBitu; i++) {
            if (bits.get(i) == 1) {
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


    private static double[][] askZ(double[][] za,int dlugoscBitu, double fn, double A1, double A2, double[] t,double Tb) {
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
        double H = 4.5;
        double[][] doubleTablica = new double[B][t.length];
        int temp = 0;
        for(int i = 0; i < B; i++){
            double sum ;
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
    private static double[][] pskZ(double[][] psk,int dlugoscBitu,  double fn, double A1, double A2, double[] t,double Tb) {
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

    private static double[][] fskZ(double[][] zf,int dlugoscBitu, double fn1, double[] t,double Tb){
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

                signal[i][j] = fskZ2[i][j] + fskZ[i][j];
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

    private static double[][] dodawanieSzumu(double[][] signal, double alpha, double Tb , double[] t) {
        int temp = 0;
        double[][] noisySignal = new double[signal.length][signal[0].length];
        for (int i = 0; i < signal.length; i++) {
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                double randomValue = -3 + (Math.random() * 6);
                noisySignal[i][j] = signal[i][j] + (alpha * randomValue);
            }
            temp += (Tb * t.length);
        }
        return noisySignal;
    }


    private static double[][] dodawanieSzumuZadanie3(double[][] signal, double B, double[] t, double Tb) {
        double[][] noisySignal = new double[signal.length][signal[0].length];
        int temp = 0;
        for (int i = 0; i < signal.length; i++) {
            for (int j = temp; j < temp + (Tb * t.length); j++) {
                noisySignal[i][j] = signal[i][j] * Math.pow(Math.E, -B * (t[j-temp]));
            }
            temp += (Tb * t.length);
        }
        return noisySignal;
    }




    private static void sprawdzenieBERCombo(int doIlu, double ile, double[][] askKodowaneOG, double[][] fskKodowaneOG, double[][] pskKodowaneOG, double t[],
                                            double Tb , int B ,ArrayList<Integer> bityZakodowane, double fn,double A1,double A2,double fn2,double fn1, boolean flag){
        for(double k = 0;k<doIlu;k+=ile) {

            double[][] askKodowane;
            double[][]  fskKodowane;
            double[][]  pskKodowane;
            if(flag) {
                askKodowane = dodawanieSzumuZadanie3(askKodowaneOG, k*10, t, Tb);
                fskKodowane = dodawanieSzumuZadanie3(fskKodowaneOG, k*10, t, Tb);
                pskKodowane = dodawanieSzumuZadanie3(pskKodowaneOG, k*10, t, Tb);


                askKodowane = dodawanieSzumu(askKodowane, k, Tb, t);
                fskKodowane = dodawanieSzumu(fskKodowane, k, Tb, t);
                pskKodowane = dodawanieSzumu(pskKodowane, k, Tb, t);
            }else{


                askKodowane = dodawanieSzumu(askKodowaneOG, k, Tb, t);
                fskKodowane = dodawanieSzumu(fskKodowaneOG, k, Tb, t);
                pskKodowane = dodawanieSzumu(pskKodowaneOG, k, Tb, t);

                askKodowane = dodawanieSzumuZadanie3(askKodowane, k*10, t, Tb);
                fskKodowane = dodawanieSzumuZadanie3(fskKodowane, k*10, t, Tb);
                pskKodowane = dodawanieSzumuZadanie3(pskKodowane, k*10, t, Tb);

            }

            double[][] askZ = askZ(askKodowane, B, fn, A1, A2, t, Tb);
            double[][] pskZ = pskZ(pskKodowane, B,  fn, A1, A2, t, Tb);
            double[][] fskZ = fskZ(fskKodowane, B,  fn2, t, Tb);
            double[][] fskZ2 = fskZ(fskKodowane, B,  fn1, t, Tb);

            double[][] askP = askP(B, askZ, Tb, t);
            double[][] pskP = pskP(B, pskZ, Tb, t);
            double[][] fskPZ = fskP(B, fskZ, Tb, t);
            double[][] fskPZ2 = fskP(B, fskZ2, Tb, t);


            double[][] fskPoOdjeciu = odejmijFsk(B, fskPZ, fskPZ2, Tb, t);





            double[][] askC = askC(askP, B, Tb, t);
            double[][] pskC = pskC(pskP, B, Tb, t);
            double[][] fskC = fskC(fskPoOdjeciu, B, Tb, t);




            int[] askCToBits = changeSignalToBits(B, Tb, askC, t);
            int[] pskCToBits = changeSignalToBits(B, Tb,pskC, t);
            int[] fskCToBits = changeSignalToBits(B, Tb,fskC, t);



            ArrayList<Integer> bityzdekodowaneAsk = new ArrayList<>();
            ArrayList<Integer> bityzdekodowanePsk = new ArrayList<>();
            ArrayList<Integer> bityzdekodowaneFsk = new ArrayList<>();


            int[] tempArrAsk = new int[7];
            int[] tempArrPsk = new int[7];
            int[] tempArrFsk = new int[7];
            int j = 0;
            int c =0;
            for(int i = 0; i < askCToBits.length;i++){
                tempArrAsk[j] = askCToBits[c];
                tempArrFsk[j] = fskCToBits[c];
                tempArrPsk[j] = pskCToBits[c];
                c++;
                j++;
                if(j%7==0 && j!=0 ){
                    int[] codeBitsAsk = dekoderHamminga74(tempArrAsk);
                    int[] codeBitsPsk = dekoderHamminga74(tempArrPsk);
                    int[] codeBitsFsk = dekoderHamminga74(tempArrFsk);
                    for(int l = 0; l<codeBitsFsk.length;l++){
                        bityzdekodowaneAsk.add(codeBitsAsk[l]);
                        bityzdekodowaneFsk.add(codeBitsFsk[l]);
                        bityzdekodowanePsk.add(codeBitsPsk[l]);
                    }
                    tempArrAsk = new int[7];
                    tempArrPsk = new int[7];
                    tempArrFsk = new int[7];

                    j = 0;
                }
            }


            double berAsk = calculateBER(bityZakodowane, askCToBits);
            double berPsk = calculateBER(bityZakodowane, pskCToBits);
            double berFsk = calculateBER(bityZakodowane, fskCToBits);
            System.out.println();

            System.out.println("Bit Error Rate for ASK: " + berAsk + " dla "+k);
            System.out.println("Bit Error Rate for PSK: " + berPsk+ " dla "+k);
            System.out.println("Bit Error Rate for FSK: " + berFsk+ " dla "+k);


            BERPlotter.createBarChart(k, berAsk, berPsk, berFsk);
        }
    }


    private static void displayFull(int B, double[] t, double[][] signal,String string){
        XYSeriesCollection dataset_fskC = new XYSeriesCollection();
        for (int i = 0; i < B; i++) {
            dataset_fskC.addSeries(createSeries(t, signal[i], "bit" + i));
        }
        JFreeChart chart_fskC = ChartFactory.createXYLineChart(string, "Time (s)", "Amplitude", dataset_fskC);

        displayChart(chart_fskC);

    }




    private static void sprawdzenieBER(int doIlu, double ile, double[][] askKodowaneOG, double[][] fskKodowaneOG, double[][] pskKodowaneOG, double t[],
                                       double Tb , int B ,ArrayList<Integer> bityZakodowane, double fn,double A1,double A2,double fn2,double fn1, Boolean flag){
        for(double k = 0;k<=doIlu;k+=ile) {

            double[][] askKodowane ;
            double[][] fskKodowane ;
            double[][] pskKodowane ;
            if(flag){
                askKodowane = dodawanieSzumuZadanie3(askKodowaneOG, k, t, Tb);
                fskKodowane = dodawanieSzumuZadanie3(fskKodowaneOG, k, t, Tb);
                pskKodowane = dodawanieSzumuZadanie3(pskKodowaneOG, k, t, Tb);
            }else {
                askKodowane = dodawanieSzumu(askKodowaneOG, k, Tb,t);
                fskKodowane = dodawanieSzumu(fskKodowaneOG, k, Tb,t);
                pskKodowane = dodawanieSzumu(pskKodowaneOG, k, Tb,t);
            }


            double[][] askZ = askZ(askKodowane, B, fn, A1, A2, t, Tb);
            double[][] pskZ = pskZ(pskKodowane, B,  fn, A1, A2, t, Tb);
            double[][] fskZ = fskZ(fskKodowane, B,  fn2, t, Tb);
            double[][] fskZ2 = fskZ(fskKodowane, B,  fn1, t, Tb);

            double[][] askP = askP(B, askZ, Tb, t);
            double[][] pskP = pskP(B, pskZ, Tb, t);
            double[][] fskPZ = fskP(B, fskZ, Tb, t);
            double[][] fskPZ2 = fskP(B, fskZ2, Tb, t);


            double[][] fskPoOdjeciu = odejmijFsk(B, fskPZ, fskPZ2, Tb, t);





            double[][] askC = askC(askP, B, Tb, t);
            double[][] pskC = pskC(pskP, B, Tb, t);
            double[][] fskC = fskC(fskPoOdjeciu, B, Tb, t);




            int[] askCToBits = changeSignalToBits(B, Tb, askC, t);
            int[] pskCToBits = changeSignalToBits(B, Tb,pskC, t);
            int[] fskCToBits = changeSignalToBits(B, Tb,fskC, t);



            ArrayList<Integer> bityzdekodowaneAsk = new ArrayList<>();
            ArrayList<Integer> bityzdekodowanePsk = new ArrayList<>();
            ArrayList<Integer> bityzdekodowaneFsk = new ArrayList<>();


            int[] tempArrAsk = new int[7];
            int[] tempArrPsk = new int[7];
            int[] tempArrFsk = new int[7];
            int j = 0;
            int c =0;
            for(int i = 0; i < askCToBits.length;i++){
                tempArrAsk[j] = askCToBits[c];
                tempArrFsk[j] = fskCToBits[c];
                tempArrPsk[j] = pskCToBits[c];
                c++;
                j++;
                if(j%7==0 && j!=0 ){
                    int[] codeBitsAsk = dekoderHamminga74(tempArrAsk);
                    int[] codeBitsPsk = dekoderHamminga74(tempArrPsk);
                    int[] codeBitsFsk = dekoderHamminga74(tempArrFsk);
                    for(int l = 0; l<codeBitsFsk.length;l++){
                        bityzdekodowaneAsk.add(codeBitsAsk[l]);
                        bityzdekodowaneFsk.add(codeBitsFsk[l]);
                        bityzdekodowanePsk.add(codeBitsPsk[l]);
                    }
                    tempArrAsk = new int[7];
                    tempArrPsk = new int[7];
                    tempArrFsk = new int[7];

                    j = 0;
                }
            }


            double berAsk = calculateBER(bityZakodowane, askCToBits);
            double berPsk = calculateBER(bityZakodowane, pskCToBits);
            double berFsk = calculateBER(bityZakodowane, fskCToBits);
            System.out.println();

            System.out.println("Bit Error Rate for ASK: " + berAsk + " dla "+k);
            System.out.println("Bit Error Rate for PSK: " + berPsk+ " dla "+k);
            System.out.println("Bit Error Rate for FSK: " + berFsk+ " dla "+k);


            BERPlotter.createBarChart(k, berAsk, berPsk, berFsk);
        }
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


    private static double calculateBER(ArrayList<Integer> originalBits, int[] decodedBits) {
        int errorCount = 0;
        int totalBits = originalBits.size();
        for (int i = 0; i < totalBits; i++) {
            if (originalBits.get(i) != decodedBits[i]) {
                errorCount++;
            }
        }
        return (double) errorCount / totalBits;
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


        String string = stringBuilder.toString();

        ArrayList<Integer> bits = new ArrayList<>();

        for (int i = 0; i < stringBuilder.length(); i++) {
            if(string.charAt(i) == 49) {
                bits.add(1);
            }else if ( string.charAt(i) == 48){
                bits.add(0);
            }
        }
        while(bits.size() % 4 !=0){
            bits.add(0);
        }

        ArrayList<Integer> bityZakodowane = new ArrayList<>();
        int[] tempArr = new int[4];
        int j = 0;
        int k =0;
        for(int i = 0; i < bits.size();i++){
            tempArr[j] = bits.get(k);
            k++;
            j++;
            if(j%4 == 0 && j!=0 ){
                int[] codeBits = kodowanieHamminga74(tempArr);
                for(int bitsCoded : codeBits){
                    bityZakodowane.add(bitsCoded);
                }
                tempArr = new int[4];
                j = 0;
            }
        }

        int dlugosc = bityZakodowane.size();

        int B = dlugosc; // dlugosc bitow
        double Tb = Tc / B;

        double fn = W * 1 / Tb;
        double A1 = 1.0;
        double A2 = 0.5;

        double fn1 = (W + 1) / Tb;
        double fn2 = (W + 2) / Tb;

        double[][] askKodowaneOG = za(dlugosc, bityZakodowane, fn, A1, A2, t, Tb);
        double[][] fskKodowaneOG = zf(dlugosc, bityZakodowane, fn2, fn1, t, Tb);
        double[][] pskKodowaneOG = zp(dlugosc, bityZakodowane, fn, t, Tb);




   //   sprawdzenieBER(10,1,askKodowaneOG,fskKodowaneOG,pskKodowaneOG,t,Tb,B,bityZakodowane,fn,A1,A2,fn2,fn1,true);
       // sprawdzenieBER(1,0.1,askKodowaneOG,fskKodowaneOG,pskKodowaneOG,t,Tb,B,bityZakodowane,fn,A1,A2,fn2,fn1,false);

        //true == tłumienie - > szum
        //false == szum - > tłumienie

        sprawdzenieBERCombo(1,0.1,askKodowaneOG,fskKodowaneOG,pskKodowaneOG,t,Tb,B,bityZakodowane,fn,A1,A2,fn2,fn1,false);
        sprawdzenieBERCombo(1,0.1,askKodowaneOG,fskKodowaneOG,pskKodowaneOG,t,Tb,B,bityZakodowane,fn,A1,A2,fn2,fn1,true);

    }


}