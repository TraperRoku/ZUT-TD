import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;

public class Main {

   static double[] B3dB = new double[6];
   static double[] B6dB = new double[6];
   static double[] B12dB = new double[6];
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

    private static void displayChart(ChartPanel chartPanel,String tytul){
        JFrame frame = new JFrame(tytul);
        frame.setLayout(new BorderLayout());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);
    }


    private static XYSeries createXY(double[]re , double[]img, double fn,  int N, String tytul,int licznikDoBXdB) {
        double[] widmoAmplitudowe = new double[N / 2];

        double[] czestotliwosci = new double[N / 2];


        for (int k = 0; k < N / 2; k++) {
            double M = Math.sqrt((re[k] * re[k]) + (img[k] * img[k]));
            double MDecebel = 10 * Math.log10(M);
            double fk = Math.log(((k * fn) / N) + 1);
            czestotliwosci[k] = fk;
            widmoAmplitudowe[k] = MDecebel;

        }
        XYSeries widmo = new XYSeries(tytul);

        for (int i = 0; i < N / 2; i++) {
            widmo.add(czestotliwosci[i], widmoAmplitudowe[i]);
        }
        B3dB[licznikDoBXdB] = obliczBdB(widmoAmplitudowe,czestotliwosci,3);
        B6dB[licznikDoBXdB] = obliczBdB(widmoAmplitudowe,czestotliwosci,6);
        B12dB[licznikDoBXdB] = obliczBdB(widmoAmplitudowe,czestotliwosci,12);

        return widmo;
    }

    private static void createChart(XYSeries series, String tytul) {

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createXYLineChart(
                tytul,
                "f[Hz]",
                "Amplituda",
                dataset
        );
        ChartPanel chartPanel = new ChartPanel(chart);
        displayChart(chartPanel, tytul);
    }



    private static double obliczBdB(double[] widmoAmplitudowe, double[] freq, double dB) {
        double maxAmp = Double.MIN_VALUE;
        double fMin = freq[0];
        double fMax = freq[freq.length - 1];


        for (int i = 0; i < widmoAmplitudowe.length; i++) {
            if (widmoAmplitudowe[i] > maxAmp) {
                maxAmp = widmoAmplitudowe[i] - 10 * Math.log10(dB);
            }
        }

        for (int i = 0; i < widmoAmplitudowe.length; i++) {
            if (widmoAmplitudowe[i] > maxAmp) {
                fMin = freq[i];
                break;
            }
        }
        for (int i = widmoAmplitudowe.length - 1; i >= 0; i--) {
            if (widmoAmplitudowe[i] > maxAmp) {
                fMax = freq[i];
                break;
            }
        }
        return fMax - fMin;
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


         XYSeriesCollection dataset_za = new XYSeriesCollection();
        XYSeriesCollection dataset_zf = new XYSeriesCollection();
 XYSeriesCollection dataset_zp = new XYSeriesCollection();


 for (int i = 0; i < B; i++) {
     dataset_za.addSeries(createSeries(t, za[i],"bit" + i ));
     dataset_zp.addSeries(createSeries(t, zp[i],"bit" + i ));
     dataset_zf.addSeries(createSeries(t, zf[i],"bit" + i ));
 }



 JFreeChart chart_za = ChartFactory.createXYLineChart("Modulated Signals - za", "Time (s)", "Amplitude", dataset_za);
 JFreeChart chart_zf = ChartFactory.createXYLineChart("Modulated Signals - zf", "Time (s)", "Amplitude", dataset_zf);
 JFreeChart chart_zp = ChartFactory.createXYLineChart("Modulated Signals - zp", "Time (s)", "Amplitude", dataset_zp);


 displayChart(chart_za);
 displayChart(chart_zf);
 displayChart(chart_zp);

 //-------------------------------------------
 //-------------------3-----------------------
        double[] rezA = new double[N];
        double[] imgzA = new double[N];

        double[] rezP= new double[N];
        double[] imgzP = new double[N];

        double[] rezF = new double[N];
        double[] imgzF = new double[N];

        for(int n = 0; n<N;n++) {
            double reWartosczA = 0, imgWartosczA = 0;
            double reWartosczP = 0, imgWartosczP = 0;
            double reWartosczF = 0, imgWartosczF = 0;

        for(int i = 0; i < dlugoscBitu; i++){
            for (int k = 0; k < N; k++) {
                reWartosczA += za[i][k] * Math.cos((2*Math.PI*n*k)/N);
                imgWartosczA += za[i][k] * -Math.sin((2*Math.PI*n*k)/N);

                reWartosczP+= zp[i][k] * Math.cos((2*Math.PI*n*k)/N);
                imgWartosczP += zp[i][k] * -Math.sin((2*Math.PI*n*k)/N);

                reWartosczF+= zf[i][k] * Math.cos((2*Math.PI*n*k)/N);
                imgWartosczF += zf[i][k] * -Math.sin((2*Math.PI*n*k)/N);
            }
            rezA[n] = reWartosczA;
            imgzA[n] = imgWartosczA;

            rezP[n] = reWartosczP;
            imgzP[n] = imgWartosczP;

            rezF[n] = reWartosczF;
            imgzF[n] = imgWartosczF;
        }
        }
        XYSeries widmozA = createXY(rezA, imgzA, fn, N, "widmozA",0);
        createChart(widmozA,"widmozA");

        XYSeries widmozF = createXY(rezF, imgzF, fn, N, "widmozF",1);
        createChart(widmozF,"widmozF");


        XYSeries widmozP = createXY(rezP, imgzP, fn, N, "widmozP",2);
        createChart(widmozP,"widmozP");

//---------------------------------------------------------------------------
// ---------------------------------4----------------------------------------
        String stringtab[] = {"widmo zA","widmo zP", "widmo zF"};

for(int i = 0; i < 3; i++) {
   System.out.println(B3dB[i] + " B3dB szerokosc pasma dla "+stringtab[i]);
   System.out.println(B6dB[i] + " B6dB szerokosc pasma dla "+stringtab[i]);
   System.out.println(B12dB[i] + " B12dB szerokosc pasma dla "+stringtab[i]);
}
    }


    }