import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;

public class Main {
    static double[] B3dB = new double[3];
    static double[] B6dB = new double[3];
    static double[] B12dB = new double[3];
    static int i =0;
    public static void main(String[] args) {
        int fm = 10;
        int fn = 100;
       double Tc = 1;
         int fs = 1000;


        double Ts = 1.0 / fs;
        //  int fm = 10;
       // int fn = 1000;

        //fs/2

        double[] kA= {0.5, 5, 25};
        double[] kP = {0.5, Math.PI / 2, 2 * Math.PI + 1};
        double[] kF = {0.5, Math.PI / 2, 2 * Math.PI + 1};


        int N = (int) (Tc * fs);
        double[] t = new double[N];

        double[] zA = new double[N];
        double[] zP = new double[N];
        double[] zF = new double[N];

        for(int i = 0; i < N;i++ ) {
            t[i] = i * Ts;
            double m = Math.sin(2 * Math.PI * fm *t[i] );

            int wariant = 0;
            zA[i] = (kA[wariant] * m + 1) * Math.cos(2 * Math.PI * fn * t[i]);
            zP[i] = Math.cos(2 * Math.PI * fn * t[i] + kP[wariant] * m);
            zF[i] = Math.cos(2 * Math.PI * fn * t[i] + ((kF[wariant]/fm) * m));
        }


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

            for (int k = 0; k < N; k++) {
                reWartosczA += zA[k] * Math.cos((2*Math.PI*n*k)/N);
                imgWartosczA += zA[k] * -Math.sin((2*Math.PI*n*k)/N);

                reWartosczP+= zP[k] * Math.cos((2*Math.PI*n*k)/N);
                imgWartosczP += zP[k] * -Math.sin((2*Math.PI*n*k)/N);

                reWartosczF+= zF[k] * Math.cos((2*Math.PI*n*k)/N);
                imgWartosczF += zF[k] * -Math.sin((2*Math.PI*n*k)/N);
            }
            rezA[n] = reWartosczA;
            imgzA[n] = imgWartosczA;

            rezP[n] = reWartosczP;
            imgzP[n] = imgWartosczP;

            rezF[n] = reWartosczF;
            imgzF[n] = imgWartosczF;

        }


        XYSeries widmozA = createXY(rezA, imgzA, fn, N, "widmozA");
        createChart(widmozA,"widmozA");

        XYSeries widmozF = createXY(rezF, imgzF, fn, N, "widmozF");
        createChart(widmozF,"widmozF");


        XYSeries widmozP = createXY(rezP, imgzP, fn, N, "widmozP");
        createChart(widmozP,"widmozP");


        createChart(zA,"zA",fn);
        createChart(zP,"zP",fn);
        createChart(zF,"zF",fn);

        String[] strinArr = {"Widmo zA","Widmo zP","Widmo zF" };

for(int i = 0; i<3; i ++) {
    double v1 = B3dB[i];
    double v2 = B6dB[i];
    double v3 = B12dB[i];

    System.out.println(strinArr[i] + ":\nDla B3dB " + v1 + '\n' + "Dla B6dB " + v2 + '\n'+ "Dla B12dB " + v3 +'\n');
}
    }

    private static XYSeries createXY(double[]re , double[]img, int fn,  int N, String tytul){
        double[] widmoAmplitudowe = new double[N/2];

        double[] czestotliwosci = new double[N/2];



        for(int k = 0; k< N/2; k++){
            double M = Math.sqrt((re[k] * re[k]) + (img[k] * img[k]));
            double MDecebel = 10 * Math.log10(M);
            double fk = Math.log(((k * fn) / N) + 1) ;
            czestotliwosci[k] = fk ;
            widmoAmplitudowe[k] = MDecebel;

        }
        XYSeries widmo = new XYSeries(tytul);

        for(int i = 0; i < N/2; i++){
            widmo.add(czestotliwosci[i],widmoAmplitudowe[i]);
        }
        B3dB[i] = obliczBdB(widmoAmplitudowe, czestotliwosci, 3);
        B6dB[i] = obliczBdB(widmoAmplitudowe, czestotliwosci, 6);
        B12dB[i++] = obliczBdB(widmoAmplitudowe, czestotliwosci, 12);
        return  widmo;
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

    private static void createChart(double[] data, String tytul, int fn) {
        XYSeries series = new XYSeries(tytul);
        for (int i = 0; i < data.length; i++) {
            series.add((double)i/fn, data[i]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createXYLineChart(
                tytul,
                "Czas",
                "Amplituda",
                dataset
        );
        ChartPanel chartPanel = new ChartPanel(chart);
        doFrame(chartPanel, tytul);
    }

    private static void createChart(XYSeries series, String tytul) {

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createXYLineChart(
                tytul,
                "f",
                "Amplituda",
                dataset
        );
        ChartPanel chartPanel = new ChartPanel(chart);
        doFrame(chartPanel, tytul);
    }



    private static void doFrame(ChartPanel chartPanel,String tytul){
        JFrame frame = new JFrame(tytul);
        frame.setLayout(new BorderLayout());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);
    }



}