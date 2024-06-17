import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;



public  class Main {

    static double f = 1000;
    static double fi = 0;

    static double tc = 1; // sekundy
    static double fs = 2205; // Hz
    static int N = (int) (tc * fs);

    public static void main(String[] args) {

        double[] tab = new double[N];
        double[] tab2 = new double[N];
        double[] tk = new double[N];

        double[][] bk = new double[3][N];
        int H1 = 5;
        int H2 = 20;
        int H3 = 50;




        double[] t = new double[N];

        for (int n = 0; n < N; n++) {
            t[n] = n / fs;
            double tx = n / fs;
            bk[0][n] = CalcEpsilonThing(t[n],H1);
            bk[1][n] = CalcEpsilonThing(t[n],H2);
            bk[2][n] = CalcEpsilonThing(t[n], H3);

            tab2[n] = func(tx,f,fi); //x
        }


        double[] y = generateFunctionY(t);
        double[] z = generateFunctionZ(t, y);
        double[] v = generateFunctionV(t, y, z);




        double[] re = new double[N];
        double[] img = new double[N];



        //Re and Imagin
        for(int n = 0; n<N;n++) {

          //  double[] syg = bk[2];

           // double t1 = n / fs;
            //  tab[n] = funcU(t1);

            double reWartosc = 0;

             double imgWartosc = 0;
            for (int k = 0; k < N; k++) {

             //  reWartosc += syg[k] * Math.cos((2*Math.PI*n*k)/N); // funckja dla  H
               // imgWartosc += syg[k] * -Math.sin((2*Math.PI*n*k)/N); // funckja dla  H

             //   reWartosc += tab[k] * Math.cos((2*Math.PI*n*k)/N); // funckja dla U
             //   imgWartosc += tab[k] * -Math.sin((2*Math.PI*n*k)/N); // funkcja dla U


                // reWartosc += y[k] * Math.cos((2*Math.PI*n*k)/N); // funkcja dla y
               // imgWartosc += y[k] * -Math.sin((2*Math.PI*n*k)/N);

            //    reWartosc += z[k] * Math.cos((2*Math.PI*n*k)/N); // funkcja dla z
              //   imgWartosc += z[k] * -Math.sin((2*Math.PI*n*k)/N);

                reWartosc += tab2[k] * Math.cos((2*Math.PI*n*k)/N); // funkcja dla x
                imgWartosc += tab2[k] * -Math.sin((2*Math.PI*n*k)/N);

            }
            re[n] = reWartosc;
            img[n] = imgWartosc;
        }


        double[] widmoAmplitudowe = new double[N/2];

        double[] czestotliwosci = new double[N/2];

        for(int k = 0; k< N/2; k++){
            double M = Math.sqrt((re[k] * re[k]) + (img[k] * img[k]));
            double MDecebel = 10 * Math.log10(M);
            double fk = k * fs/N;

            czestotliwosci[k] = fk;

            widmoAmplitudowe[k] = MDecebel;
        }

        XYSeries series = new XYSeries("Funkcja V");

        for(int i = 0; i < N/2; i++){
            series.add(czestotliwosci[i],widmoAmplitudowe[i]);
        }

        XYSeriesCollection dane = new XYSeriesCollection(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Funkcja V",
                "Czestotliowsc",
                "Amplituda w db",
                dane
        );
        ChartPanel chartPanel = new ChartPanel(chart);

        JFrame frame  = new JFrame("Funckja V");
        frame.setLayout(new BorderLayout());
       frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);

    }

    private static double CalcEpsilonThing(double t, int H) {
        double wynikZEpsilona = 0;

        for( int h = 1; h <= H; h++){
            wynikZEpsilona+= ((Math.pow(-1,h)/ h) * Math.sin(h * Math.PI * 2 * t));
        }
        return wynikZEpsilona * (2/Math.PI);

    }



    private static double func(double t, double f, double fi) {
        return Math.sin(2 * Math.PI * f * t * Math.cos(3 * Math.PI * t) + t * fi);
    }
    private static double funcU(double t){
        if(t < 1.8 && t>=0){
            return Math.sin(12*Math.cos(Math.PI*t)+Math.pow(t,2));
        }
        if(t < 2.3 && t>=1.8){
            return 3 * ( t-1.7) * Math.sin(3*Math.PI * t) * Math.cos(20*Math.pow(t,2));
        }
        if(t < 3 && t>=2.3){
            return (Math.pow(t,3)/16) * Math.sin(8* Math.PI * t);
        }
        if(t < 3.5 && t >= 3){
            return Math.log(t) /(2 + Math.sin(4 *Math.PI * t));
        }
        else return 0;
    }


    private static double[] generateFunctionY(double[] t) {
        double[] y = new double[ N];
        for (int i = 0; i <  N; i++) {
            y[i] = ( func(t[i], f, fi) * t[i]) / (3 + Math.cos(20 * Math.PI * t[i]));
        }
        return y;
    }

    private static double[] generateFunctionZ(double[] t, double[] y) {
        double[] z = new double[ N];
        for (int i = 0; i <  N; i++) {

            z[i] = Math.pow(t[i], 2) * Math.abs(func(t[i], f, fi) * y[i] - (2 / (10 + y[i])));
        }
        return z;
    }

    private static double[] generateFunctionV(double[] t, double[] y, double[] z) {
        double[] v = new double[N];
        for (int i = 0; i < N; i++) {
            v[i] = Math.pow(z[i], 3) + 3 * Math.sin(z[i] * y[i]) * Math.abs(y[i] - func(t[i],f,fi));
        }
        return v;
    }

}