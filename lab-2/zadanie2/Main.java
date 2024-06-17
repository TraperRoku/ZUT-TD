import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;

public class Main {
    public static void main(String[] args) {
        double[] tab = new double[4];
        tab[0] = -2;
        tab[1] = 3;
        tab[2] = 0;
        tab[3] = 1;

        double fs = 1000;

        double[] re = new double[4];
        double[] img = new double[4];

        //Re and Imagin
        for(int n = 0; n<tab.length;n++) {
             double reWartosc = 0;
             double imgWartosc = 0;
            for (int k = 0; k < tab.length; k++) {

               reWartosc += tab[k] * Math.cos((2*Math.PI*n*k)/tab.length);
                imgWartosc += tab[k] * -Math.sin((2*Math.PI*n*k)/tab.length);
            }
            re[n] = reWartosc;
            img[n] = imgWartosc;
        }

        int N = tab.length;
        double[] widmoAmplitudowe = new double[N/2];

        double[] czestotliwosci = new double[N/2];

        for(int k = 0; k< N/2; k++){
            double M = Math.sqrt((re[k] * re[k]) + (img[k] * img[k]));
            double MDecebel = 10 * Math.log10(M);
            double fk = k * fs/N;

            czestotliwosci[k] = fk;

            widmoAmplitudowe[k] = MDecebel;
        }

        XYSeries series = new XYSeries("widmo amplitudowe");

        for(int i = 0; i < N/2; i++){
            series.add(czestotliwosci[i],widmoAmplitudowe[i]);
        }

        XYSeriesCollection dane = new XYSeriesCollection(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Widmo amplitudowe",
                "Czestotliowsc",
                "Amplituda w db",
                dane
        );

        ChartPanel chartPanel = new ChartPanel(chart);

        JFrame frame  = new JFrame("Widmo amplitudowe");
        frame.setLayout(new BorderLayout());
       frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);
    }
}