package org.example;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;

public class Main extends JFrame {
    double f = 1000;
    double fi = 0;
    double tc = 1; // sekundy
    double fs = 2205; // Hz
    int N = (int) (tc * fs);
    double tab[] = new double[N];

    double[][] bk = new double[3][N];
    int H1 = 5;
    int H2 = 20;
    int H3 = 50;



    public Main() {
        generateData();
        initUI();
    }

    private void generateData() {
        double [] t = new double[N];

        for (int n = 0; n < N; n++) {
            t[n] = n / fs;
            //bk[0][n] = CalcEpsilonThing(t[n],H1);
           //bk[1][n] = CalcEpsilonThing(t[n],H2);
             bk[2][n] = CalcEpsilonThing(t[n],H3);

        }
    }

    private void initUI() {
        setSize(800, 600);
        setDefaultCloseOperation(EXIT_ON_CLOSE);

        JPanel chartPanel = createChartPanel();
        add(chartPanel);

        setLocationRelativeTo(null);
    }

    private JPanel createChartPanel() {
        double[] t = new double[N];




            XYSeries series = new XYSeries("Funkcja h3(t)");
            for (int i = 0; i < N; i++) {
                t[i] = i / fs;
                series.add(t[i], bk[2][i]);

        }



        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Wykres funkcji h3(t)",
                "Czas [s]",
                "Wartość",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYPlot plot = chart.getXYPlot();
        NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setAutoRangeIncludesZero(false);

        return new ChartPanel(chart);
    }

    private static double CalcEpsilonThing(double t, int H) {
        double wynikZEpsilona = 0;

      for( int h = 1; h <= H; h++){
          wynikZEpsilona+= ((Math.pow(-1,h)/ h) * Math.sin(h * Math.PI * 2 * t));
      }
      return wynikZEpsilona * (2/Math.PI);

    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            Main ex = new Main();
            ex.setVisible(true);
        });
    }
}
