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
    double fs = 1000; // Hz
    int N = (int) (tc * fs);
    double tab[] = new double[N];

    public Main() {
        generateData();
        initUI();
    }

    private void generateData() {
        for (int n = 0; n < N; n++) {
            double t = n / fs;
            tab[n] = func(t, 1000, 0);
        }
    }

    private void initUI() {
        setTitle("Sygnał sinusoidalny z modulacją amplitudy");
        setSize(800, 600);
        setDefaultCloseOperation(EXIT_ON_CLOSE);

        JPanel chartPanel = createChartPanel();
        add(chartPanel);

        setLocationRelativeTo(null);
    }

    private JPanel createChartPanel() {
        double[] t = new double[N];
        for (int i = 0; i < N; i++) {
            t[i] = i / fs;
        }

        //XYSeries seriesY = new XYSeries("FUNKCJA y");
      // XYSeries seriesZ = new XYSeries("FUNKCJA z");
        XYSeries seriesV = new XYSeries("FUNKCJA v");

        double[] y = generateFunctionY(t);
        double[] z = generateFunctionZ(t, y);
        double[] v = generateFunctionV(t, y, z);

        for (int i = 0; i < N; i++) {
            //seriesY.add(t[i], y[i]);
        //    seriesZ.add(t[i], z[i]);
            seriesV.add(t[i], v[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        //dataset.addSeries(seriesY);
        //dataset.addSeries(seriesZ);
        dataset.addSeries(seriesV);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Wykresy funkcji v(t)",
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

    private static double func(double t, double f, double fi) {
        return Math.sin(2 * Math.PI * f * t * Math.cos(3 * Math.PI * t) + t * fi);
    }

    private double[] generateFunctionY(double[] t) {
        double[] y = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            y[i] = ( func(t[i], f, fi) * t[i]) / (3 + Math.cos(20 * Math.PI * t[i]));
        }
        return y;
    }

    private double[] generateFunctionZ(double[] t, double[] y) {
        double[] z = new double[t.length];
        for (int i = 0; i < t.length; i++) {

            z[i] = Math.pow(t[i], 2) * Math.abs(func(t[i], f, fi) * y[i] - (2 / (10 + y[i])));
        }
        return z;
    }

    private double[] generateFunctionV(double[] t, double[] y, double[] z) {
        double[] v = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            v[i] = Math.pow(z[i], 3) + 3 * Math.sin(z[i] * y[i]) * Math.abs(y[i] - func(t[i],f,fi));
        }
        return v;
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            Main ex = new Main();
            ex.setVisible(true);
        });
    }
}
