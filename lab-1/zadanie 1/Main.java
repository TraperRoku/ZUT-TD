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


public class Main extends JFrame {
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
            tab[n] = func(t, 500, 0);
        }
    }

    private void initUI() {
        setTitle("Zadanie 1");
        setSize(800, 600);
        setDefaultCloseOperation(EXIT_ON_CLOSE);

        JPanel chartPanel = createChartPanel();
        add(chartPanel);

        setLocationRelativeTo(null);
    }

    private JPanel createChartPanel() {
        XYSeries series = new XYSeries("FUNKCJA 5");
        for (int i = 0; i < N; i++) {
            series.add(i / fs, tab[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Wykres",
                "Czas  [s]",
                "Wartosc ",
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

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            Main ex = new Main();
            ex.setVisible(true);
        });
    }
}
