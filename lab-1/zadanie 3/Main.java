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



        XYSeries series = new XYSeries("FUNKCJA u");

        for (int i = 0; i < N; i++) {
            series.add(t[i], tab[i]);

        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Wykresy funkcji u(t)",
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



    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            Main ex = new Main();
            ex.setVisible(true);
        });
    }
}
