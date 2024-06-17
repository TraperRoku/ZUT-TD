package org.example;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import org.jfree.chart.ChartUtils;

public class Main {
    public static void main(String[] args) {
        XYSeries s = new XYSeries("Funkcja");

        for (double x = -20; x <= 20; x += 0.1) {
            s.add(x, funkcja(x));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(s);
        JFreeChart chart = ChartFactory.createXYLineChart(
                "FUNKCJA", "x", "y", dataset
        );

        chart.setBackgroundPaint(Color.white);
        chart.getXYPlot().setBackgroundPaint(Color.lightGray);
        chart.getXYPlot().setDomainGridlinePaint(Color.white);
        chart.getXYPlot().setRangeGridlinePaint(Color.white);


        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));


        JFrame frame = new JFrame("FUNKCJA");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);

        try {
            ChartUtils.saveChartAsPNG(new File("f.png"), chart, 800, 600);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private static double funkcja(double x) {
        return x * x + 4 * x - 10;
    }
}
