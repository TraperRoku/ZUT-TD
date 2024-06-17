import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

import javax.swing.*;

public class BERPlotter {

    public static void createBarChart(double k, double berAsk, double berPsk, double berFsk) {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        dataset.addValue(berAsk, "BER", "ASK");
        dataset.addValue(berPsk, "BER", "PSK");
        dataset.addValue(berFsk, "BER", "FSK");

        JFreeChart barChart = ChartFactory.createBarChart(
                "BER for ASK, PSK, FSK (k=" + k + ")",
                "Modulation Type",
                "BER",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);

        ChartPanel chartPanel = new ChartPanel(barChart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame("BER Chart");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.getContentPane().add(chartPanel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
}
