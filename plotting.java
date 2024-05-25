
/**
 * @author ALPEREN ULUKAYA
 */
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.DefaultXYDataset;

import javax.swing.*;
import java.awt.*;

public class plotting {

    private static final int SIZE = 10;
    private static final double CELL_LENGTH = 0.1;
    private static final double K = 8.99e9;
    private static final double CHARGE = getLastDigitCharge();
    private static final double EMF_BATTERY = 12.0;

    private static double getLastDigitCharge() {
        String idNumber = "20220808006";
        int lastDigit = Integer.parseInt(idNumber.substring(idNumber.length() - 1));
        return lastDigit * 1e-9;
    }

    public static void main(String[] args) {
        double[][] potentials = calculatePotentials(SIZE, CHARGE, K, CELL_LENGTH);

        printPotentials(potentials);

        JFreeChart chart2D = create2DChart(potentials);
        JFreeChart chartVvsX = createVvsXChart(potentials);
        JFreeChart chartVvsDiagonal = createVvsDiagonalChart(potentials);
        JFreeChart chartEquipotential = createEquipotentialChart(potentials);

        ChartPanel chartPanel2D = new ChartPanel(chart2D);
        ChartPanel chartPanelVvsX = new ChartPanel(chartVvsX);
        ChartPanel chartPanelVvsDiagonal = new ChartPanel(chartVvsDiagonal);
        ChartPanel chartPanelEquipotential = new ChartPanel(chartEquipotential);

        JFrame frame = new JFrame("Electric Potential Distribution");
        frame.setLayout(new GridLayout(2, 2));

        frame.add(chartPanel2D);
        frame.add(chartPanelVvsX);
        frame.add(chartPanelVvsDiagonal);
        frame.add(chartPanelEquipotential);

        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(1200, 800);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);

        System.out.println();
        solveCapacitorBatteryProblem();
    }

    private static double[][] calculatePotentials(int size, double charge, double k, double cellLength) {
        double[][] potentials = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double distance = Math.sqrt(i * cellLength * i * cellLength + j * cellLength * j * cellLength);
                if (distance == 0) {
                    potentials[i][j] = Double.POSITIVE_INFINITY;
                } else {
                    potentials[i][j] = (k * charge) / distance;
                }
            }
        }
        return potentials;
    }

    private static void printPotentials(double[][] potentials) {
        System.out.println("***1***\nPotentials:\n");
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                System.out.printf("%.2f\t", potentials[i][j]);
            }
            System.out.println();
        }
    }

    private static JFreeChart create2DChart(double[][] potentials) {
        DefaultXYDataset dataset = new DefaultXYDataset();
        int size = potentials.length;
        double[][] data = new double[2][size * size];
        int index = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                data[0][index] = i * CELL_LENGTH; // x coordinate
                data[1][index] = j * CELL_LENGTH; // y coordinate
                index++;
            }
        }
        dataset.addSeries("Potential", data);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "Electric Potential Distribution",
                "x (m)",
                "y (m)",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);
        return chart;
    }

    private static JFreeChart createVvsXChart(double[][] potentials) {
        DefaultXYDataset dataset = new DefaultXYDataset();
        int size = potentials.length;
        double[][] data = new double[2][size - 1];
        for (int i = 1; i < size; i++) {
            data[0][i - 1] = i * CELL_LENGTH;
            data[1][i - 1] = potentials[i][0];
        }
        dataset.addSeries("Potential vs x", data);
        JFreeChart chart = ChartFactory.createXYLineChart(
                "Electric Potential vs x (j=0)",
                "x (m)",
                "Potential (V)",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);
        return chart;
    }

    private static JFreeChart createVvsDiagonalChart(double[][] potentials) {
        DefaultXYDataset dataset = new DefaultXYDataset();
        int size = potentials.length;
        double[][] data = new double[2][size - 1];
        for (int i = 1; i < size; i++) {
            data[0][i - 1] = i * CELL_LENGTH * Math.sqrt(2);
            data[1][i - 1] = potentials[i][i];
        }
        dataset.addSeries("Potential vs r (Diagonal)", data);
        JFreeChart chart = ChartFactory.createXYLineChart(
                "Electric Potential vs r (Diagonal)",
                "r (m)",
                "Potential (V)",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);
        return chart;
    }

    private static JFreeChart createEquipotentialChart(double[][] potentials) {
        DefaultXYDataset dataset = new DefaultXYDataset();
        int size = potentials.length;

        double maxPotential = Double.NEGATIVE_INFINITY;
        double minPotential = Double.POSITIVE_INFINITY;

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double potential = potentials[i][j];
                if (potential > maxPotential) {
                    maxPotential = potential;
                }
                if (potential < minPotential) {
                    minPotential = potential;
                }
            }
        }

        int numLevels = 5;
        double delta = (maxPotential - minPotential) / (numLevels + 1);

        for (int level = 1; level <= numLevels; level++) {
            double[][] data = new double[2][size * size];
            int index = 0;
            double targetPotential = minPotential + level * delta;

            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    double potential = potentials[i][j];

                    if (Math.abs(potential - targetPotential) < delta / 2) {
                        data[0][index] = i * CELL_LENGTH;
                        data[1][index] = j * CELL_LENGTH;
                        index++;
                    }
                }
            }
            dataset.addSeries("Equipotential " + level, data);
        }

        JFreeChart chart = ChartFactory.createScatterPlot(
                "Equipotential Lines",
                "x (m)",
                "y (m)",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);
        return chart;
    }

    private static void solveCapacitorBatteryProblem() {
        // Given data
        double dOriginal = 0.02;
        double dNew = dOriginal / 2.0;
        double EMF_BATTERY = 12.0;

        double COriginal = calculateCapacitance(dOriginal);
        double CNew = calculateCapacitance(dNew);

        double VOld = EMF_BATTERY;
        double VNew = VOld * (COriginal / CNew);

        double ratio = VNew / VOld;

        System.out.println(
                "***2***\nPhysics Problem : A parallel plate capacitor is connected with wires of negligible resistance to a battery having emf E until the\r\n"
                        +
                        "capacitor is fully charged. The battery is then disconnected from the circuit and the plates of the capacitor are\r\n"
                        +
                        "moved to half of their original separation using insulated gloves. Let Vnew be the potential difference across the\r\n"
                        +
                        "capacitor plates when the plates are moved together. Let V old be the potential difference across the capacitor\r\n"
                        +
                        "plates when connected to the battery.\r\n" +
                        "V new / V old = ?\r\n");

        System.out.println("Solution: ");

        System.out.printf("Original Separation (d0): %.2f m", dOriginal);
        System.out.println();
        System.out.printf("New Separation (d_new): %.2f m", dNew);
        System.out.println();
        System.out.printf("Original Capacitance (C0): %.2e F", COriginal);
        System.out.println();
        System.out.printf("New Capacitance (C_new): %.2e F", CNew);
        System.out.println();

        System.out.printf("Potential difference when connected to battery (V_old): %.2f V", VOld);
        System.out.println();
        System.out.printf("Potential difference when plates are moved together (V_new): %.2f V\n", VNew);
        System.out.println();

        System.out.println("(Ratio) V_new / V_old : " + ratio);
    }

    // calculate the capacitance of a capacitor
    private static double calculateCapacitance(double separation) {

        double epsilon_0 = 8.85e-12;
        double A = 1.0;
        double C = (epsilon_0 * A) / separation;
        return C;
    }
}
