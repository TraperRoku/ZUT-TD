

public class Main {
    public static void main(String[] args) {
        double[] tab = new double[4];
        tab[0] = -2;
        tab[1] = 3;
        tab[2] = 0;
        tab[3] = 1;


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


    }
}