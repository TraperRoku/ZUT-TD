

public class Main {


    //4 liczba bitow niosa informacje
    private static int[] kodowanieHamminga74(int [] bity) {
        int x3 = bity[0];
        int x5 = bity[1];
        int x6 = bity[2];
        int x7 = bity[3];


        int x1 = ((x3 + x5) % 2 + x7) % 2;
        int x2 = ((x3 + x6) % 2 + x7) % 2;
        int x4 = ((x5 + x6) % 2 + x7) % 2;

        return new int[]{x1, x2, x3, x4, x5, x6, x7};
    }


    private static int[] dekoderHamminga74(int[] bity){
        int x1 = bity[0];
        int x2 = bity[1];
        int x3 = bity[2];
        int x4 = bity[3];
        int x5 = bity[4];
        int x6 = bity[5];
        int x7 = bity[6];

        int[] tabDoZwrotu = new int[4];


        int x1prim = ((x3+x5)%2+x7)%2;
        int x2prim = ((x3+x6)%2+x7)%2;
        int x4prim = ((x5+x6)%2+x7)%2;

        int x1zPodlogaUp = (x1 + x1prim)%2;
        int x2zPodlogaUp = (x2 + x2prim)%2;
        int x4zPodlogaUp = (x4 + x4prim)%2;

        int syndrom= x1zPodlogaUp + x2zPodlogaUp*2 + x4zPodlogaUp *4 ;

        if(syndrom != 0){
            if(bity[syndrom-1] == 1){
                bity[syndrom-1] = 0;
            }else {
                bity[syndrom - 1] = 1;
            }
            }

        tabDoZwrotu[0] = bity[2];
        tabDoZwrotu[1] = bity[4];
        tabDoZwrotu[2] = bity[5];
        tabDoZwrotu[3] = bity[6];

        return tabDoZwrotu;
    }


    private static int[][] generateP() {
        int k = 11;
        int r = 4; // m = 4 dla kodu Hamminga (15, 11)

        int[][] P = new int[k][r];

        for (int i = 0; i < k; i++) {
            int bitValue = i + 1;
            for (int j = 0; j < r; j++) {
                P[i][j] = (bitValue >> j) & 1;
            }
        }
        return P;
    }

    private static int[][] generateG() {
        int k = 11;
        int n = 15;
        int[][] I = new int[k][k];
        int[][] P = generateP();

        for (int i = 0; i < k; i++) {
            I[i][i] = 1;
        }

        int[][] G = new int[k][n];
        for (int i = 0; i < k; i++) {
            System.arraycopy(I[i], 0, G[i], 0, k);
            System.arraycopy(P[i], 0, G[i], k, n - k);
        }
        return G;
    }

    private static int[][] generateH(int[][] G) {
        int k = G.length;
        int n = G[0].length;
        int r = n - k;

        int[][] P = new int[k][r];
        for (int i = 0; i < k; i++) {
            System.arraycopy(G[i], k, P[i], 0, r);
        }

        int[][] PTranspose = new int[r][k];
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < r; j++) {
                PTranspose[j][i] = P[i][j];
            }
        }

        int[][] I = new int[r][r];
        for (int i = 0; i < r; i++) {
            I[i][i] = 1;
        }

        int[][] H = new int[r][n];
        for (int i = 0; i < r; i++) {
            System.arraycopy(I[i], 0, H[i], 0, r);
            System.arraycopy(PTranspose[i], 0, H[i], r, k);
        }

        return H;
    }

    private static int[] kodowanieHamminga1511(int[] bity, int[][] G) {
        int[] kod = new int[15];
        for (int i = 0; i < 15; i++) {
            kod[i] = 0;
            for (int j = 0; j < 11; j++) {
                kod[i] += bity[j] * G[j][i];
            }
            kod[i] %= 2;
        }
        return kod;
    }

    private static int[] dekoderHamminga1511(int[] bity, int[][] H) {
        int[] syndrom = new int[4];
        for (int i = 0; i < 4; i++) {
            syndrom[i] = 0;
            for (int j = 0; j < 15; j++) {
                syndrom[i] += bity[j] * H[i][j];
            }
            syndrom[i] %= 2;
        }

        int syndromWartosc = syndrom[0] + syndrom[1] * 2 + syndrom[2] * 4 + syndrom[3] * 8;

        if (syndromWartosc != 0) {
            if(bity[syndromWartosc-1] == 1) {
                bity[syndromWartosc - 1] = (bity[syndromWartosc - 1] == 1) ? 0 : 1;
            }
        }

        int[] wynik = new int[11];
        for (int i = 0; i < 11; i++) {
            wynik[i] = bity[i];
        }

        return wynik;
    }


    public static void main(String[] args) {
        int[] podstawowaTablica = {1, 1, 0, 1};
        int[] tab = kodowanieHamminga74(podstawowaTablica);

        System.out.println("Kodowanie wartosci bitowej 1, 1 ,0 ,1");
        for(int i : tab) {
            System.out.print(i+" ");
        }
        System.out.println();
        int[] tab2 = dekoderHamminga74(tab);
        System.out.println("dekodowanie ");
        for(int i : tab2){
            System.out.print(i +" ");
        }

        int[] tablicka = {1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0};

        int[][] G = generateG();
        int[][] H = generateH(G);



        int[] zakodowane = kodowanieHamminga1511(tablicka, G);

        System.out.println("Zakodowana wiadomość:");
        for (int i : zakodowane) {
            System.out.print(i + " ");
        }
        System.out.println();


        int[] zdekodowane = dekoderHamminga1511(zakodowane, H);

        System.out.println("Zdekodowana wiadomość:");
        for (int i : zdekodowane) {
            System.out.print(i + " ");
        }
    }
}