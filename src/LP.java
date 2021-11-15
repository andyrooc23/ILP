import java.util.Arrays;

public class LP {

    public static void main(String[] args) {

        double[][] a = new double[3][3];
        a[0][0] = 32.0;
        a[0][1] = 87.0;
        a[0][2] = 59.0;
        a[1][0] = 24.0;
        a[1][1] = 31.0;
        a[1][2] = 11.0;
        a[2][0] = 97.0;
        a[2][1] = 61.0;
        a[2][2] = 30.0;
        double[] b = new double[3];
        b[0] = 50.0;
        b[1] = 129.0;
        b[2] = 20.0;
        double[] s = new double[3];
        s[0] = 37.0;
        s[1] = 38.0;
        s[2] = 10.0;
        double[] ct = new double[3];
        ct[0] = 82.0;
        ct[1] = 97.0;
        ct[2] = 11.0;

        System.out.println(Arrays.toString(solve(a, b, s, ct)));
    }

    public static double[] solve(double[][] a, double[] b, double[] s, double[] ct) {

        int n = b.length;
        double[] ls = new double[n];

        for (int i = 0; i < n; i++) {
            ls[i] = b[i] - s[i];
        }

        //invert a
        double[][] invA = invert(a);

        // xvals = dot product of (inverse of a and ls)
        double[] xvals = dot(invA, ls);

        double biggest = xvals[0];
        for (int i = 0; i < xvals.length; i++) {
            if (ldot(ct, xvals[i]) > ldot(ct, biggest))
                biggest = xvals[i];
        }

        return xvals;
    }

    public static double[] dot(double a[][], double ls[]) {
        double[] prodA = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            double product = 0;
            for (int j = 0; j < a[i].length; j++) {
                product = product + a[i][j] * ls[i];
            }
            prodA[i] = product;
        }
        return prodA;
    }

    public static double ldot(double a[], double ls) {
        double product = 0;
        for (int i = 0; i < a.length; i++) {
            product = product + a[i] * ls;
        }
        return product;
    }

    // excerpt from java library for inverting a matrix
    public static double[][] invert(double[][] a) {
        int n = a.length;
        double[][] x = new double[n][n];
        double[][] b = new double[n][n];
        int[] index = new int[n];

        for (int i = 0; i < n; i++) {
            b[i][i] = 1;
        }

        gaussian(a, index);

        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                            -= a[index[j]][i]*b[index[i]][k];

        // Perform backward substitutions
        for (int i=0; i<n; ++i)
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j)
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k)
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }

    // excerpt from java library for gaussian transformation
    public static void gaussian(double a[][], int index[]) {
        int n = index.length;
        double c[] = new double[n];

        // Initialize the index
        for (int i=0; i<n; ++i)
            index[i] = i;

        // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i)
        {
            double c1 = 0;
            for (int j=0; j<n; ++j)
            {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1)
                    c1 = c0;
            }
            c[i] = c1;
        }

        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j)
        {
            double pi1 = 0;
            for (int i=j; i<n; ++i) {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1) {
                    pi1 = pi0;
                    k = i;
                }
            }

            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i)
            {
                double pj = a[index[i]][j]/a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;

                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }
}
