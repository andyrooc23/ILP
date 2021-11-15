import java.util.Arrays;

public class ILP {

    public static void main(String[] args) {
        int[][] a = new int[3][3];
        a[0][0] = 32;
        a[0][1] = 87;
        a[0][2] = 59;
        a[1][0] = 24;
        a[1][1] = 31;
        a[1][2] = 11;
        a[2][0] = 97;
        a[2][1] = 61;
        a[2][2] = 30;
        int[] b = new int[3];
        b[0] = 50;
        b[1] = 129;
        b[2] = 20;
        int[] s = new int[3];
        s[0] = 37;
        s[1] = 38;
        s[2] = 10;
        int[] ct = new int[3];
        ct[0] = 82;
        ct[1] = 97;
        ct[2] = 11;

        System.out.println(Arrays.toString(solve(a, b, s, ct)));
    }

    public static int[] solve(int[][] a, int[] b, int[] s, int[] ct) {

        int n = b.length;
        int[] ls = new int[n];

        for (int i = 0; i < n; i++) {
            ls[i] = b[i] - s[i];
        }

        //invert a
        int[][] invA = invert(a);

        // xvals = dot product of (inverse of a and ls)
        int[] xvals = dot(invA, ls);

        int biggest = xvals[0];
        for (int i = 0; i < xvals.length; i++) {
            if (ldot(ct, xvals[i]) > ldot(ct, biggest))
                biggest = xvals[i];
        }

        return xvals;
    }

    public static int[] dot(int a[][], int ls[]) {
        int[] prodA = new int[a.length];
        for (int i = 0; i < a.length; i++) {
            int product = 0;
            for (int j = 0; j < a[i].length; j++) {
                product = product + a[i][j] * ls[i];
            }
            prodA[i] = product;
        }
        return prodA;
    }

    public static int ldot(int a[], int ls) {
        int product = 0;
        for (int i = 0; i < a.length; i++) {
            product = product + a[i] * ls;
        }
        return product;
    }

    // excerpt from java library for inverting a matrix
    public static int[][] invert(int[][] a) {
        int n = a.length;
        int[][] x = new int[n][n];
        int[][] b = new int[n][n];
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
    public static void gaussian(int a[][], int index[]) {
        int n = index.length;
        int c[] = new int[n];

        // Initialize the index
        for (int i=0; i<n; ++i)
            index[i] = i;

        // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i)
        {
            int c1 = 0;
            for (int j=0; j<n; ++j)
            {
                int c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }

        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j)
        {
            int pi1 = 0;
            for (int i=j; i<n; ++i)
            {
                int pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1)
                {
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
                int pj = a[index[i]][j]/a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;

                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }

}
