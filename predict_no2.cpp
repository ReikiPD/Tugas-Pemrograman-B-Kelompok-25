#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>    // for system()

using namespace std;

// ----------------------
// 1. Lagrange Interpolation
// ----------------------
double lagrange(const vector<double>& x, const vector<double>& y, double xi) {
    int n = x.size();
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        double term = y[i];
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                term *= (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }
    return result;
}

// ----------------------
// 2. Natural Cubic Spline
// ----------------------
struct Spline {
    vector<double> a, b, c, d, xData;
};

Spline buildSpline(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<double> a = y, b(n), c(n), d(n), h(n-1), alpha(n), l(n), mu(n), z(n);

    // compute h and alpha
    for (int i = 0; i < n-1; ++i) h[i] = x[i+1] - x[i];
    alpha[0] = alpha[n-1] = 0;
    for (int i = 1; i < n-1; ++i)
        alpha[i] = (3.0/h[i])*(a[i+1]-a[i]) - (3.0/h[i-1])*(a[i]-a[i-1]);

    // solve tridiagonal
    l[0]=1; mu[0]=z[0]=0;
    for (int i = 1; i < n-1; ++i) {
        l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i]  = (alpha[i] - h[i-1]*z[i-1]) / l[i];
    }
    l[n-1]=1; z[n-1]=c[n-1]=0;
    for (int j = n-2; j >= 0; --j) {
        c[j] = z[j] - mu[j]*c[j+1];
        b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3.0;
        d[j] = (c[j+1]-c[j]) / (3*h[j]);
    }

    return Spline{a,b,c,d,x};
}

double evalSpline(const Spline& s, double xi) {
    int n = s.xData.size(), i = n-2;
    if      (xi <= s.xData.front())    i = 0;
    else if (xi >= s.xData.back())     i = n-2;
    else {
        for (int j=0; j<n-1; ++j)
            if (xi>=s.xData[j] && xi<=s.xData[j+1]) { i=j; break; }
    }
    double dx = xi - s.xData[i];
    return s.a[i] + s.b[i]*dx + s.c[i]*dx*dx + s.d[i]*dx*dx*dx;
}

// ----------------------
// 3. Least-Squares Quadratic Regression
// ----------------------
struct QuadCoeff { double A, B, C; };

QuadCoeff quadraticRegression(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    double Sx=0, Sx2=0, Sx3=0, Sx4=0, Sy=0, Sxy=0, Sx2y=0;
    for (int i=0; i<n; ++i) {
        double xi=x[i], yi=y[i];
        Sx+=xi; Sx2+=xi*xi; Sx3+=xi*xi*xi; Sx4+=xi*xi*xi*xi;
        Sy+=yi; Sxy+=xi*yi; Sx2y+=xi*xi*yi;
    }
    double M[3][4] = {
        { (double)n, Sx,  Sx2,  Sy   },
        { Sx,        Sx2, Sx3,  Sxy  },
        { Sx2,       Sx3, Sx4,  Sx2y }
    };
    // Gaussian elimination
    for (int k=0; k<3; ++k) {
        double piv = M[k][k];
        for (int j=k; j<4; ++j) M[k][j]/=piv;
        for (int i=k+1; i<3; ++i) {
            double f = M[i][k];
            for (int j=k; j<4; ++j)
                M[i][j] -= f * M[k][j];
        }
    }
    double sol[3];
    for (int i=2; i>=0; --i) {
        sol[i] = M[i][3];
        for (int j=i+1; j<3; ++j)
            sol[i] -= M[i][j] * sol[j];
    }
    return {sol[0], sol[1], sol[2]};
}

double evalQuad(const QuadCoeff& qc, double xi) {
    return qc.A + qc.B*xi + qc.C*xi*xi;
}

// ----------------------
// Main
// ----------------------
int main() {
    vector<double> waktu = {8,10,12,14,16}, no2 = {40,55,65,70,60};
    Spline spline = buildSpline(waktu,no2);
    QuadCoeff qc  = quadraticRegression(waktu,no2);

    cout<<fixed<<setprecision(4);
    cout<<"=== Prediksi Konsentrasi NO2 ===\n"
        <<"Masukkan waktu (jam, misal 9.5), atau non-numeric untuk keluar.\n";
    while (true) {
        cout<<"\n>> Waktu = ";
        double t;
        if (!(cin>>t)) break;
        cout<<"Lagrange interp.    : "<<lagrange(waktu,no2,t)<<" ppb\n"
            <<"Natural cubic spline: "<<evalSpline(spline,t)  <<" ppb\n"
            <<"Quadratic regresi   : "<<evalQuad(qc,t)      <<" ppb\n";
    }

    // tulis CSV
    ofstream fout("predict_no2_data.csv");
    fout<<"time,lagrange,spline,quadratic\n";
    for (double t=8.0; t<=16.0; t+=0.25) {
        fout<<t<<","<<lagrange(waktu,no2,t)<<","
                  <<evalSpline(spline,t)<<","
                  <<evalQuad(qc,t)<<"\n";
    }
    fout.close();
    cout<<"\nData CSV tersimpan di predict_no2_data.csv\n";

    // tulis skrip Gnuplot
    ofstream gp("plot_predict_no2.gp");
    gp<<"set datafile separator \",\"\n"
      <<"set title 'Prediksi Konsentrasi NO_{2} terhadap Waktu'\n"
      <<"set xlabel 'Waktu (jam)'\n"
      <<"set ylabel 'NO2 (ppb)'\n"
      <<"set grid\n"
      <<"set xrange [8:16]\n"
      <<"plot \\\n"
      <<"  'predict_no2_data.csv' using 1:2 with lines lw 2 title 'Lagrange', \\\n"
      <<"  ''                        using 1:3 with lines lw 2 title 'Spline Kubik', \\\n"
      <<"  ''                        using 1:4 with lines lw 2 title 'Regresi Kuadrat'\n"
      <<"pause -1\n";
    gp.close();
    cout<<"Skrip gnuplot tersimpan di plot_predict_no2.gp\n";

    // panggil Gnuplot dengan path penuh (Windows)
    const char* gnuplotPath = "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"";
    string cmd = string(gnuplotPath) + " -p plot_predict_no2.gp";
    if (system(cmd.c_str()) != 0) {
        cerr<<"Warning: gagal menjalankan gnuplot. "
               "Pastikan path sudah benar dan file plot_predict_no2.gp ada di sini.\n";
    }

    cout<<"Program selesai.\n";
    return 0;
}
