#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>  

using namespace std;

void risolviPALU2x2(const double A[2][2], const double b[2], double x[2]) {
    double LU[2][2] = {
        { A[0][0], A[0][1] },
        { A[1][0], A[1][1] }
    };

    // Vettore di permutazione per tenere traccia degli swap
    int p[2] = { 0, 1 };

    // Controllo per fare pivoting
    if (fabs(LU[1][0]) > fabs(LU[0][0])) {
        swap(LU[0][0], LU[1][0]);
        swap(LU[0][1], LU[1][1]);
        swap(p[0], p[1]);
    }

    double L21 = LU[1][0] / LU[0][0];
    // Salvo L21 in posizione (1,0)
    LU[1][0] = L21;
    // Aggiorno l'elemento in U(1,1)
    LU[1][1] = LU[1][1] - L21 * LU[0][1];             

    // Costruisco Pb 
    double bp[2] = { b[p[0]], b[p[1]] };

    double y[2];
    y[0] = bp[0];
    y[1] = bp[1] - L21 * y[0];

    // Risolviamo Ux = y  con sostituzione all'indietro.
    x[1] = y[1] / LU[1][1];
    x[0] = (y[0] - LU[0][1] * x[1]) / LU[0][0];
}

// Per la decomposizione QR sfrutto il metodo di Gram-Schmidt.
void risolviQR2x2(const double A[2][2], const double b[2], double x[2]) {
    // a0 è la prima colonna e a1 la seconda.
    double a0[2] = { A[0][0], A[1][0] };
    double a1[2] = { A[0][1], A[1][1] };

    double norm_a0 = sqrt(a0[0]*a0[0] + a0[1]*a0[1]);
    // Check per eventuale prima colonna nulla
    if (norm_a0 == 0) {
        cerr << "Errore: Prima colonna nulla nella decomposizione QR." << endl;
        return;
    }
    double q0[2] = { a0[0] / norm_a0, a0[1] / norm_a0 };
    double r0 = norm_a0;

    // r1 è la proiezione di a1 su q0.
    double r1 = q0[0] * a1[0] + q0[1] * a1[1];

    double u[2] = { a1[0] - r1 * q0[0], a1[1] - r1 * q0[1] };
    double r2 = sqrt(u[0]*u[0] + u[1]*u[1]);

    // Normalizzo u per ottenere q1
    double q1[2];
    // Utilizzo questo check per evitare di dividere per zero o quasi per zero
    if (r2 > 1e-12) {
        q1[0] = u[0] / r2;
        q1[1] = u[1] / r2;
    } else {
        // Se r2 è quasi zero, a tutti gli effetti la seconda colonna è linearmente dipendente.
        q1[0] = 0.0; 
        q1[1] = 0.0;
    }

    double y0 = q0[0] * b[0] + q0[1] * b[1];
    double y1 = q1[0] * b[0] + q1[1] * b[1];

    // Risolvo Rx = y con sostituzione all'indietro
    double x1_val;
    if (fabs(r2) > 1e-12) {
        x1_val = y1 / r2;
    } else {
        x1_val = 0.0;
    }
    
    double x0_val = (y0 - r1 * x1_val) / r0;
    x[0] = x0_val;
    x[1] = x1_val;
}

// Funzione per il calcolo dell'errore relativo
double erroreRelativo(const double x[2]) {
    double diff0 = x[0] + 1.0;
    double diff1 = x[1] + 1.0;
    double err = sqrt(diff0 * diff0 + diff1 * diff1);
    double normAtteso = sqrt(1.0 * 1.0 + 1.0 * 1.0); // sqrt(2)
    return err / normAtteso;
}

int main() {
    cout << fixed << setprecision(15);
    
    // Sistema 1:
    double A1[2][2] = {
        { 5.547001962252291e-01, -3.770900990025203e-02 },
        { 8.320502943378437e-01, -9.992887623566787e-01 }
    };
    double b1[2] = { -5.169911863249772e-01, 1.672384680188350e-01 };
    double x1_PALU[2], x1_QR[2];
    
    risolviPALU2x2(A1, b1, x1_PALU);
    risolviQR2x2(A1, b1, x1_QR);
    
    cout << "Sistema 1:" << endl;
    cout << "Soluzione PALU: (" << x1_PALU[0] << ", " << x1_PALU[1] << ")" << endl;
    cout << "Errore relativo PALU: " << erroreRelativo(x1_PALU) << endl;
    cout << "Soluzione QR  : (" << x1_QR[0] << ", " << x1_QR[1] << ")" << endl;
    cout << "Errore relativo QR  : " << erroreRelativo(x1_QR) << endl;
    cout << endl;
    
    // Sistema 2:
    double A2[2][2] = {
        { 5.547001962252291e-01, -5.540607316466765e-01 },
        { 8.320502943378437e-01, -8.324762492991313e-01 }
    };
    double b2[2] = { -6.394645785530173e-04, 4.259549612877223e-04 };
    double x2_PALU[2], x2_QR[2];
    
    risolviPALU2x2(A2, b2, x2_PALU);
    risolviQR2x2(A2, b2, x2_QR);
    
    cout << "Sistema 2:" << endl;
    cout << "Soluzione PALU: (" << x2_PALU[0] << ", " << x2_PALU[1] << ")" << endl;
    cout << "Errore relativo PALU: " << erroreRelativo(x2_PALU) << endl;
    cout << "Soluzione QR  : (" << x2_QR[0] << ", " << x2_QR[1] << ")" << endl;
    cout << "Errore relativo QR  : " << erroreRelativo(x2_QR) << endl;
    cout << endl;
    
    // Sistema 3:
    double A3[2][2] = {
        { 5.547001962252291e-01, -5.547001955851905e-01 },
        { 8.320502943378437e-01, -8.320502947645361e-01 }
    };
    double b3[2] = { -6.400391328043042e-10, 4.266924591433963e-10 };
    double x3_PALU[2], x3_QR[2];
    
    risolviPALU2x2(A3, b3, x3_PALU);
    risolviQR2x2(A3, b3, x3_QR);
    
    cout << "Sistema 3:" << endl;
    cout << "Soluzione PALU: (" << x3_PALU[0] << ", " << x3_PALU[1] << ")" << endl;
    cout << "Errore relativo PALU: " << erroreRelativo(x3_PALU) << endl;
    cout << "Soluzione QR  : (" << x3_QR[0] << ", " << x3_QR[1] << ")" << endl;
    cout << "Errore relativo QR  : " << erroreRelativo(x3_QR) << endl;
    cout << endl;
    
    return 0;
}
