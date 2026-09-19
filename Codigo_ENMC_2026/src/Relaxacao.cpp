// Relaxacao.cpp

#include "Relaxacao.h"

#include <cmath>
#include <iostream>
#include <string>

#include "NBI.h"

using namespace std;
using Eigen::MatrixXd;


namespace {

// ----------------------------------------------------------------------------
// Diagonal de T para N = 2
// ----------------------------------------------------------------------------
// Com N = 2, G+ e G- sao escalares em cada nodo (gma[k], gme[k]). O elemento i da diagonal de T e
//     diag[i] = (1 - w) + sum_{j>i} [ prod_{k=i+1}^{j-1} (w gma[k])^2 ] * w^2 * gme[i] * gme[j]
// Os somatorios sao calculados de tras para frente por recorrencia (custo O(n), nao O(n^2)):
//     T[i] = gme[i+1] + (w gma[i+1])^2 * T[i+1]
//     A[i] = gme[i+1] + (w gma[i+1])^2 * (A[i+1] + T[i+1])      (A e usado na derivada)
// Se "deriv" nao for nulo, tambem devolve d(diag[i])/dw.
void diagonalN2(const vector<double>& gma, const vector<double>& gme, double w,
                vector<double>& diag, vector<double>* deriv) {
    const int n = gma.size();
    const double w2 = w * w;

    diag[n - 1] = 1 - w;
    if (deriv) (*deriv)[n - 1] = -1;

    double T = 0, A = 0;
    for (int i = n - 2; i >= 0; i--) {
        const double f = w2 * gma[i + 1] * gma[i + 1];
        const double Tn = gme[i + 1] + f * T;
        const double An = gme[i + 1] + f * (A + T);
        T = Tn;
        A = An;

        diag[i] = 1 - w + w2 * gme[i] * T;
        if (deriv) (*deriv)[i] = -1 + (2.0 / w) * w2 * gme[i] * A;
    }
}


// G+ e G- (escalares, N = 2) de cada regiao
struct GN2 {
    vector<double> mais, menos;
};

GN2 calculaGN2(const Problema& p) {
    Problema p2 = p.comOrdem(2);
    NBI nbi2(p2);

    GN2 g;
    for (int r = 0; r < p.nI; r++) {
        g.mais.push_back(nbi2.R()[r](0, 0));
        g.menos.push_back(nbi2.R()[r](0, 1));
    }
    return g;
}

// Valores de G+ e G- de cada nodo do problema
void gPorNodo(const Problema& p, const GN2& g, vector<double>& gma, vector<double>& gme) {
    gma.resize(p.ntn);
    gme.resize(p.ntn);
    for (int k = 0; k < p.ntn; k++) {
        gma[k] = g.mais[p.regiaoDoNodo[k]];
        gme[k] = g.menos[p.regiaoDoNodo[k]];
    }
}


// ----------------------------------------------------------------------------
// Funcoes auxiliares
// ----------------------------------------------------------------------------

int indiceMaiorModulo(const vector<double>& v) {
    int k = 0;
    for (int i = 1; i < (int)v.size(); i++)
        if (fabs(v[i]) > fabs(v[k]))
            k = i;
    return k;
}

double mediaModulos(const vector<double>& v) {
    double s = 0;
    for (double x : v) s += fabs(x);
    return s / v.size();
}

double varianciaModulos(const vector<double>& v) {
    const double m = mediaModulos(v);
    double s = 0;
    for (double x : v) s += (fabs(x) - m) * (fabs(x) - m);
    return s / v.size();
}


// Procura w em (1, 2) que minimiza uma funcao f(w), a partir da sua derivada:
//  1) avanca w de "passo" em "passo" ate a derivada ficar positiva (o minimo foi ultrapassado);
//  2) refina por bissecao no intervalo [w - passo, w] ate o desvio relativo ser < 1e-4.
// Se nao houver minimo em (1, 2), avisa e devolve w = 1.
template <class Derivada>
double bissecao(Derivada derivada, double passo, const string& onde) {
    double w0 = 1;
    double d;
    do {
        w0 += passo;
        d = derivada(w0);
    } while (d <= 0 && w0 < 2.0);

    if (d <= 0) {
        cout << "AVISO: nao foi encontrado minimo em 1 < w < 2 (" << onde << "). Usando w = 1." << endl;
        return 1.0;
    }

    double a = w0 - passo, b = w0, wc = 0, desvio = 1;
    for (int it = 0; desvio > 1e-4 && it < 1000; it++) {
        const double wAnterior = wc;
        wc = (a + b) / 2;

        if (derivada(wc) > 0)
            b = wc;
        else
            a = wc;

        desvio = fabs((wc - wAnterior) / wc);
    }
    return wc;
}

const double H = 1e-4; // Passo das derivadas por diferencas finitas

} // namespace


// ----------------------------------------------------------------------------
// Estrategias
// ----------------------------------------------------------------------------

vector<double> wMaiorDiagonal(const Problema& p) {
    vector<double> gma, gme;
    gPorNodo(p, calculaGN2(p), gma, gme);

    vector<double> diag(p.ntn), deriv(p.ntn);
    auto derivada = [&](double w) { // Derivada do maior modulo da diagonal
        diagonalN2(gma, gme, w, diag, &deriv);
        const int k = indiceMaiorModulo(diag);
        return diag[k] > 0 ? deriv[k] : -deriv[k];
    };

    return vector<double>(p.nI, bissecao(derivada, 0.01, "diagonal global"));
}


vector<double> wVariancia(const Problema& p) {
    vector<double> gma, gme;
    gPorNodo(p, calculaGN2(p), gma, gme);

    // A diagonal de T tem 2(ntn+1) termos: os ntn primeiros dependem de G+ e G-, os demais sao 1 - w
    const int tam = 2 * (p.ntn + 1);
    vector<double> diag(p.ntn), deriv(p.ntn), d(tam), dd(tam);

    auto derivada = [&](double w) { // Derivada da variancia dos modulos da diagonal
        diagonalN2(gma, gme, w, diag, &deriv);
        for (int i = 0; i < tam; i++) {
            d[i] = (i < p.ntn) ? diag[i] : 1 - w;
            dd[i] = (i < p.ntn) ? deriv[i] : -1;
        }

        const double media = mediaModulos(d);
        double soma = 0;
        for (int i = 0; i < tam; i++)
            soma += (fabs(d[i]) - media) * (d[i] > 0 ? dd[i] : -dd[i]);
        return soma / (p.ntn + 1);
    };

    return vector<double>(p.nI, bissecao(derivada, 0.01, "variancia global"));
}


vector<double> wMaiorDiagonalPorRegiao(const Problema& p) {
    const GN2 g = calculaGN2(p);
    vector<double> w(p.nI);

    for (int r = 0; r < p.nI; r++) {
        const int n = p.nNodos[r];
        const vector<double> gma(n, g.mais[r]), gme(n, g.menos[r]);
        vector<double> diag(n);

        auto f = [&](double x) { // Maior modulo da diagonal da regiao
            diagonalN2(gma, gme, x, diag, nullptr);
            return fabs(diag[indiceMaiorModulo(diag)]);
        };
        auto derivada = [&](double x) { return (f(x + H) - f(x)) / H; };

        w[r] = bissecao(derivada, 0.1, "regiao " + to_string(r + 1));
    }
    return w;
}


vector<double> wVarianciaPorRegiao(const Problema& p) {
    const GN2 g = calculaGN2(p);
    vector<double> w(p.nI);

    for (int r = 0; r < p.nI; r++) {
        const int n = p.nNodos[r];

        if (n == 1) { // Todos os termos da diagonal sao 1 - w: variancia nula, nao ha o que otimizar
            w[r] = 1.0;
            continue;
        }

        const vector<double> gma(n, g.mais[r]), gme(n, g.menos[r]);
        vector<double> diag(n), d(2 * (n + 1));

        auto f = [&](double x) { // Variancia da diagonal da regiao (2(n+1) termos)
            diagonalN2(gma, gme, x, diag, nullptr);
            for (int i = 0; i < (int)d.size(); i++)
                d[i] = (i < n) ? diag[i] : 1 - x;
            return varianciaModulos(d);
        };
        auto derivada = [&](double x) { return (f(x + H) - f(x)) / H; };

        w[r] = bissecao(derivada, 0.1, "regiao " + to_string(r + 1));
    }
    return w;
}


// ----------------------------------------------------------------------------
// w otimo: menor raio espectral da matriz de iteracao T
// ----------------------------------------------------------------------------

namespace {

// Matriz de iteracao do SOR-PNBI:  T = (I - C)^-1 D
// C: parte "ja atualizada" (varredura em curso) | D: parte da iteracao anterior (inclui 1 - w)
MatrixXd matrizT(const Problema& p, const vector<MatrixXd>& R, double w) {
    const int N = p.N, odq = p.odq, ntn = p.ntn;
    const int tam = N * (ntn + 1);

    MatrixXd C = MatrixXd::Zero(tam, tam);
    MatrixXd D = MatrixXd::Zero(tam, tam);
    D.diagonal().setConstant(1 - w);

    for (int l = 1; l <= ntn; l++) { // Fluxos positivos
        const MatrixXd& Rr = R[p.regiaoDoNodo[l - 1]];
        C.block(odq * l, odq * (l - 1), odq, odq) = w * Rr.topLeftCorner(odq, odq);
        D.block(odq * l, tam - odq * (l + 1), odq, odq) = w * Rr.topRightCorner(odq, odq);
    }
    for (int l = 0; l < ntn; l++) { // Fluxos negativos
        const MatrixXd& Rr = R[p.regiaoDoNodo[l]];
        C.block(tam - odq * (l + 1), tam - odq * (l + 2), odq, odq) = w * Rr.topLeftCorner(odq, odq);
        C.block(tam - odq * (l + 1), odq * l, odq, odq) = w * Rr.topRightCorner(odq, odq);
    }

    // I - C e triangular inferior, entao basta resolver por substituicao
    const MatrixXd A = MatrixXd::Identity(tam, tam) - C;
    return A.triangularView<Eigen::Lower>().solve(D);
}

double raioEspectral(const MatrixXd& T) {
    Eigen::EigenSolver<MatrixXd> solver(T, false); // So autovalores
    return solver.eigenvalues().cwiseAbs().maxCoeff();
}

} // namespace


double wOtimoRaioEspectral(const Problema& p) {
    const NBI nbi(p); // R nao depende de w
    auto raio = [&](double w) { return raioEspectral(matrizT(p, nbi.R(), w)); };

    // 1) Varredura grossa: avanca de 0.01 em 0.01 enquanto o raio espectral diminui
    const double passo = 0.01;
    double melhor = 1.0, fMelhor = raio(melhor);
    while (melhor + passo < 2.0) {
        const double f = raio(melhor + passo);
        if (f >= fMelhor)
            break;
        melhor += passo;
        fMelhor = f;
    }

    // 2) Refino por busca da secao aurea em [melhor - passo, melhor + passo]
    double a = max(1.0, melhor - passo), b = min(2.0, melhor + passo);
    const double phi = (sqrt(5.0) - 1) / 2;
    double w1 = b - phi * (b - a), w2 = a + phi * (b - a);
    double f1 = raio(w1), f2 = raio(w2);

    while (b - a > 1e-5) {
        if (f1 < f2) {
            b = w2; w2 = w1; f2 = f1;
            w1 = b - phi * (b - a);
            f1 = raio(w1);
        } else {
            a = w1; w1 = w2; f1 = f2;
            w2 = a + phi * (b - a);
            f2 = raio(w2);
        }
    }
    return (a + b) / 2;
}
