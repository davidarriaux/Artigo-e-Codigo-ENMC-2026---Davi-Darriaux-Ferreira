// NBI.cpp

#include "NBI.h"

#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;


NBI::NBI(const Problema& p) : p_(p) {
    const int N = p.N, odq = p.odq, nI = p.nI;

    // ---- Por zona material: autovalores/autovetores e fluxo particular ----
    // Problema de autovalores:  A v = lambda v,  com  A = diag(1/mu) * M
    // Fluxo particular:         M psi = Q/2,     com  M = sigma_t*I - (sigma_s/2) * 1 * peso^T
    std::vector<VectorXd> lambda(p.nZM);
    std::vector<MatrixXd> V(p.nZM);
    std::vector<VectorXd> psiPart(p.nZM);

    for (int z = 0; z < p.nZM; z++) {
        MatrixXd M(N, N);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                M(i, j) = (i == j ? p.sigmaT[z] : 0.0) - 0.5 * p.sigmaS[z] * p.peso[j];

        MatrixXd A = M;
        for (int i = 0; i < N; i++)
            A.row(i) /= p.mu[i];

        Eigen::EigenSolver<MatrixXd> solver(A);
        lambda[z] = solver.eigenvalues().real(); // Os autovalores deste problema sao reais
        V[z] = solver.eigenvectors().real();

        psiPart[z] = M.partialPivLu().solve(VectorXd::Constant(N, 0.5 * p.fonte[z]));
    }

    // ---- Por regiao: matrizes R e S ----
    // R relaciona os fluxos de saida aos de entrada da regiao; S e a contribuicao da fonte.
    R_.resize(nI);
    Sbase_.resize(nI);

    for (int r = 0; r < nI; r++) {
        const int z = p.mapaZona[r];
        const double h = p.espRegiao[r]; // A regiao ocupa [0, h]; so a espessura importa

        MatrixXd A(N, N), B(N, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                const double lam = lambda[z](j);
                const double v = V[z](i, j);
                const double gama = (lam > 0) ? 0.0 : h; // Referencia da exponencial (evita overflow)
                const double eIni = exp(lam * gama);         // exp(-lam * (0 - gama))
                const double eFim = exp(-lam * (h - gama));

                if (i < odq) {
                    A(i, j) = v * eIni;
                    B(i, j) = v * eFim;
                } else {
                    A(i, j) = v * eFim;
                    B(i, j) = v * eIni;
                }
            }
        }

        R_[r] = B * A.inverse();
        Sbase_[r] = (MatrixXd::Identity(N, N) - R_[r]) * psiPart[z];
    }

    // ---- Espaco de trabalho ----
    Gmais_.assign(nI, MatrixXd::Zero(odq, odq));
    Gmenos_.assign(nI, MatrixXd::Zero(odq, odq));
    S_.assign(nI, VectorXd::Zero(N));
    umMenosW_.assign(nI, 0.0);

    psiD_ = MatrixXd::Zero(odq, p.ntn + 1);
    psiE_ = MatrixXd::Zero(odq, p.ntn + 1);
    phi_.assign(p.ntn + 1, 0.0);

    aux_ = VectorXd::Zero(odq);
    entrada_ = VectorXd::Zero(N);
    saida_ = VectorXd::Zero(N);
}


void NBI::aplicaRelaxacao(const std::vector<double>& w) {
    const int odq = p_.odq;

    for (int r = 0; r < p_.nI; r++) {
        Gmais_[r] = w[r] * R_[r].topLeftCorner(odq, odq);
        Gmenos_[r] = w[r] * R_[r].topRightCorner(odq, odq);
        S_[r] = w[r] * Sbase_[r];
        umMenosW_[r] = 1.0 - w[r];
    }
}


void NBI::iniciaFluxos() {
    psiD_.setZero();
    psiE_.setZero();
    psiD_.col(0).setConstant(p_.ccEsq);        // Entrada pelo contorno esquerdo
    psiE_.col(p_.ntn).setConstant(p_.ccDir);   // Entrada pelo contorno direito

    phi_.assign(p_.ntn + 1, 0.0);
    divergiu_ = false;
}


// Full NBI: em cada nodo, os fluxos de saida (nos dois sentidos) saem de R * entradas + S.
void NBI::varreduraFullNBI() {
    const int odq = p_.odq;

    for (int k = 0; k < p_.ntn; k++) {
        const int r = p_.regiaoDoNodo[k];

        entrada_.head(odq) = psiD_.col(k);
        entrada_.tail(odq) = psiE_.col(k + 1);
        saida_.noalias() = R_[r] * entrada_;

        psiD_.col(k + 1) = saida_.head(odq) + Sbase_[r].head(odq);
        psiE_.col(k) = saida_.tail(odq) + Sbase_[r].tail(odq);
    }
}


// Partial NBI com sobre-relaxacao (SOR): uma varredura para a direita (fluxos positivos)
// seguida de uma para a esquerda (fluxos negativos), ja usando os valores recem-atualizados.
//   psi_novo = G+ * psi_in1 + G- * psi_in2 + S + (1 - w) * psi_antigo
void NBI::varreduraPNBI() {
    const int odq = p_.odq;

    for (int k = 0; k < p_.ntn; k++) { // Para a direita
        const int r = p_.regiaoDoNodo[k];

        aux_.noalias() = Gmais_[r] * psiD_.col(k);
        aux_.noalias() += Gmenos_[r] * psiE_.col(k + 1);
        psiD_.col(k + 1) = aux_ + S_[r].head(odq) + umMenosW_[r] * psiD_.col(k + 1);
    }

    for (int k = p_.ntn - 1; k >= 0; k--) { // Para a esquerda
        const int r = p_.regiaoDoNodo[k];

        aux_.noalias() = Gmais_[r] * psiE_.col(k + 1);
        aux_.noalias() += Gmenos_[r] * psiD_.col(k);
        psiE_.col(k) = aux_ + S_[r].tail(odq) + umMenosW_[r] * psiE_.col(k);
    }
}


// Fluxo escalar phi = soma(peso * psi). Convergiu quando o maior desvio relativo entre duas
// iteracoes consecutivas (ou absoluto, se phi for quase zero) fica <= erro.
bool NBI::convergiu() {
    const int odq = p_.odq;
    double maiorDesvio = 0.0;

    for (int i = 0; i <= p_.ntn; i++) {
        double phi = 0.0;
        for (int j = 0; j < odq; j++)
            phi += psiD_(j, i) * p_.peso[j] + psiE_(j, i) * p_.peso[j + odq];

        if (!std::isfinite(phi)) {
            divergiu_ = true;
            return false;
        }

        const double dif = fabs(phi - phi_[i]);
        const double desvio = (fabs(phi) > 1e-8) ? dif / fabs(phi) : dif;
        if (desvio > maiorDesvio)
            maiorDesvio = desvio;

        phi_[i] = phi;
    }

    return maiorDesvio <= p_.erro;
}
