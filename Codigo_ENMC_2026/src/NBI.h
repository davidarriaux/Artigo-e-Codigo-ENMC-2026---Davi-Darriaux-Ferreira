// NBI.h
// Nucleo numerico: matrizes resposta de cada regiao, varreduras (Full NBI e Partial NBI / SOR-PNBI)
// e criterio de convergencia do fluxo escalar.

#ifndef NBI_H
#define NBI_H

#include <vector>

#include "Eigen/Dense"
#include "Problema.h"

class NBI {
public:
    // Calcula autovalores/autovetores, fluxo particular e as matrizes R e S de cada regiao.
    explicit NBI(const Problema& problema);

    // Multiplica G+, G- e S pelo fator de relaxacao w de cada regiao (w = 1 -> Partial NBI puro).
    void aplicaRelaxacao(const std::vector<double>& w);

    void iniciaFluxos();     // Zera os fluxos e aplica as condicoes de contorno
    void varreduraFullNBI(); // Uma iteracao do Full NBI (usa R e S, sem relaxacao)
    void varreduraPNBI();    // Uma iteracao do Partial NBI / SOR-PNBI (usa G+, G-, S e w)

    // Atualiza o fluxo escalar e informa se o desvio relativo maximo ficou abaixo do erro pedido.
    bool convergiu();
    bool divergiu() const { return divergiu_; } // true se apareceu NaN/infinito

    const std::vector<Eigen::MatrixXd>& R() const { return R_; } // Matriz R (N x N) de cada regiao
    const std::vector<double>& fluxoEscalar() const { return phi_; } // ntn + 1 valores

private:
    const Problema& p_;

    std::vector<Eigen::MatrixXd> R_;      // Matriz resposta de cada regiao
    std::vector<Eigen::VectorXd> Sbase_;  // Termo de fonte de cada regiao, sem relaxacao

    std::vector<Eigen::MatrixXd> Gmais_;  // w * R(0:odq, 0:odq)
    std::vector<Eigen::MatrixXd> Gmenos_; // w * R(0:odq, odq:N)
    std::vector<Eigen::VectorXd> S_;      // w * Sbase
    std::vector<double> umMenosW_;        // 1 - w de cada regiao

    Eigen::MatrixXd psiD_; // Fluxo angular no sentido +x (odq linhas, ntn+1 colunas)
    Eigen::MatrixXd psiE_; // Fluxo angular no sentido -x
    std::vector<double> phi_; // Fluxo escalar em cada aresta de nodo
    bool divergiu_ = false;

    Eigen::VectorXd aux_, entrada_, saida_; // Vetores de trabalho (evitam alocar dentro dos lacos)
};

#endif
