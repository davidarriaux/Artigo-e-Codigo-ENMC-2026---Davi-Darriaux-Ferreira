// Problema.cpp

#include "Problema.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "Eigen/Dense"

using namespace std;


// Nos (mu) e pesos da quadratura de Gauss-Legendre de ordem N (metodo de Golub-Welsch):
// os nos sao os autovalores da matriz de Jacobi e os pesos vem do 1o componente dos autovetores.
static void quadraturaGaussLegendre(int N, vector<double>& mu, vector<double>& peso) {
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(N, N);
    for (int j = 1; j < N; j++)
        J(j - 1, j) = J(j, j - 1) = sqrt(1.0 / (4.0 - 1.0 / (double(j) * j)));

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J); // Autovalores em ordem crescente

    int odq = N / 2;
    mu.assign(N, 0.0);
    peso.assign(N, 0.0);
    for (int i = 0; i < odq; i++) {
        mu[i] = -es.eigenvalues()(i);          // Direcoes positivas
        mu[i + odq] = -mu[i];                  // Direcoes negativas
        peso[i] = 2.0 * pow(es.eigenvectors()(0, i), 2);
        peso[i + odq] = peso[i];
    }
}


void Problema::calculaDerivados() {
    odq = N / 2;

    ntn = 0;
    for (int n : nNodos)
        ntn += n;

    regiaoDoNodo.clear();
    for (int r = 0; r < nI; r++)
        regiaoDoNodo.insert(regiaoDoNodo.end(), nNodos[r], r);

    espRegiao.assign(nI, 0.0);
    for (int r = 0; r < nI; r++)
        espRegiao[r] = tamRegiao[r] / nNodos[r];

    quadraturaGaussLegendre(N, mu, peso);
}


Problema Problema::comOrdem(int novoN) const {
    Problema p = *this;
    p.N = novoN;
    p.calculaDerivados();
    return p;
}


// ----------------------------------------------------------------------------
// Leitura do arquivo
// ----------------------------------------------------------------------------

// O arquivo tem 10 pares de linhas: uma linha de titulo e, logo abaixo, a linha com os valores.
// Linhas em branco sao ignoradas.

static vector<string> linhasNaoVazias(const string& caminho) {
    ifstream arq(caminho);
    if (!arq.is_open())
        throw runtime_error("Nao foi possivel abrir o arquivo '" + caminho + "'. Confira o nome/caminho "
                            "(ex.: Problemas/Prob_Ex_2.txt) e se o programa esta sendo executado a partir da pasta do projeto.");

    vector<string> linhas;
    string linha;
    while (getline(arq, linha))
        if (linha.find_first_not_of(" \t\r\n") != string::npos)
            linhas.push_back(linha);
    return linhas;
}

template <class T>
static vector<T> leValores(const string& linha, int quantidade, const string& titulo) {
    istringstream ss(linha);
    vector<T> v(quantidade);
    for (int i = 0; i < quantidade; i++)
        if (!(ss >> v[i]))
            throw runtime_error("Na linha de dados de '" + titulo + "' eram esperados " +
                                to_string(quantidade) + " valor(es), mas a linha e: '" + linha + "'.");
    return v;
}

Problema lerProblema(const string& caminho) {
    vector<string> linhas = linhasNaoVazias(caminho);
    if (linhas.size() < 20)
        throw runtime_error("O arquivo '" + caminho + "' esta incompleto: sao esperadas 20 linhas "
                            "(10 titulos + 10 linhas de valores).");

    // Titulo do item k: linhas[2k] | valores do item k: linhas[2k+1]
    auto titulo = [&](int k) { return linhas[2 * k]; };
    auto dados = [&](int k) { return linhas[2 * k + 1]; };

    Problema p;
    p.N = leValores<int>(dados(0), 1, titulo(0))[0];
    p.nI = leValores<int>(dados(1), 1, titulo(1))[0];
    p.nZM = leValores<int>(dados(2), 1, titulo(2))[0];

    if (p.N < 2 || p.N % 2 != 0)
        throw runtime_error("A ordem da quadratura (N) deve ser um numero par >= 2.");
    if (p.nI < 1 || p.nZM < 1)
        throw runtime_error("O numero de regioes e de zonas materiais deve ser >= 1.");

    p.mapaZona = leValores<int>(dados(3), p.nI, titulo(3));
    p.tamRegiao = leValores<double>(dados(4), p.nI, titulo(4));
    p.nNodos = leValores<int>(dados(5), p.nI, titulo(5));

    vector<double> cc = leValores<double>(dados(6), 2, titulo(6));
    p.ccEsq = cc[0];
    p.ccDir = cc[1];

    p.sigmaT = leValores<double>(dados(7), p.nZM, titulo(7));
    p.sigmaS = leValores<double>(dados(8), p.nZM, titulo(8));
    p.fonte = leValores<double>(dados(9), p.nZM, titulo(9));

    for (int r = 0; r < p.nI; r++) {
        if (p.mapaZona[r] < 0 || p.mapaZona[r] >= p.nZM)
            throw runtime_error("Mapeamento das regioes: a regiao " + to_string(r + 1) +
                                " aponta para uma zona material inexistente (use valores de 0 a " +
                                to_string(p.nZM - 1) + ").");
        if (p.nNodos[r] < 1)
            throw runtime_error("Cada regiao precisa ter pelo menos 1 nodo.");
        if (p.tamRegiao[r] <= 0)
            throw runtime_error("O tamanho de cada regiao deve ser positivo.");
    }
    for (int z = 0; z < p.nZM; z++)
        if (p.sigmaT[z] <= 0)
            throw runtime_error("A secao de choque total deve ser positiva.");

    p.calculaDerivados();
    return p;
}
