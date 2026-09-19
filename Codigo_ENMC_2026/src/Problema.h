// Problema.h
// Dados de entrada do problema (lidos de um .txt) e grandezas derivadas deles
// (quadratura de Gauss-Legendre, nodos, espessuras).

#ifndef PROBLEMA_H
#define PROBLEMA_H

#include <string>
#include <vector>

struct Problema {

    // ---- Lidos do arquivo ----
    int N = 0;                    // Ordem da quadratura (numero par de direcoes)
    int nI = 0;                   // Numero de regioes
    int nZM = 0;                  // Numero de zonas materiais
    std::vector<int> mapaZona;    // Zona material de cada regiao
    std::vector<double> tamRegiao;// Tamanho de cada regiao (cm)
    std::vector<int> nNodos;      // Numero de nodos de cada regiao
    double ccEsq = 0;             // Fluxo incidente no contorno esquerdo
    double ccDir = 0;             // Fluxo incidente no contorno direito
    std::vector<double> sigmaT;   // Secao de choque total de cada zona
    std::vector<double> sigmaS;   // Secao de choque de espalhamento de cada zona
    std::vector<double> fonte;    // Fonte externa de cada zona

    // ---- Parametros do metodo iterativo ----
    double erro = 1e-10;          // Criterio de convergencia (desvio relativo do fluxo escalar)
    int maxIter = 100000;         // Numero maximo de iteracoes

    // ---- Calculados a partir dos dados acima ----
    int odq = 0;                       // Metade da ordem da quadratura (N/2)
    int ntn = 0;                       // Numero total de nodos
    std::vector<int> regiaoDoNodo;     // Regiao a que cada nodo pertence
    std::vector<double> espRegiao;     // Espessura de cada nodo em cada regiao
    std::vector<double> mu;            // Direcoes discretas (mu[0..odq-1] > 0, o resto e o oposto)
    std::vector<double> peso;          // Pesos da quadratura

    void calculaDerivados();           // Preenche os campos "calculados"
    Problema comOrdem(int novoN) const;// Copia do problema com outra ordem de quadratura
};

// Le o arquivo do problema. Em caso de erro lanca std::runtime_error com mensagem explicativa.
Problema lerProblema(const std::string& caminho);

#endif
