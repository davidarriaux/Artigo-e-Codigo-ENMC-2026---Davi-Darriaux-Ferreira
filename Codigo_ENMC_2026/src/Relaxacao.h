// Relaxacao.h
// Estrategias para escolher o fator de relaxacao w do SOR-PNBI.
// Todas devolvem um vetor com um w por regiao.

#ifndef RELAXACAO_H
#define RELAXACAO_H

#include <vector>

#include "Problema.h"

// Metodo 3: w unico que minimiza o maior modulo da diagonal de T (N = 2, bissecao)
std::vector<double> wMaiorDiagonal(const Problema& p);

// Metodo 4: w unico que minimiza a variancia da diagonal de T (N = 2, bissecao)
std::vector<double> wVariancia(const Problema& p);

// Metodo 5: um w por regiao que minimiza o maior modulo da diagonal (N = 2, bissecao)
std::vector<double> wMaiorDiagonalPorRegiao(const Problema& p);

// Metodo 6: um w por regiao que minimiza a variancia da diagonal (N = 2, bissecao)
std::vector<double> wVarianciaPorRegiao(const Problema& p);

// Metodo 0: w unico que minimiza o raio espectral da matriz T (caso otimo).
// Cria matrizes densas de tamanho N*(ntn+1): e LENTO para problemas grandes.
double wOtimoRaioEspectral(const Problema& p);

#endif
