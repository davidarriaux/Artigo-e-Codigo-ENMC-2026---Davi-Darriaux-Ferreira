// Metodos.h
// Executa um dos metodos iterativos ate a convergencia e mede o tempo.

#ifndef METODOS_H
#define METODOS_H

#include <string>
#include <vector>

#include "Problema.h"

// Numeracao dos metodos
//   0 = SOR-PNBI com w otimo (menor raio espectral - lento)   5 = SOR-PNBI, w por regiao (diagonal)
//   1 = Full NBI                                              6 = SOR-PNBI, w por regiao (variancia)
//   2 = Partial NBI (w = 1)                                   7 = SOR-PNBI com w informado pelo usuario
//   3 = SOR-PNBI, w pela diagonal (maior modulo)
//   4 = SOR-PNBI, w pela variancia
const int METODO_MIN = 0;
const int METODO_MAX = 7;

std::string nomeMetodo(int metodo);

struct Resultado {
    double tempo = 0;              // Tempo medio (s)
    int iteracoes = 0;
    bool convergiu = false;        // false se atingiu o maximo de iteracoes ou divergiu
    bool divergiu = false;
    std::vector<double> w;         // Fator de relaxacao de cada regiao
    std::vector<double> fluxoEscalar;
};

// Executa o metodo "execucoes" vezes. O tempo medio descarta a mais rapida e a mais lenta
// quando ha mais de 2 execucoes. "wFixo" so e usado nos metodos 0 e 7.
Resultado executaMetodo(const Problema& p, int metodo, int execucoes, double wFixo = 1.0);

#endif
