// Metodos.cpp

#include "Metodos.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "NBI.h"
#include "Relaxacao.h"

using namespace std;


string nomeMetodo(int metodo) {
    switch (metodo) {
        case 0: return "SOR-PNBI (w otimo, menor raio espectral)";
        case 1: return "Full NBI";
        case 2: return "Partial NBI (w = 1)";
        case 3: return "SOR-PNBI (w pela diagonal)";
        case 4: return "SOR-PNBI (w pela variancia)";
        case 5: return "SOR-PNBI (w por regiao - diagonal)";
        case 6: return "SOR-PNBI (w por regiao - variancia)";
        case 7: return "SOR-PNBI (w informado)";
    }
    return "Metodo invalido";
}


// Uma execucao completa: prepara matrizes, escolhe w e itera ate convergir.
// O tempo medido inclui tudo isso, inclusive o calculo de w.
static Resultado executaUmaVez(const Problema& p, int metodo, double wFixo) {
    const auto inicio = chrono::steady_clock::now();

    NBI nbi(p);

    vector<double> w(p.nI, 1.0);
    switch (metodo) {
        case 0:
        case 7: w.assign(p.nI, wFixo); break;
        case 3: w = wMaiorDiagonal(p); break;
        case 4: w = wVariancia(p); break;
        case 5: w = wMaiorDiagonalPorRegiao(p); break;
        case 6: w = wVarianciaPorRegiao(p); break;
        default: break; // 1 e 2: w = 1
    }
    nbi.aplicaRelaxacao(w);
    nbi.iniciaFluxos();

    int it = 0;
    bool convergiu = false;
    while (!convergiu && !nbi.divergiu() && it < p.maxIter) {
        if (metodo == 1)
            nbi.varreduraFullNBI();
        else
            nbi.varreduraPNBI();

        convergiu = nbi.convergiu();
        it++;
    }

    Resultado res;
    res.tempo = chrono::duration<double>(chrono::steady_clock::now() - inicio).count();
    res.iteracoes = it;
    res.convergiu = convergiu;
    res.divergiu = nbi.divergiu();
    res.w = w;
    res.fluxoEscalar = nbi.fluxoEscalar();
    return res;
}


Resultado executaMetodo(const Problema& p, int metodo, int execucoes, double wFixo) {
    if (metodo < METODO_MIN || metodo > METODO_MAX)
        throw runtime_error("Metodo invalido: " + to_string(metodo));
    if (execucoes < 1)
        execucoes = 1;

    vector<double> tempos;
    Resultado res;
    ostringstream descarte;
    streambuf* saidaOriginal = cout.rdbuf();
    for (int i = 0; i < execucoes; i++) {
        if (i == 1) cout.rdbuf(descarte.rdbuf()); // Avisos iguais nao precisam se repetir a cada execucao
        res = executaUmaVez(p, metodo, wFixo);
        tempos.push_back(res.tempo);
    }
    cout.rdbuf(saidaOriginal);

    sort(tempos.begin(), tempos.end());
    double soma = 0;
    int ini = (execucoes > 2) ? 1 : 0;                // Descarta menor e maior tempo
    int fim = (execucoes > 2) ? execucoes - 1 : execucoes;
    for (int i = ini; i < fim; i++)
        soma += tempos[i];

    res.tempo = soma / (fim - ini);
    return res;
}
