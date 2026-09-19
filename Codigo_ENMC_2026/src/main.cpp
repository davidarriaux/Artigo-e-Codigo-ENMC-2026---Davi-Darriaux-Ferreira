// main.cpp
// Programa principal: le o problema, pergunta o metodo e imprime/salva os resultados.
//
// Uso interativo:   ./sor_pnbi                       (o programa faz as perguntas)
// Uso direto:       ./sor_pnbi <arquivo> <metodo> [execucoes] [w]
//   Ex.:            ./sor_pnbi Problemas/Prob_Ex_2.txt 4 5

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Metodos.h"
#include "Problema.h"
#include "Relaxacao.h"

using namespace std;

const int METODO_TODOS = 8; // Roda os metodos 1 a 6 e mostra uma tabela comparativa


// ----------------------------------------------------------------------------
// Entrada de dados
// ----------------------------------------------------------------------------

// Le uma linha do teclado, removendo espacos, aspas (comuns ao arrastar um arquivo para o terminal)
// e marca BOM. Devolve "padrao" se a linha ficar vazia.
static string leLinha(const string& pergunta, const string& padrao) {
    cout << pergunta << " [Enter = " << padrao << "]: ";
    string s;
    getline(cin, s);

    if (s.compare(0, 3, "\xEF\xBB\xBF") == 0)
        s.erase(0, 3);
    const string lixo = " \t\r\n\"'";
    size_t ini = s.find_first_not_of(lixo);
    size_t fim = s.find_last_not_of(lixo);
    s = (ini == string::npos) ? "" : s.substr(ini, fim - ini + 1);

    return s.empty() ? padrao : s;
}

static string perguntaTexto(const string& pergunta, const string& padrao) {
    return leLinha(pergunta, padrao);
}

static double perguntaNumero(const string& pergunta, double padrao) {
    ostringstream ss;
    ss << padrao;
    try {
        return stod(leLinha(pergunta, ss.str()));
    } catch (const exception&) {
        throw runtime_error("Valor numerico invalido.");
    }
}

static void mostraMenu() {
    cout << endl << "Metodos disponiveis:" << endl;
    for (int m = 0; m <= 7; m++)
        cout << "  " << m << " - " << nomeMetodo(m) << endl;
    cout << "  " << METODO_TODOS << " - Todos os metodos de 1 a 6 (comparacao)" << endl << endl;
}

// Tenta abrir o arquivo como digitado e, se nao achar, a partir de pastas acima
// (assim funciona mesmo quando a IDE executa o programa dentro da pasta de build).
static string localizaArquivo(const string& caminho) {
    vector<string> candidatos = {caminho, "../" + caminho, "../../" + caminho};
#ifdef RAIZ_PROJETO
    candidatos.push_back(string(RAIZ_PROJETO) + "/" + caminho);
#endif
    for (const string& c : candidatos) {
        ifstream f(c);
        if (f.is_open())
            return c;
    }
    return caminho; // lerProblema mostrara a mensagem de erro
}


// ----------------------------------------------------------------------------
// Saida de dados
// ----------------------------------------------------------------------------

static string listaW(const vector<double>& w) {
    ostringstream ss;
    ss << setprecision(6) << fixed;

    bool todosIguais = true;
    for (double x : w)
        if (x != w[0]) todosIguais = false;

    if (todosIguais) {
        ss << "w = " << w[0];
    } else {
        for (size_t r = 0; r < w.size(); r++)
            ss << (r ? " / " : "") << "w_R" << r + 1 << " = " << w[r];
    }
    return ss.str();
}

// Mostra todos os fluxos se forem poucos; senao os 3 primeiros e os 3 ultimos
static void imprimeFluxo(const vector<double>& phi) {
    cout << "Fluxo escalar: |  " << scientific << setprecision(4);
    if (phi.size() <= 7) {
        for (double x : phi) cout << x << "  ";
    } else {
        for (int i = 0; i < 3; i++) cout << phi[i] << "  ";
        cout << "...  ";
        for (size_t i = phi.size() - 3; i < phi.size(); i++) cout << phi[i] << "  ";
    }
    cout << "|" << endl;
}

static void imprimeResultado(int metodo, const Resultado& r) {
    cout << endl << "--- " << nomeMetodo(metodo) << " ---" << endl;
    if (r.divergiu)
        cout << "ATENCAO: o metodo divergiu (o fluxo virou infinito ou NaN)." << endl;
    else if (!r.convergiu)
        cout << "ATENCAO: limite de iteracoes atingido sem convergir." << endl;

    cout << "Iteracoes: " << r.iteracoes << " | Tempo medio: " << scientific << setprecision(3) << r.tempo
         << " s | " << listaW(r.w) << endl;
    imprimeFluxo(r.fluxoEscalar);
}

// Salva "x  fluxo escalar" de cada aresta de nodo
static string salvaFluxo(const Problema& p, const string& arquivoProblema, int metodo, const Resultado& r) {
    size_t ini = arquivoProblema.find_last_of("/\\");
    string nome = arquivoProblema.substr(ini == string::npos ? 0 : ini + 1);
    nome = nome.substr(0, nome.find_last_of('.'));
    string saida = "resultado_" + nome + "_metodo" + to_string(metodo) + ".txt";

    ofstream f(saida);
    if (!f.is_open())
        return "";

    f << "# " << nomeMetodo(metodo) << " | " << listaW(r.w) << " | iteracoes: " << r.iteracoes << endl;
    f << "# x(cm)    fluxo_escalar" << endl;
    f << scientific << setprecision(10);

    double x = 0;
    f << x << "  " << r.fluxoEscalar[0] << endl;
    for (int k = 0; k < p.ntn; k++) {
        x += p.espRegiao[p.regiaoDoNodo[k]];
        f << x << "  " << r.fluxoEscalar[k + 1] << endl;
    }
    return saida;
}


// ----------------------------------------------------------------------------
// Programa
// ----------------------------------------------------------------------------

int main(int argc, char** argv) {
    const bool interativo = (argc < 2);
    int retorno = 0;

    try {
        cout << "=== SOR-PNBI: transporte de particulas em 1D (Full NBI, Partial NBI e SOR-PNBI) ===" << endl << endl;

        string arquivo = (argc > 1) ? argv[1] : perguntaTexto("Arquivo do problema", "Problemas/Prob_Ex_2.txt");
        arquivo = localizaArquivo(arquivo);
        Problema p = lerProblema(arquivo);

        cout << endl << "Problema lido: N = " << p.N << " | regioes = " << p.nI << " | nodos = " << p.ntn << endl;

        int metodo;
        if (argc > 2) {
            metodo = stoi(argv[2]);
        } else {
            mostraMenu();
            metodo = (int)perguntaNumero("Escolha o metodo", 4);
        }

        int execucoes = (argc > 3) ? stoi(argv[3]) : (interativo ? (int)perguntaNumero("Numero de execucoes (para media do tempo)", 5) : 5);

        double wFixo = 1.0;
        if (metodo == 7)
            wFixo = (argc > 4) ? stod(argv[4]) : perguntaNumero("Valor de w", 1.1);

        if (metodo == METODO_TODOS) {
            struct Linha { int metodo; Resultado r; };
            vector<Linha> tabela;
            for (int m = 1; m <= 6; m++) {
                Resultado r = executaMetodo(p, m, execucoes);
                imprimeResultado(m, r);
                tabela.push_back({m, r});
            }

            cout << endl << "=== Resumo ===" << endl;
            cout << left << setw(42) << "Metodo" << right << setw(12) << "Iteracoes" << setw(14) << "Tempo (s)" << endl;
            for (const Linha& l : tabela)
                cout << left << setw(42) << nomeMetodo(l.metodo) << right << setw(12) << l.r.iteracoes
                     << setw(14) << scientific << setprecision(3) << l.r.tempo << endl;
        } else {
            if (metodo == 0) {
                cout << endl << "Procurando o w otimo (menor raio espectral). Pode demorar em problemas grandes..." << endl;
                wFixo = wOtimoRaioEspectral(p);
                cout << "w otimo encontrado: " << fixed << setprecision(5) << wFixo << endl;
            }

            Resultado r = executaMetodo(p, metodo, execucoes, wFixo);
            imprimeResultado(metodo, r);

            string saida = salvaFluxo(p, arquivo, metodo, r);
            if (!saida.empty())
                cout << endl << "Fluxo escalar completo salvo em: " << saida << endl;
        }
    } catch (const exception& e) {
        cerr << endl << "ERRO: " << e.what() << endl;
        retorno = 1;
    }

    if (interativo) {
        cout << endl << "Pressione Enter para sair...";
        string s;
        getline(cin, s);
    }
    return retorno;
}
