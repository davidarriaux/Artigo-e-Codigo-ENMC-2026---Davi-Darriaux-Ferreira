# SOR-PNBI — Transporte de partículas em 1D (Full NBI, Partial NBI e SOR-PNBI)

Programa (pronto para usar, sem precisar compilar nada) que acompanha o artigo **"TÉCNICA DE SOBRE-RELAXAÇÃO SUCESSIVA PARA O ESQUEMA ITERATIVO DE INVERSÃO NODAL PARCIAL EM SIMULAÇÕES DE TRANSPORTE S<sub>N</sub> DE PARTÍCULAS NEUTRAS EM GEOMETRIA UNIDIMENSIONAL CARTESIANA"**.

Repositório (artigo e programa): https://github.com/davidarriaux/Artigo-e-Codigo-ENMC-2026---Davi-Darriaux-Ferreira

Ele resolve a equação de transporte de partículas (ordenadas discretas, geometria plana 1D, meio com várias regiões, fonte fixa) pelos métodos iterativos:

| Nº | Método |
|----|--------|
| 0 | SOR-PNBI com o ***w* ótimo** (menor raio espectral de T) — **lento**, só para problemas pequenos |
| 1 | **Full NBI** |
| 2 | **Partial NBI** (equivale ao SOR-PNBI com *w* = 1) |
| 3 | SOR-PNBI com **um *w* geral**, escolhido pelo maior módulo da diagonal de T |
| 4 | SOR-PNBI com **um *w* geral**, escolhido pela variância da diagonal de T |
| 5 | SOR-PNBI com **um *w* por região** (maior módulo da diagonal) |
| 6 | SOR-PNBI com **um *w* por região** (variância da diagonal) |
| 7 | SOR-PNBI com um ***w* informado por você** (para testar valores de *w* à vontade) |
| 8 | Roda os métodos 1 a 6 e mostra uma tabela comparando iterações e tempo |

**Você não precisa instalar nada, nem compilar nada.** Basta baixar esta pasta e rodar o programa já pronto (executável) do seu sistema.

---

## Sumário

1. [O que tem neste repositório](#1-o-que-tem-neste-repositório)
2. [Passo 1 — Baixar](#2-passo-1--baixar)
3. [Passo 2 — Executar](#3-passo-2--executar)
4. [Como saber se está tudo certo (valores de referência)](#4-como-saber-se-está-tudo-certo-valores-de-referência)
5. [Criando o seu próprio problema (qualquer modelo 1D de fonte fixa)](#5-criando-o-seu-próprio-problema-qualquer-modelo-1d-de-fonte-fixa)
6. [Problemas comuns](#6-problemas-comuns)
7. [Licença e créditos](#7-licença-e-créditos)

---

## 1. O que tem neste repositório

```
.
├── README.md                      <- este arquivo
├── LICENSE_Eigen_MPL2.txt         <- licença da biblioteca usada internamente (Eigen)
├── Problemas/
│   ├── Problema_Modelo.txt        <- modelo em branco, pronto para você editar (1 região, roda na hora)
│   ├── Prob_Ex_1.txt              <- Exemplo 1 do artigo (1 região, N = 32, 100 nodos)
│   └── Prob_Ex_2.txt              <- Exemplo 2 do artigo (4 regiões, N = 16, 10 nodos)
└── executaveis/
    ├── windows/sor_pnbi.exe       <- programa para Windows
    ├── linux/sor_pnbi             <- programa para Linux
    └── macos/sor_pnbi             <- programa para macOS
```

Este programa resolve qualquer problema de transporte de partículas **unidimensional, de fonte fixa** (não é só os dois exemplos do artigo): quantas regiões você quiser, com os materiais, espessuras e fonte que você quiser. Veja a seção 5.

---

## 2. Passo 1 — Baixar

1. Abra a página do repositório no GitHub no seu navegador.
2. Clique no botão verde **`<> Code`** (fica no alto da lista de arquivos, à direita).
3. No menu que abrir, clique em **`Download ZIP`**.
4. **Extraia** o `.zip` (normalmente baixado na pasta *Downloads*):
   * **Windows:** clique com o botão direito no arquivo → **Extrair tudo…** → **Extrair**.
   * **macOS:** dê dois cliques no arquivo.
   * **Linux:** botão direito → *Extrair aqui* (ou `unzip nome-do-arquivo.zip` no terminal).
5. Dentro da pasta extraída você deve ver `README.md`, `Problemas/` e `executaveis/`. **Essa é a pasta do projeto** — é nela que você vai trabalhar.

> Dica: evite deixar a pasta em um caminho com muitos acentos/espaços. Não é obrigatório, mas evita dores de cabeça.

---

## 3. Passo 2 — Executar

Escolha o executável da pasta `executaveis/` de acordo com o seu sistema: `windows/sor_pnbi.exe`, `linux/sor_pnbi` ou `macos/sor_pnbi`.

### Windows

* Dê **dois cliques** em `executaveis\windows\sor_pnbi.exe` — abre uma janela preta com o menu. A janela só fecha depois que você apertar Enter no fim.
* Ou, pelo terminal (PowerShell/cmd) dentro da pasta do projeto:
  ```
  .\executaveis\windows\sor_pnbi.exe
  ```

### Linux e macOS

Antes da primeira vez, é preciso autorizar o arquivo a ser executado (motivo de segurança do sistema, só precisa fazer uma vez):

```
chmod +x executaveis/linux/sor_pnbi      # no Linux
chmod +x executaveis/macos/sor_pnbi      # no macOS
```

> **macOS:** ao tentar abrir, o Gatekeeper provavelmente vai bloquear com a mensagem *"não é possível verificar o desenvolvedor"* (porque o programa não é assinado). Para liberar: clique com o **botão direito** (ou Control+clique) em `sor_pnbi` → **Abrir** → confirme **Abrir** na caixa de diálogo. É preciso fazer isso só uma vez; depois disso ele abre normalmente, inclusive pelo terminal.

Depois, pelo terminal, dentro da pasta do projeto:

```
./executaveis/linux/sor_pnbi        # Linux
./executaveis/macos/sor_pnbi        # macOS
```

### Modo interativo (o programa faz as perguntas)

```
=== SOR-PNBI: transporte de particulas em 1D (Full NBI, Partial NBI e SOR-PNBI) ===

Arquivo do problema [Enter = Problemas/Prob_Ex_2.txt]:          <- Enter, ou digite o caminho do seu .txt
Problema lido: N = 16 | regioes = 4 | nodos = 10

Metodos disponiveis:
  ...
Escolha o metodo [Enter = 4]:                                   <- digite um número de 0 a 8
Numero de execucoes (para media do tempo) [Enter = 5]:          <- Enter
```

1. **Arquivo do problema:** o caminho do `.txt` (veja a seção 5 para usar o seu próprio). Para o Exemplo 1 do artigo digite `Problemas/Prob_Ex_1.txt`. Você também pode arrastar o arquivo para dentro da janela do terminal.
2. **Método:** um dos números da tabela do início deste README. O **8** é ótimo para comparar os métodos 1 a 6 de uma vez. Se escolher o **7**, o programa pergunta qual valor de *w* usar (entre 1 e 2, tipicamente).
3. **Número de execuções:** o problema é resolvido esse número de vezes para obter um tempo médio confiável. Use `1` se só quiser o resultado.

Saída típica (método 4, Exemplo 2):

```
--- SOR-PNBI (w pela variancia) ---
Iteracoes: 143 | Tempo medio: 6.074e-04 s | w = 1.147578
Fluxo escalar: |  1.4988e+01  3.2418e+02  3.9701e+02  ...  3.2635e-07  1.6279e-13  5.4174e-19  |

Fluxo escalar completo salvo em: resultado_Prob_Ex_2_metodo4.txt
```

* **Iteracoes:** quantas varreduras foram necessárias até a convergência (desvio relativo do fluxo escalar ≤ 10⁻¹⁰).
* **Tempo medio:** tempo de uma execução completa em segundos.
* **w:** fator de relaxação usado (se houver um por região, aparecem todos: `w_R1 = … / w_R2 = …`).
* **Fluxo escalar:** os 3 primeiros e os 3 últimos valores (ou todos, se forem poucos), um por aresta de nodo.
* O **fluxo escalar completo** é gravado no arquivo `resultado_<problema>_metodo<n>.txt`, na pasta onde você executou o programa. Tem duas colunas — posição *x* (cm) e fluxo escalar — e pode ser aberto no Excel, Python, Origin, gnuplot etc. para fazer gráficos.

Mensagens que você pode ver:

* `ATENCAO: o metodo divergiu…` — o *w* escolhido fez o método divergir.
* `ATENCAO: limite de iteracoes atingido sem convergir.` — passou de 100 000 iterações.
* `AVISO: nao foi encontrado minimo em 1 < w < 2 (...). Usando w = 1.` — a estratégia de escolha de *w* não achou mínimo nesse problema, então o programa usa *w* = 1 (Partial NBI).

### Modo direto (sem perguntas, bom para automatizar)

```
sor_pnbi <arquivo do problema> <método> [execuções] [w]
```

Exemplos (Windows; em Linux/macOS troque `sor_pnbi.exe` por `./sor_pnbi` e o caminho pelo seu sistema):

```
executaveis\windows\sor_pnbi.exe Problemas\Prob_Ex_2.txt 4          (método 4, 5 execuções)
executaveis\windows\sor_pnbi.exe Problemas\Prob_Ex_1.txt 8 3        (métodos 1 a 6, 3 execuções cada)
executaveis\windows\sor_pnbi.exe Problemas\Prob_Ex_2.txt 7 1 1.3    (método 7, 1 execução, w = 1.3)
```

### Sobre o método 0 (w ótimo)

Para achar o *w* ótimo o programa monta a matriz de iteração **T** inteira, de tamanho `N × (nodos + 1)`, e calcula seus autovalores várias vezes. No Exemplo 2 (T de 176 × 176) leva menos de 1 segundo; no Exemplo 1 original (N = 32, T de 3232 × 3232) leva muito mais. Use-o em problemas pequenos, ou use o método 7 com um *w* que você já conheça.

---

## 4. Como saber se está tudo certo (valores de referência)

Rode com **5 execuções** e compare o número de **iterações** (o tempo depende do seu computador).

**Exemplo 2 — `Problemas/Prob_Ex_2.txt`** (N = 16, 4 regiões)

| Método | Iterações | *w* |
|:------:|:---------:|-----|
| 1 (Full NBI) | 204 | — |
| 2 (Partial NBI) | 197 | 1 |
| 3 | 197 | 1.000078 |
| 4 | 143 | 1.147578 |
| 5 | 197 | 1.000098 / 1.191504 / 1.015137 / 1.000098 |
| 6 | 138 | 1.164941 / 1.191504 / 1.015137 / 1.000000 |
| 0 (ótimo) | 40 | 1.4901 |
| 7 (com *w* = 1.3) | 99 | 1.3 |

**Exemplo 1 — `Problemas/Prob_Ex_1.txt`** (N = 32, 1 região, 100 nodos)

| Método | Iterações | *w* |
|:------:|:---------:|-----|
| 1 (Full NBI) | 2742 | — |
| 2 (Partial NBI) | 471 | 1 |
| 3 e 5 | 236 | 1.0411 |
| 4 e 6 | 97 | 1.0640 |

Os fluxos escalares finais do Exemplo 1 começam em `1.5691e+01  2.3122e+01  2.9003e+01 …`.

---

## 5. Criando o seu próprio problema (qualquer modelo 1D de fonte fixa)

O programa não está limitado aos exemplos do artigo: ele resolve **qualquer problema de transporte de partículas em geometria plana 1D, com fonte fixa**, com o número de regiões, materiais, espessuras e fontes que você quiser. Basta descrever o problema em um arquivo `.txt`.

Comece copiando `Problemas/Problema_Modelo.txt` — ele já é um problema válido de 1 região, pronto para rodar, com cada linha explicando o que colocar nela. Abra-o com qualquer editor de texto simples (Bloco de Notas, TextEdit, VS Code, etc.):

```
Ordem da quadratura N (numero PAR: 2, 4, 6, 8, 16, 32...):
8
Numero de Regioes do problema (nI):
1
Numero de Zonas Materiais distintas (nZM):
1
Mapeamento regiao -> material (nI valores; material comeca em 0):
0
Tamanho de cada regiao, em cm (nI valores):
10
Numero de nodos em que cada regiao e dividida (nI valores):
10
Condicoes de contorno: fluxo incidente na esquerda e na direita (xb xc):
0 0
Secao de choque total de cada material - sigma_t (nZM valores):
1
Secao de choque de espalhamento de cada material - sigma_s (nZM valores; deve ser menor que sigma_t):
0.5
Fonte externa fixa de cada material - Q (nZM valores):
1
```

O arquivo sempre tem **10 itens**, cada um ocupando **duas linhas**: uma linha de título (o texto não importa para o programa — pode editar ou deixar como está) e, logo abaixo, a linha com os valores, separados por espaço. Linhas totalmente em branco são ignoradas. Use **ponto** como separador decimal (`0.999`, não `0,999`).

Para um problema com mais de uma região (como o Exemplo 2, de 4 regiões), as linhas de **"Mapeamento"**, **"Tamanho"** e **"Numero de nodos"** precisam ter **um valor por região** (`nI` valores); as de **sigma_t**, **sigma_s** e **Fonte** precisam ter **um valor por material** (`nZM` valores, que podem se repetir entre regiões via o mapeamento). Veja `Prob_Ex_2.txt` como referência de um caso com 4 regiões e 4 materiais.

Regras físicas a observar: `sigma_s` (espalhamento) deve ser menor ou igual a `sigma_t` (total) em cada material; a quadratura `N` deve ser par. O programa confere a formatação e mostra uma mensagem de erro em português dizendo o que está errado, caso algo não bata.

Depois de editar, salve com outro nome (ex.: `Problemas/MeuProblema.txt`) e informe esse caminho ao programa — na pergunta "Arquivo do problema" (modo interativo) ou como primeiro argumento (modo direto, seção 3).

---

## 6. Problemas comuns

| Sintoma | O que fazer |
|---------|-------------|
| `Nao foi possivel abrir o arquivo…` | Confira o nome/caminho do `.txt` e execute o programa a partir da pasta do projeto (ou digite o caminho completo). |
| Linux/macOS: `Permission denied` ao executar | Rode `chmod +x` no arquivo (veja a seção 3) antes de tentar de novo. |
| macOS: *"não é possível verificar o desenvolvedor"* / *"malicious software"* | Clique com o botão direito no arquivo → **Abrir** → confirme **Abrir** (só precisa uma vez). Veja a seção 3. |
| Windows: antivírus/SmartScreen bloqueia o `.exe` | É comum com executáveis baixados da internet sem assinatura digital. Clique em **Mais informações** → **Executar assim mesmo** (no aviso do SmartScreen), ou libere no antivírus. |
| `Metodo invalido` | Use um número de 0 a 8. |
| Erro dizendo que faltam valores numa linha | Confira se a linha correspondente tem a quantidade certa de números (`nI` ou `nZM` valores — veja a seção 5). |

---

## 7. Licença e créditos

* A biblioteca **Eigen**, usada internamente pelo programa (compilada dentro dos executáveis), é distribuída sob a licença MPL 2.0 — veja `LICENSE_Eigen_MPL2.txt` e https://eigen.tuxfamily.org.
* Ao usar este programa em trabalhos acadêmicos, cite o artigo: *Técnica de sobre-relaxação sucessiva para o esquema iterativo de inversão nodal parcial em simulações de transporte S<sub>N</sub> de partículas neutras em geometria unidimensional cartesiana*.
* Artigo e programa: https://github.com/davidarriaux/Artigo-e-Codigo-ENMC-2026---Davi-Darriaux-Ferreira
* Autor: **Davi Darriaux Ferreira** — davi.darriaux@iprj.uerj.br (IPRJ/UERJ).
