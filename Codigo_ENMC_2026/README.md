# SOR-PNBI — Transporte de partículas em 1D (Full NBI, Partial NBI e SOR-PNBI)

Código em C++ que acompanha o artigo **"TÉCNICA DE SOBRE-RELAXAÇÃO SUCESSIVA PARA O ESQUEMA ITERATIVO DE INVERSÃO NODAL PARCIAL EM SIMULAÇÕES DE TRANSPORTE S<sub>N</sub> DE PARTÍCULAS NEUTRAS EM GEOMETRIA UNIDIMENSIONAL CARTESIANA"**.

Repositório (artigo e código): https://github.com/davidarriaux/Artigo-e-Codigo-ENMC-2026---Davi-Darriaux-Ferreira

Ele resolve a equação de transporte de partículas (ordenadas discretas, geometria plana 1D, meio com várias regiões) pelos métodos iterativos:

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

Você **não precisa instalar nenhuma biblioteca**: a biblioteca de álgebra linear usada (Eigen) já vem dentro da pasta `src/Eigen`. Só é preciso ter um **compilador C++**.

> Testado com g++ 13.1 (MinGW-w64) no Windows 11, compilando pelo terminal e pelo CMake. O código usa apenas C++17 padrão e deve compilar igualmente com g++/clang no Linux e no macOS.

---

## Sumário

1. [O que tem neste repositório](#1-o-que-tem-neste-repositório)
2. [Passo 1 — Baixar o código](#2-passo-1--baixar-o-código)
3. [Passo 2 — Instalar um compilador C++](#3-passo-2--instalar-um-compilador-c)
4. [Passo 3 — Compilar](#4-passo-3--compilar)
5. [Passo 4 — Executar](#5-passo-4--executar)
6. [Como saber se está tudo certo (valores de referência)](#6-como-saber-se-está-tudo-certo-valores-de-referência)
7. [Usando o seu próprio problema (formato do arquivo)](#7-usando-o-seu-próprio-problema-formato-do-arquivo)
8. [Problemas comuns](#8-problemas-comuns)
9. [Organização do código](#9-organização-do-código)
10. [Licença e créditos](#10-licença-e-créditos)

---

## 1. O que tem neste repositório

```
.
├── README.md              <- este arquivo
├── CMakeLists.txt         <- receita de compilação (para quem usa CMake / CLion)
├── LICENSE_Eigen_MPL2.txt <- licença da biblioteca Eigen
├── Problemas/
│   ├── Prob_Ex_1.txt      <- Exemplo 1 do artigo (1 região, N = 32, 100 nodos)
│   └── Prob_Ex_2.txt      <- Exemplo 2 do artigo (4 regiões, N = 16, 10 nodos)
└── src/                   <- código-fonte
    ├── main.cpp           <- menu e impressão de resultados
    ├── Problema.h/.cpp    <- leitura do arquivo do problema e quadratura
    ├── NBI.h/.cpp         <- matrizes resposta, varreduras e convergência
    ├── Relaxacao.h/.cpp   <- escolha do fator de relaxação w
    ├── Metodos.h/.cpp     <- executa um método e mede o tempo
    └── Eigen/             <- biblioteca Eigen (não precisa mexer)
```

---

## 2. Passo 1 — Baixar o código

1. Abra a página do repositório no GitHub no seu navegador.
2. Clique no botão verde **`<> Code`** (fica no alto da lista de arquivos, à direita).
3. No menu que abrir, clique em **`Download ZIP`**.
4. O navegador vai baixar um arquivo `.zip` (normalmente para a pasta *Downloads*).
5. **Extraia** o `.zip`:
   * **Windows:** clique com o botão direito no arquivo → **Extrair tudo…** → **Extrair**.
   * **macOS:** dê dois cliques no arquivo.
   * **Linux:** botão direito → *Extrair aqui* (ou `unzip nome-do-arquivo.zip` no terminal).
6. Entre na pasta que foi criada e, dentro dela, na subpasta **`Codigo_ENMC_2026`** (o ZIP também traz o artigo). Nessa subpasta você deve ver `README.md`, `CMakeLists.txt`, a pasta `src` e a pasta `Problemas`. **Essa é a "pasta do projeto"** — vamos usá-la nos próximos passos.

> Dica: evite deixar a pasta em um caminho com muitos acentos/espaços (ex.: `C:\Users\Maria\Downloads\SOR` é melhor que um caminho muito longo). Não é obrigatório, mas evita dores de cabeça.

---

## 3. Passo 2 — Instalar um compilador C++

Um *compilador* é o programa que transforma o código-fonte (`.cpp`) em um programa executável. **Se você já tem um compilador C++ (g++, clang, CLion, Code::Blocks…), pule para o Passo 3.**

Para saber se você já tem o g++, abra um terminal (veja o Passo 3 para saber como) e digite:

```
g++ --version
```

Se aparecer algo como `g++ (...) 13.1.0`, você já tem. Se aparecer *"não é reconhecido como um comando"* ou *"command not found"*, instale conforme seu sistema abaixo. É preciso um compilador que aceite **C++17** (qualquer g++ a partir da versão 7 serve).

### Windows — g++ pelo MSYS2

1. Acesse **https://www.msys2.org** e baixe o instalador (`msys2-x86_64-….exe`).
2. Execute o instalador e clique em *Avançar* até o fim, mantendo a pasta padrão `C:\msys64`. Ao terminar, deixe marcada a opção *Run MSYS2 now* (ou abra depois pelo Menu Iniciar o **"MSYS2 UCRT64"**).
3. Na janela preta que abrir, cole o comando abaixo e pressione Enter. Quando perguntar `Proceed with installation? [Y/n]`, digite `Y` e Enter:
   ```
   pacman -S --needed mingw-w64-ucrt-x86_64-gcc
   ```
4. Agora é preciso avisar o Windows onde o g++ ficou (colocá-lo no **PATH**):
   1. Aperte a tecla **Windows**, digite **variáveis de ambiente** e abra *"Editar as variáveis de ambiente do sistema"*.
   2. Clique em **Variáveis de Ambiente…**.
   3. Na parte de cima (*Variáveis de usuário*), selecione **Path** e clique em **Editar…**.
   4. Clique em **Novo** e digite: `C:\msys64\ucrt64\bin`
   5. Clique **OK** em todas as janelas.
5. **Feche e abra de novo** qualquer terminal que estivesse aberto e teste:
   ```
   g++ --version
   ```
   Deve aparecer a versão do g++.

> Alternativa mais simples para quem nunca programou: instalar o **CLion** (https://www.jetbrains.com/clion/ — gratuito para uso não comercial/estudantes). Ele já traz um compilador (MinGW) embutido. Veja a opção *B* do Passo 3.

### Linux (Ubuntu/Debian)

```
sudo apt update
sudo apt install g++
```

(Em Fedora: `sudo dnf install gcc-c++`.)

### macOS

Abra o **Terminal** e rode:

```
xcode-select --install
```

Aceite a instalação na janela que aparecer. Isso instala o compilador `clang++` e também o comando `g++` (que no macOS é o clang).

---

## 4. Passo 3 — Compilar

Escolha **uma** das opções.

### Opção A — Pelo terminal (uma linha só; funciona em qualquer sistema)

**1) Abra um terminal dentro da pasta do projeto.**

* **Windows:** abra a pasta do projeto no Explorador de Arquivos, clique na **barra de endereço** (onde aparece o caminho), apague o texto, digite `cmd` e pressione **Enter**. Abre-se uma janela preta já dentro da pasta certa.
* **macOS:** abra o Terminal, digite `cd ` (com um espaço no final), arraste a pasta do projeto para dentro da janela do Terminal e pressione Enter.
* **Linux:** abra o terminal, `cd` até a pasta do projeto (ou, no gerenciador de arquivos, *botão direito → Abrir no terminal*).

Confirme que está no lugar certo listando os arquivos: `dir` (Windows) ou `ls` (macOS/Linux). Deve aparecer `README.md`, `src`, `Problemas`…

**2) Rode o comando de compilação.**

Windows:

```
g++ -O3 -std=c++17 -static src/main.cpp src/Problema.cpp src/NBI.cpp src/Relaxacao.cpp src/Metodos.cpp -o sor_pnbi
```

Linux e macOS (sem o `-static`):

```
g++ -O3 -std=c++17 src/main.cpp src/Problema.cpp src/NBI.cpp src/Relaxacao.cpp src/Metodos.cpp -o sor_pnbi
```

O que significa cada parte: `-O3` liga as otimizações (deixa o programa muito mais rápido — não esqueça!), `-std=c++17` escolhe a versão da linguagem, os cinco arquivos `.cpp` são o código, e `-o sor_pnbi` é o nome do programa que será criado.

**3) Espere.** A compilação leva de **30 segundos a 1–2 minutos** (a biblioteca Eigen é grande). Enquanto o cursor não voltar, está compilando. Se não aparecer nenhuma mensagem de erro, terminou: surgiu na pasta o arquivo `sor_pnbi.exe` (Windows) ou `sor_pnbi` (Linux/macOS).

### Opção B — CLion

1. Abra o CLion → **Open** → selecione a **pasta do projeto** (a que contém o `CMakeLists.txt`) → **OK** → *Trust Project*.
2. Aguarde o CLion terminar de carregar o CMake (barra de progresso no canto inferior direito).
3. No canto superior direito, escolha o alvo **`sor_pnbi`** e clique no botão verde ▶ **Run** (ou aperte `Shift+F10`).
4. O programa abre na aba **Run** embaixo; clique dentro dela para digitar as respostas do menu.

(O `CMakeLists.txt` já força a otimização máxima, inclusive no perfil *Debug*.)

### Opção C — Pelo CMake, no terminal

Com o CMake instalado (https://cmake.org), dentro da pasta do projeto:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

O programa será criado dentro da pasta `build` (ou `build/Release` no Visual Studio).

### Outras IDEs (Code::Blocks, Dev-C++, Visual Studio…)

Crie um projeto de console C++ vazio, **adicione os 5 arquivos `.cpp` da pasta `src`** (`main`, `Problema`, `NBI`, `Relaxacao`, `Metodos`), selecione o padrão **C++17** e a configuração **Release** (otimizada). Não é preciso configurar nenhum caminho de include: o `Eigen` é encontrado sozinho porque está ao lado dos arquivos. Ao executar, a *pasta de trabalho* deve ser a pasta do projeto (ou informe o caminho completo do arquivo do problema quando o programa perguntar).

---

## 5. Passo 4 — Executar

### Modo interativo (o programa faz as perguntas)

No mesmo terminal, dentro da pasta do projeto:

* **Windows:** `sor_pnbi.exe`   (no PowerShell: `.\sor_pnbi.exe`)
* **Linux / macOS:** `./sor_pnbi`

(No Windows também dá para dar **dois cliques** em `sor_pnbi.exe`; a janela só fecha depois que você apertar Enter no fim.)

O programa pergunta, e para cada pergunta **basta apertar Enter para aceitar o valor entre colchetes**:

```
=== SOR-PNBI: transporte de particulas em 1D (Full NBI, Partial NBI e SOR-PNBI) ===

Arquivo do problema [Enter = Problemas/Prob_Ex_2.txt]:          <- Enter
Problema lido: N = 16 | regioes = 4 | nodos = 10

Metodos disponiveis:
  1 - Full NBI
  ...
Escolha o metodo [Enter = 4]:                                   <- digite um número de 0 a 8
Numero de execucoes (para media do tempo) [Enter = 5]:          <- Enter
```

1. **Arquivo do problema:** o caminho do `.txt`. Para o Exemplo 1 do artigo digite `Problemas/Prob_Ex_1.txt`. Você também pode arrastar o arquivo para dentro da janela do terminal.
2. **Método:** um dos números da tabela do início deste README. O **8** é ótimo para comparar os métodos 1 a 6 de uma vez. Se escolher o **7**, o programa pergunta qual valor de *w* usar (entre 1 e 2, tipicamente).
3. **Número de execuções:** o problema é resolvido esse número de vezes para obter um tempo médio confiável (descartando o mais rápido e o mais lento quando são mais de 2). Use `1` se só quiser o resultado.

Saída típica (método 4, Exemplo 2):

```
--- SOR-PNBI (w pela variancia) ---
Iteracoes: 143 | Tempo medio: 6.074e-04 s | w = 1.147578
Fluxo escalar: |  1.4988e+01  3.2418e+02  3.9701e+02  ...  3.2635e-07  1.6279e-13  5.4174e-19  |

Fluxo escalar completo salvo em: resultado_Prob_Ex_2_metodo4.txt
```

* **Iteracoes:** quantas varreduras foram necessárias até a convergência (desvio relativo do fluxo escalar ≤ 10⁻¹⁰).
* **Tempo medio:** tempo de uma execução completa em segundos, **incluindo o cálculo do *w*** nos métodos que o calculam.
* **w:** fator de relaxação usado (se houver um por região, aparecem todos: `w_R1 = … / w_R2 = …`).
* **Fluxo escalar:** os 3 primeiros e os 3 últimos valores (ou todos, se forem poucos), um por aresta de nodo.
* O **fluxo escalar completo** é gravado no arquivo `resultado_<problema>_metodo<n>.txt`, na pasta onde você executou o programa. Ele tem duas colunas — posição *x* (cm) e fluxo escalar — e pode ser aberto no Excel, Python, Origin, gnuplot etc. para fazer gráficos.

Mensagens que você pode ver:

* `ATENCAO: o metodo divergiu…` — o *w* escolhido fez o método divergir.
* `ATENCAO: limite de iteracoes atingido sem convergir.` — passou de 100 000 iterações.
* `AVISO: nao foi encontrado minimo em 1 < w < 2 (...). Usando w = 1.` — a estratégia de escolha de *w* não achou mínimo nesse problema, então o programa usa *w* = 1 (Partial NBI).
### Modo direto (sem perguntas, bom para automatizar)

```
sor_pnbi <arquivo do problema> <método> [execuções] [w]
```

Exemplos:

```
sor_pnbi Problemas/Prob_Ex_2.txt 4          (método 4, 5 execuções)
sor_pnbi Problemas/Prob_Ex_1.txt 8 3        (métodos 1 a 6, 3 execuções cada)
sor_pnbi Problemas/Prob_Ex_2.txt 7 1 1.3    (método 7, 1 execução, w = 1.3)
```

(No Linux/macOS use `./sor_pnbi`; no PowerShell, `.\sor_pnbi.exe`.)

### Sobre o método 0 (w ótimo)

Para achar o *w* ótimo o programa monta a matriz de iteração **T** inteira, de tamanho `N × (nodos + 1)`, e calcula seus autovalores várias vezes. No Exemplo 2 (T de 176 × 176) leva menos de 1 segundo; com N = 8 e 100 nodos (T de 808 × 808) levou cerca de 50 s; no Exemplo 1 original (N = 32, T de 3232 × 3232) leva **muito** mais. Use-o em problemas pequenos, ou use o método 7 com um *w* que você já conheça.

---

## 6. Como saber se está tudo certo (valores de referência)

Rode com **5 execuções** e compare o número de **iterações** (o tempo depende do seu computador). Estes valores foram obtidos com o código deste repositório:

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

## 7. Usando o seu próprio problema (formato do arquivo)

O arquivo do problema é um `.txt` com **10 itens**. Cada item ocupa **duas linhas**: uma linha de título (o programa **ignora** o texto do título, mas a linha precisa existir) e, logo abaixo, a linha com os valores. Linhas totalmente em branco são ignoradas. Use **ponto** como separador decimal (`0.999`, não `0,999`).

Exemplo comentado (`Prob_Ex_2.txt`). Os textos após `←` **não** fazem parte do arquivo:

```
Ordem da quadratura:                     ← título
16                                       ← N: ordem da quadratura (número PAR: 2, 4, 6, 8, 16, 32…)
Número de Regiões:
4                                        ← nI: quantidade de regiões do meio
Número de Zonas Materiais:
4                                        ← nZM: quantidade de materiais diferentes
Mapeamento regiões:
0 1 2 3                                  ← material de cada região (começa em 0): região 1 usa o material 0, etc.
Tamanho das Regiões:
40 30 20 10                              ← espessura de cada região, em cm
Número de Nodos:
4 3 2 1                                  ← em quantos nodos cada região é dividida
Condições de contorno (xb,xc):
0 0                                      ← fluxo incidente no contorno esquerdo e no direito
Seção de choque total:
2 2 1.5 1.2                              ← σ_t de cada material
Seção de choque de espalhamento:
1.998 1.8 0.75 0.36                      ← σ_s de cada material
Fonte Externa:
1 0.8 0 0                                ← fonte Q de cada material
```

Regras: as linhas de "Mapeamento", "Tamanho" e "Número de Nodos" precisam ter **`nI` valores**; as de σ_t, σ_s e Fonte precisam ter **`nZM` valores**. O programa confere tudo e mostra uma mensagem de erro em português dizendo o que está errado.

Para criar um problema novo, copie `Prob_Ex_2.txt`, edite os valores e rode o programa apontando para o novo arquivo.

---

## 8. Problemas comuns

| Sintoma | O que fazer |
|---------|-------------|
| `g++ não é reconhecido…` / `command not found` | O compilador não está instalado ou não está no PATH. Refaça o Passo 2 (no Windows, **feche e reabra o terminal** depois de mexer no PATH). |
| `fatal error: Eigen/Dense: No such file or directory` | Você compilou de fora da pasta do projeto ou moveu os arquivos. Rode o comando **na pasta que contém `src`** e mantenha a pasta `src/Eigen` junto dos `.cpp`. |
| A compilação parece travada | É normal: leva de 30 s a 2 min por causa do Eigen. Espere o cursor voltar. |
| `Nao foi possivel abrir o arquivo…` | Confira o nome/caminho do `.txt` e execute o programa a partir da pasta do projeto (ou digite o caminho completo). |
| Ao dar dois cliques no `.exe`, aparece erro de DLL (`libstdc++-6.dll`…) | Compile no Windows com a opção `-static`, como no comando do Passo 3. |
| O programa está muito lento | Confirme que compilou com `-O3` (Opção A) ou com o `CMakeLists.txt` deste repositório. O método 0 é lento por natureza (veja a seção do método 0). |
| `Metodo invalido` | Use um número de 0 a 8. |

---

## 9. Organização do código

Para quem quiser ler ou modificar:

* **`Problema.cpp`** — lê e valida o `.txt`; calcula a quadratura de Gauss-Legendre (método de Golub-Welsch), a lista de nodos e as espessuras.
* **`NBI.cpp`** — para cada região calcula autovalores/autovetores, o fluxo particular e as matrizes resposta **R** e **S**; contém as varreduras do Full NBI e do Partial NBI/SOR-PNBI e o critério de convergência.
* **`Relaxacao.cpp`** — as estratégias de escolha do *w* (métodos 0 e 3–6). A diagonal de **T** para N = 2 é calculada por recorrência em O(n).
* **`Metodos.cpp`** — executa o método escolhido (várias vezes, se pedido) e mede o tempo.
* **`main.cpp`** — menu, leitura do teclado e impressão/gravação dos resultados.

Parâmetros de convergência (tolerância `1e-10` e máximo de `100000` iterações) ficam no início de `Problema.h`.

---

## 10. Licença e créditos

* A biblioteca **Eigen** (pasta `src/Eigen`) é distribuída sob a licença MPL 2.0 — veja `LICENSE_Eigen_MPL2.txt` e https://eigen.tuxfamily.org.
* Ao usar este código em trabalhos acadêmicos, cite o artigo: *Técnica de sobre-relaxação sucessiva para o esquema iterativo de inversão nodal parcial em simulações de transporte S<sub>N</sub> de partículas neutras em geometria unidimensional cartesiana*.
* Artigo e código: https://github.com/davidarriaux/Artigo-e-Codigo-ENMC-2026---Davi-Darriaux-Ferreira
* Autor: **Davi Darriaux Ferreira** — davi.darriaux@iprj.uerj.br (IPRJ/UERJ).
