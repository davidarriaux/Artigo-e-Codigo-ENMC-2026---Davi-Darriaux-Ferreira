# Técnica de Sobre-Relaxação Sucessiva para o Esquema Iterativo de Inversão Nodal Parcial

Código e artigo referentes ao desenvolvimento de uma técnica de Sobre-Relaxação Sucessiva (SOR) para o esquema iterativo de Inversão Nodal Parcial (PNBI) em simulações de transporte (S_N) de partículas neutras em geometria unidimensional cartesiana.

O trabalho propõe o esquema iterativo SOR-PNBI, aplicado em conjunto com o método Matriz Resposta (RM), além de duas estratégias para a estimação do fator de relaxação ω. Os resultados numéricos mostram reduções expressivas no número de iterações necessárias para a convergência, quando comparado ao esquema PNBI convencional.

Trabalho apresentado no XXIX ENMC / XVII ECTM (Bento Gonçalves - RS, 2026).

## Conteúdo deste repositório

- **Artigo**: versão completa do trabalho submetido ao evento, com a fundamentação teórica, a metodologia e os resultados numéricos obtidos.
- **`Codigo ENMC 2026/`**: pasta com a implementação numérica do esquema SOR-PNBI e das estratégias de estimativa do fator de relaxação ω. Dentro dela há um README próprio, explicando como instalar as dependências e executar o código.

## Autores

- Davi D. Ferreira
- Leonardo R. C. Moraes
- Ricardo C. Barros

Universidade do Estado do Rio de Janeiro (UERJ), Instituto Politécnico – Nova Friburgo, RJ.
