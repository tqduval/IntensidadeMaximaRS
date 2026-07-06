---
title: Funcionamento das funções
---
## Funções do modelo d-GEV
As funções com sufixo `dgev` correspondem às funções em que a distribuição GEV é construída baseada na incorporação da invariância de escala e do modelo de Koutsoyiannis et al. (1998).

Três RScripts no total compõem as aplicações do modelo, alterando somente a forma como os intervalos de confiança são estimados. O modelo principal é o que integra o método mais robusto do perfil de verossimilhança.

As funções dos RScripts [fitdelta_dgev.R](R/funcoes/fitdelta_dgev.R) e [fitboot_dgev.R](R/funcoes/fitboot_dgev.R) foram integradas em um único script [fit_dgev.R](R/funcoes/fit_dgev.R) para compor uma única aplicação que permita a utilização de ambas as formas de estimar intervalos de confiança. A aplicação de [fitprof_dgev.R](R/funcoes/fitprof_dgev.R) ocorre em um script separado, pois ela será utilizada para gerar os resultados finais.
### d-GEV com perfil de verossimilhança
As seguintes funções compõem o Rscript `fitprof_dgev.R`:
1. `fdgev()`: realiza o ajuste da distribuição d-GEV utilizando máxima verossimilhança (convencional ou generalizada). Carrega internamente as funções `nll.dgev()` --- log-verossimilhança negativa para minimização pela função `optim()`  --- e `grad.ll.dgev()` --- funções gradientes da log-verossimilhança;
2. `prof.nll.dgev()`: responsável por selecionar um dos parâmetros da distribuição e passá-lo 