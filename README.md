# Intensidades máximas anuais no Rio Grande do Sul

## Sobre

Este repositório compila um conjunto de códigos em R produzidos no contexto do mestrado no Programa de Pós-Graduação em Tecnologia Ambiental e Recursos Hídricos da Universidade de Brasília [(PPGTARH-UnB)](https://ptarh.unb.br/).

A pesquisa, intitulada "*Curvas de Intensidade-Duração-Frequência baseadas em invariância de escala e informação regional*", investiga alternativas para incluir informação adicional no ajuste de curvas IDF, utilizando o estado do Rio Grande do Sul como região de estudo.

Os códigos disponíveis neste repositório são referentes à estratégia que envolve a reparametrização da distribuição Generalizada de Valores Extremos incluindo uma estrutura de dependência temporal (*duration-dependent GEV*), seguindo a metodologia de [Ulrich et al. (2020)](https://doi.org/10.3390/w12113119). **Somente os códigos em R são disponibilizados**.

Foram utilizadas séries de precipitação subdiária compiladas por [Das Neves Almeida et al. (2025)](https://doi.org/10.1080/02626667.2025.2506193).
## Códigos

Os códigos estão localizados na pasta [R](R/). Os principais seguem a seguinte ordem de aplicação:
1. [extract_agg_imax](extract_agg_imax.R): extrair intensidades máximas anuais de precipitação para diferentes durações. Utiliza os scripts auxiliares [fun_filter_set](R/funcoes/fun_filter_set.R), [fun_group_ts](R/funcoes/fun_group_ts.R) e [fun_imax_agg_wateryear](IntensidadeMaximaRS/R/funcoes/fun_imax_agg_wateryear.R);
2. [fitprof_gev](fitprof_gev.R): ajusta a distribuição GEV por duração e estima intervalos de confiança com verossimilhança perfilada para parâmetros e quantis. Utiliza o script auxiliar [fgev](IntensidadeMaximaRS/R/funcoes/fgev.R);
3. [fitprof_dgev](IntensidadeMaximaRS/R/fitprof_dgev.R): ajusta o modelo d-GEV e estima intervalos de confiança com verossimilhança perfilada para parâmetros e quantis. Utiliza o script auxiliar [fdgev](IntensidadeMaximaRS/R/funcoes/fdgev.R);
4. [visualize_subdaily](visualize_subdaily.R): avalia os resultados dos ajustes e gera figuras e tabelas.

Demais códigos são auxiliares.