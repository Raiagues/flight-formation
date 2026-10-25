# Reprodução do experimento

## Fonte selecionada

O relatório é “Filtro de Kalman Estendido Aplicado à Formação em Voo de Satélites”, de Emilly Raiane Rodrigues e Willer Gomes dos Santos, ENCITA 2023. Dentro do ZIP, a versão `original/NOVO/CODIGOS/NAO_LINEAR/tests/CODIGO_PROFESSOR/` contém os resultados salvos que correspondem à tabela. O restante do arquivo inclui versões anteriores e experimentos; não existe um identificador de versão da execução que originou a publicação.

## Parâmetros preservados

| Parâmetro | Valor |
| --- | --- |
| μ terrestre | 3,986004418 × 10¹⁴ m³/s² |
| Raio terrestre | 6.378.140 m |
| Semieixo maior do líder | 6.978.140 m |
| Excentricidade | 0,1 |
| θ inicial | 0 rad |
| (x, y, z) / a | (−0,01127; 0; 0,1) |
| (vx, vy, vz) / (n a) | (0,02; 0,02; 0) |
| Intervalo de medida | 5 s |
| σ de medida | 10√2 m por eixo |
| Dados arquivados | 4.601 amostras de t = 5 a 23.005 s |

O perigeu do líder resulta em altitude −97,814 km. Essa condição não representa uma órbita de voo livre realizável em torno de uma Terra sólida. Ela foi preservada para comparação matemática; o cenário circular na interface é explicitamente separado. A separação inicial é da ordem de 700 km, como resulta da normalização do código, e não foi reduzida para parecer uma formação próxima.

## O que foi encontrado na comparação

No `main.m` selecionado, `res = y - H*xb` é a **inovação**. `res2 = Xc(1:3) - xb(1:3)` é a **correção**. As estatísticas chamadas de RMSE e desvio padrão no arquivo de resultados são calculadas com essas séries a partir do índice MATLAB `500:end`. Elas não medem diretamente `medição − verdade` e `estimativa − verdade`.

Além disso, esse `main.m` atualiza o filtro com a primeira medição (marcada t = 5 s) antes de propagar o estado inicial de t = 0. Na implementação web, primeiro propagamos até o instante da medição. O original permanece sem alterações.

| Série / execução | X (m) | Y (m) | Z (m) |
| --- | ---: | ---: | ---: |
| Relatório: coluna de medição | 15,2502 | 16,2529 | 14,9726 |
| Relatório: coluna de estimativa | 3,1836 | 4,7768 | 2,9203 |
| Reexecução histórica: RMS da inovação | 14,9210 | 14,8623 | 14,6561 |
| Reexecução histórica: RMS da correção | 2,7319 | 2,7159 | 2,6771 |
| Execução alinhada: RMSE medição − verdade | 14,0820 | 13,9977 | 13,6991 |
| Execução alinhada: RMSE estimativa − verdade | 4,9606 | 4,7766 | 4,6374 |

**A reprodução numérica exata da tabela não foi obtida.** A reexecução histórica preserva a ordem atualização/propagação, a covariância inicial e Q/R da versão arquivada, mas usa RK4 e forma de Joseph. Os valores não devem ser apresentados como uma confirmação exata dos resultados publicados. As versões de desenvolvimento, a ausência de semente no `randn` original, a ordem temporal, o arredondamento dos `.dat` e diferenças de integrador são fatores a investigar. Não é possível atribuir a discrepância a uma única causa com o material disponível.

O relatório descreve reduções de 21%, 29% e 20%. Pelos próprios números da tabela, esses percentuais correspondem aproximadamente à **fração remanescente** (estimado / medido); as reduções correspondentes seriam aproximadamente 79%, 71% e 80%. A interface não repete esses percentuais como uma nova validação.

## Implementação web

- Dinâmica relativa não linear de 10 estados; líder e seguidor obedecem à gravitação central de dois corpos.
- RK4 com passo máximo de 5 s para estado e covariância.
- Jacobiano analítico; covariância contínua `Pdot = F P + P Fᵀ + Q`.
- Atualizações escalares sequenciais equivalentes à observação vetorial com R diagonal, em forma de Joseph.
- Q diagonal reproduzido do `edo.m` selecionado. A verdade simulada é determinística; Q representa incerteza do modelo no filtro, não uma força aleatória aplicada aos satélites.
- Gerador de medições com semente fixa 2023 (seguidor 01) e 2024 (seguidor 02), permitindo comparar o ajuste de ruído. O ruído não move a trajetória verdadeira nem é ampliado na cena.
- Ruído zero usa piso numérico de 10⁻¹² m² em R para evitar inversão singular.
- Gráfico mostra erros reais por eixo na janela móvel de até 240 amostras. O RMSE da interface usa essa janela e não o descarte inicial usado na tabela acima.
- A simulação termina em 23.005 s e pode ser reiniciada; não há salto automático para uma órbita anterior.
- Segundo seguidor: Z inicial espelhado, propagação e filtro independentes. Não há forças ou medições entre seguidores.

## Referenciais e relógios

A dinâmica original é definida no plano orbital. A transformação para ECI usa uma inclinação visual fixa de 35° e longitude do nó de 0°. No Three.js, ECI `(X,Y,Z)` é desenhado como `(X,Z,−Y)`, mantendo o polo norte no eixo vertical. A coordenada LVLH é X radial, Y tangencial e Z normal.

A textura terrestre gira pelo ângulo sideral de Greenwich calculado para UTC (aproximação de UT1). O Sol usa as [coordenadas aproximadas do USNO](https://aa.usno.navy.mil/faq/sun_approx). O horário sideral segue a expressão GMST referenciada na [documentação NOAA de referenciais](https://www.ngdc.noaa.gov/stp/space-weather/online-publications/stp_sii/spades/algorithm-theoretic-basis-documents/atbds_mag/atbd_mag07_coordinate-systems_v1-3.pdf). São aproximações visuais, sem correções completas de precessão/nutação ou consulta a efemérides de alta precisão.

O horário da Terra e do Sol é o relógio real do dispositivo. O tempo orbital inicia em zero e pode ser acelerado/pausado independentemente. Nenhum dos satélites representa rastreamento de uma missão real ou usa TLE. A rotação da Terra é realista para a hora atual; a posição orbital dos satélites é a do experimento.

## Verificações automatizadas

```bash
node scripts/verify-original.mjs
npm test
npm run reproduce
npm run build
```

Os testes verificam a dinâmica relativa contra uma integração cartesiana inercial independente por 23.000 s (erro abaixo de 0,3 m), conservação de energia e momento angular do líder, Jacobiano por diferenças finitas, trajetória contra os dados MATLAB arquivados (tolerância de 0,02 m; maior diferença observada de 0,00517 m), redução do erro pelo FKE, covariância simétrica, ruído zero, independência dos seguidores e verificações de astronomia em J2000 e solstícios.

Os resultados detalhados da reprodução ficam em [`reproduction.json`](reproduction.json). Nenhum script de reprodução modifica os arquivos em `original/`. O workflow inclui teste de navegação em Chromium, carregamento do modelo Draco, controles e layout móvel, com screenshots como artefato. A execução local do navegador/servidor não foi permitida pelo ambiente de desenvolvimento; o resultado desse teste deve ser consultado na execução do GitHub Actions.
