# Flight Formation

Simulação de formação em voo de satélites, com um líder em dinâmica de dois corpos e seguidores descritos pela dinâmica relativa não linear no referencial LVLH. O Filtro de Kalman Estendido (FKE) estima a posição relativa a partir de medições com ruído.

**[Abrir simulação](https://raiagues.github.io/flight-formation/)** · **[Relatório original](research/relatorio.pdf)** · **[Reprodução e validação](research/REPRODUCTION.md)**

Este repositório reúne o projeto de iniciação científica de **Emilly Raiane Rodrigues**, sob orientação de **Willer Gomes dos Santos**, apresentado no **XXVIII ENCITA / ITA, 2023**, e uma implementação web para explorar a dinâmica estudada.

## Visualização

- Terra com mapas NASA de superfície e luzes noturnas; rotação sideral e iluminação solar calculadas a partir do horário UTC atual.
- Modelo 3D do CubeSat MiRaTa, disponibilizado pela NASA, com estrutura, instrumentos e painéis solares.
- Vistas da órbita geocêntrica e da formação no referencial LVLH do líder.
- Pausa, reinício e velocidade orbital de 1× a 600×.
- Eixos e trajetórias que podem ser ativados ou desativados.
- Segundo seguidor funcional, **desligado inicialmente**, com condição inicial espelhada em Z e filtro próprio.
- Ruído de posição ajustável entre 0 e 50 m por eixo; gráfico do erro medido e estimado em relação à verdade simulada.

Arraste para girar a câmera; use a roda do mouse ou o gesto de pinça para aproximar. Com a cena focada, **Espaço** pausa e **R** restaura a câmera. Alterar o ruído ou o cenário reinicia o experimento com uma semente fixa.

O relógio da Terra e do Sol permanece em tempo real, inclusive durante a pausa. A velocidade selecionada afeta apenas a dinâmica orbital dos satélites. As imagens são mapas estáticos: não representam meteorologia ao vivo.

## Código original

Os **101 arquivos** do `codigo.zip` foram preservados, byte a byte, em [`original/NOVO`](original/NOVO), incluindo scripts MATLAB, dados, figuras, versões experimentais e bibliografia. O relatório fornecido está em [`research/relatorio.pdf`](research/relatorio.pdf). Os hashes de cada arquivo estão em [`research/original-sha256.json`](research/original-sha256.json).

A pasta [`CODIGO_PROFESSOR`](original/NOVO/CODIGOS/NAO_LINEAR/tests/CODIGO_PROFESSOR) contém um `RESULTADOS_FINAIS.txt` com os valores associados à tabela do relatório. As outras pastas são versões de desenvolvimento e **não devem ser adicionadas simultaneamente ao path do MATLAB**, pois repetem nomes como `main.m`, `medidas.m` e `edo.m`.

Para executar essa versão no MATLAB, a partir da raiz do repositório:

```matlab
cd('original/NOVO/CODIGOS/NAO_LINEAR/tests/CODIGO_PROFESSOR')
main
```

O script usa os arquivos `.dat` existentes e pode acrescentar conteúdo a `RESULTADOS_FINAIS.txt`. Execute em uma cópia caso queira manter o arquivo original intacto. `medidas` gera novas medições aleatórias e sobrescreve os `.dat`; não é necessário executá-lo para usar os dados arquivados. MATLAB não foi executado nesta implementação; a reprodução automatizada abaixo usa JavaScript e os dados originais.

## Executar a página localmente

Requisitos: **Node.js 20+**, **Python 3** e navegador com **WebGL 2**. As bibliotecas, os modelos e as texturas estão incluídos; não há dependência de CDN nem chave de API.

```bash
npm test
node scripts/verify-original.mjs
npm run reproduce
npm run build
python3 -m http.server 4173 --directory dist --bind 127.0.0.1
```

Abra <http://127.0.0.1:4173>. Para editar a interface, `npm run dev` serve diretamente `web/`; o link do relatório fica disponível na versão compilada em `dist/`.

O build apenas copia os arquivos estáticos e o relatório. A cada push em `main`, o workflow executa os testes, a reprodução, o build e os testes em Chromium antes de publicar no GitHub Pages. Os screenshots de desktop, formação e celular são salvos no artefato `browser-screenshots` da execução.

Para repetir os testes de navegador localmente, com o servidor ativo:

```bash
npm install --no-save --package-lock=false playwright@1.56.1
npx playwright install chromium
node scripts/browser-smoke.mjs
```

## Modelo e resultados

O estado é `[x, y, z, vx, vy, vz, r₀, ṙ₀, θ₀, θ̇₀]`, em unidades SI. A versão web integra as equações não lineares do relatório por RK4, com passo máximo de 5 s, e usa o Jacobiano analítico no FKE. As medidas são assimiladas no instante correto; a atualização da covariância usa a forma de Joseph.

Há duas limitações relevantes para interpretar a demonstração:

1. **A órbita original é um cenário matemático.** Com `a = 6.978.140 m` e `e = 0,1`, o perigeu é `r = 6.280.326 m`, isto é, **97,814 km abaixo** do raio terrestre adotado. Não foram incluídos colisão, atmosfera ou arrasto. A página mantém esses parâmetros e avisa quando o líder está sob a superfície. A opção circular a 600 km é uma adaptação separada.
2. **A tabela do relatório ainda não foi reproduzida exatamente.** O código original calcula estatísticas da inovação e da correção do filtro, que não são o mesmo que erro em relação à verdade. A implementação atual, usando os `.dat` arquivados e alinhamento temporal correto, encontra RMSE de posição de aproximadamente **14,08 / 14,00 / 13,70 m** nas medidas e **4,96 / 4,78 / 4,64 m** nas estimativas (X / Y / Z), após o descarte inicial original. Veja a [análise completa](research/REPRODUCTION.md).

A inclinação de 35° é uma escolha de visualização, pois o estudo fornece apenas a dinâmica relativa e o movimento do líder no plano orbital. Os CubeSats têm tamanho ampliado para serem visíveis. Os eixos são LVLH, sem modelo de atitude. O modelo MiRaTa representa visualmente os veículos; não identifica um satélite efetivamente simulado na pesquisa. Estrelas decorativas e tamanho/distância visuais do Sol não estão em escala astronômica. A direção solar e o terminador usam o horário real.

## Estrutura

```text
original/NOVO/          Projeto recebido, sem alterações
research/              Relatório, hashes e análise de reprodução
web/                   Interface, simulação e recursos 3D
scripts/               Build, reprodução e verificação
tests/                 Validação numérica
.github/workflows/     Testes e publicação no Pages
```

## Créditos

Pesquisa: Emilly Raiane Rodrigues (UFSM, bolsista PIBIC-CNPq) e Willer Gomes dos Santos (ITA). Apoio: ITA, grupo E2 MoC e CNPq, conforme o relatório.

Modelo 3D e imagens: NASA e colaboradores. Renderização: Three.js; compressão geométrica: Draco. Fontes, versões, licenças e alterações estão em [`ASSETS.md`](ASSETS.md). Os materiais originais e bibliográficos mantêm os direitos de seus respectivos autores; nenhuma licença geral é atribuída a esse conjunto.
