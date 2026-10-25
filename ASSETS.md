# Recursos, fontes e direitos

Todos os recursos necessários à simulação são servidos pelo próprio site, sem chamadas externas durante a execução.

| Arquivo | Origem / crédito | Uso |
| --- | --- | --- |
| `web/assets/cubesat-mirata.glb` | [NASA — CubeSat MiRaTa](https://science.nasa.gov/3d-resources/cubesat-mirata/), NASA / Christopher R. Meaney | GLB original, sem alteração; centralizado e redimensionado em tempo de execução. |
| `web/assets/earth-day.jpg` | [NASA Blue Marble](https://science.nasa.gov/resource/blue-marble/), NASA Goddard / Reto Stöckli e Robert Simmon | Mapa equiretangular original 2048×1024. |
| `web/assets/earth-night.jpg` | [NASA Earth at Night](https://earthobservatory.nasa.gov/features/NightLights), NASA / NOAA DMSP | Mapa equiretangular original 2400×1200; composição estática. |
| `web/vendor/three*.js`, loaders, controls, utilities, environment | [Three.js r180](https://github.com/mrdoob/three/tree/r180), MIT | Cópia local de módulos do pacote 0.180.0; imports reescritos para caminhos relativos. Licença em `THREE-LICENSE.txt`. |
| `web/vendor/draco/` | [Google Draco](https://github.com/google/draco), Apache 2.0 | Decodificador WASM distribuído com Three.js r180. Licença em `DRACO-LICENSE.txt`. |

Downloads diretos usados:

- [Modelo MiRaTa GLB](https://assets.science.nasa.gov/content/dam/science/cds/3d/resources/model/cubesat---mirata/CubeSat%20-%20MiRaTa.glb)
- [Blue Marble, land / ocean / ice](https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57730/land_ocean_ice_2048.jpg)
- [Earth lights](https://eoimages.gsfc.nasa.gov/images/imagerecords/55000/55167/earth_lights_lrg.jpg)

As imagens e o modelo NASA seguem as [diretrizes de uso de mídia da NASA](https://www.nasa.gov/nasa-brand-center/images-and-media/). O uso não implica endosso institucional. Não são imagens geradas por IA; o brilho do Sol e as estrelas decorativas são efeitos procedurais. O modelo é uma representação visual do MiRaTa, e não CAD de fabricação ou geometria documentada dos satélites da pesquisa.

O material de `original/` e o relatório preservam autoria e eventuais direitos de terceiros. O ZIP inclui bibliografia e código de referência; a publicação não atribui uma licença de software ao conjunto nem altera os avisos existentes.
