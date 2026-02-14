# Simulação Epidemiológica SEIRV 2D via Autômatos Celulares
Este projeto é uma aplicação web interativa desenvolvida com Streamlit para simular a propagação de uma doença infecciosa utilizando um Autômato Celular (AC) bidimensional. O modelo baseia-se na dinâmica compartimental SEIRV (Suscetível, Exposto, Infectado, Recuperado, Vacinado) e introduz a análise em tempo real da Dimensão Fractal (pelo método de Box-Counting) para quantificar a complexidade espacial do surto.

## 🎯 Objetivo
A ferramenta foi projetada para apoiar experimentos práticos em epidemiologia computacional, permitindo a visualização de como diferentes parâmetros de contágio e estratégias espaciais de vacinação afetam a curva epidêmica e a morfologia da rede de infectados ao longo do tempo.

## 🦠 Dinâmica do Modelo (SEIRV)
A grade de simulação consiste em uma malha $L \times L$, onde cada célula representa um indivíduo que pode assumir um dos cinco estados $S, E, I, R, V$. As transições de estado ocorrem de forma estocástica e são influenciadas pela vizinhança local (Vizinhança de Moore - 8 vizinhos adjacentes).
As principais transições matemáticas definidas no modelo são:
- Infecção ($S \rightarrow E$): Um indivíduo suscetível contrai a doença baseado na taxa de infecção ($\beta$) e no número de vizinhos infectados ($n$). A probabilidade de infecção é calculada por: $$P_{inf} = 1 - (1 - \beta)^n$$
- Incubação ($E \rightarrow I$): A passagem do estado latente para o estado infeccioso ocorre com probabilidade $\sigma$.
- Recuperação ($I \rightarrow R$): A transição para o estado imune natural ocorre com probabilidade $\gamma$.
## 💉 Estratégias de Vacinação

O modelo permite testar intervenções farmacológicas simulando diferentes formas de alocar vacinas (estado $V$) no espaço geométrico da população:
- Nenhuma: O vírus circula livremente sem intervenção externa.
- Aleatória: Simula campanhas de vacinação em massa não direcionadas. Uma porcentagem inicial da população suscetível é imunizada aleatoriamente no passo $t=0$.
- Cluster (Barreira): Cria cordões sanitários geográficos em forma de cruz, dividindo a malha em quatro quadrantes. O objetivo é testar a eficácia de compartimentar o surto espacialmente.
- Reativa: Simula o rastreamento de contatos dinâmico. A cada passo de tempo, suscetíveis que estão na vizinhança de um infectado têm uma chance $\eta$ de serem vacinados rapidamente, bloqueando a expansão de clusters locais.
-
- ## 📐 Dimensão Fractal (Box-Counting)
-
- Um diferencial deste modelo é o cálculo contínuo da dimensão fractal ($D_S$) da massa de infectados. O algoritmo sobrepõe grades de diferentes tamanhos ($\epsilon$) e conta o número de caixas $N(\epsilon)$ que contêm pelo menos uma célula infectada, derivando a dimensão a partir da relação:$$N(\epsilon) \propto \epsilon^{-D_S}$$Interpretação:$D_S \approx 0$: Epidemia contida, focos isolados e esparsos.$D_S \approx 1$: Propagação linear, avançando como filamentos ou finas frentes de onda.$D_S \approx 2$: Propagação densa, preenchendo maciçamente áreas inteiras do espaço geográfico.
-
- ## ⚙️ Instalação e Execução
- ### Pré-requisitos
- Certifique-se de ter o Python 3.8+ instalado. Recomenda-se o uso de um ambiente virtual (venv).
- ### Instalação das dependências
- Abra o terminal e instale as bibliotecas necessárias:
```Bash
  pip install streamlit numpy scipy pandas matplotlib
```
## Executando a aplicação

Navegue até o diretório onde o arquivo do código (ex: app.py) está salvo e execute:
```bash
streamlit run app.py
```
Isso abrirá automaticamente uma nova aba no seu navegador padrão (geralmente em http://localhost:8501) contendo a interface da simulação.
