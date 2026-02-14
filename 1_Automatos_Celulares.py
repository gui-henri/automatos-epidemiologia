import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
import pandas as pd


st.title("Automatos Celulares")
st.write("Essa página servirá para criar meus experimentos com automatos celulares para as cadeiras de epidemiologia computacional e automatos celulares.")

# --- CONFIGURAÇÕES E CONSTANTES ---
S, E, I, R, V = 0, 1, 2, 3, 4
STATES = [S, E, I, R, V]
COLORS = {
    S: [240, 240, 240], # Branco/Cinza claro (Suscetível)
    E: [255, 165, 0],   # Laranja (Exposto)
    I: [255, 0, 0],     # Vermelho (Infectado)
    R: [0, 128, 0],     # Verde (Recuperado)
    V: [0, 0, 255]      # Azul (Vacinado)
}

st.set_page_config(page_title="Simulação SEIRV - Autômato Celular", layout="wide")

# --- FUNÇÕES NUCLEARES (CORE) ---

def init_grid(size, initial_count=5, location_mode="Aleatória"):
    """
    Inicializa o grid com células S e uma quantidade fixa de I em posições controladas.
    """
    grid = np.zeros(size, dtype=int) # Tudo começa como S (0)
    rows, cols = size
    
    # Garante que não pedimos mais infectados do que celulas existem
    initial_count = min(initial_count, rows * cols)
    
    if location_mode == "Aleatória":
        # Espalha aleatoriamente pelo mapa todo
        indices = np.random.choice(rows * cols, initial_count, replace=False)
        grid.flat[indices] = I
        
    elif location_mode == "Centro":
        # Cria um bloco compacto no centro
        center_r, center_c = rows // 2, cols // 2
        # Raio aproximado para caber a quantidade
        radius = int(np.sqrt(initial_count)) + 1
        
        count = 0
        for r in range(center_r - radius, center_r + radius + 1):
            for c in range(center_c - radius, center_c + radius + 1):
                if 0 <= r < rows and 0 <= c < cols and count < initial_count:
                    grid[r, c] = I
                    count += 1
                    
    elif location_mode == "Canto (Q1)":
        # Inicia no canto superior esquerdo (ideal para testar barreiras)
        # Começa um pouco afastado da borda (offset 5) para visualização
        start_r, start_c = 5, 5
        radius = int(np.sqrt(initial_count)) + 1
        
        count = 0
        for r in range(start_r, start_r + radius + 1):
            for c in range(start_c, start_c + radius + 1):
                if 0 <= r < rows and 0 <= c < cols and count < initial_count:
                    grid[r, c] = I
                    count += 1
                    
    return grid

def count_neighbors(grid, state):
    """Conta vizinhos em um estado específico usando convolução (Moore Neighborhood)."""
    mask = (grid == state).astype(int)
    kernel = np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]])
    return convolve2d(mask, kernel, mode='same', boundary='fill', fillvalue=0)

def apply_vaccination_strategies(grid, strategy, vac_params):
    """Aplica vacinação estática (Aleatória ou Cluster)."""
    new_grid = grid.copy()
    rows, cols = grid.shape
    
    if strategy == "Aleatória":
        # Vacina % das células S aleatoriamente
        rate = vac_params['random_rate']
        mask_s = (grid == S)
        random_v = np.random.rand(rows, cols) < rate
        new_grid[mask_s & random_v] = V
        
    elif strategy == "Cluster":
        # --- ESTRATÉGIA DE CORTINA DE FOGO (BARREIRAS GEOGRÁFICAS) ---
        # Cria barreiras em cruz para dividir o grid em 4 quadrantes isolados.
        # Isso simula cordões sanitários ou barreiras naturais (como rios/montanhas).
        
        # 'cluster_size' aqui define a ESPESSURA da parede (proporcional ao tamanho do grid)
        # Um valor baixo (ex: 0.05) cria uma barreira fina, porém transponível se houver "falhas" (mas aqui é sólida).
        thickness_factor = vac_params['cluster_size'] 
        
        # Calcula a espessura em pixels (mínimo de 1 pixel)
        thickness = max(1, int(min(rows, cols) * thickness_factor * 0.2)) 
        
        center_r, center_c = rows // 2, cols // 2
        
        # Define limites da barreira Horizontal (divide Norte/Sul)
        r_start = max(0, center_r - thickness)
        r_end = min(rows, center_r + thickness + 1)
        
        # Define limites da barreira Vertical (divide Leste/Oeste)
        c_start = max(0, center_c - thickness)
        c_end = min(cols, center_c + thickness + 1)
        
        # Máscara para a Cruz (+)
        # 1. Aplica vacina na faixa horizontal inteira
        mask_h = np.zeros_like(grid, dtype=bool)
        mask_h[r_start:r_end, :] = True
        
        # 2. Aplica vacina na faixa vertical inteira
        mask_v = np.zeros_like(grid, dtype=bool)
        mask_v[:, c_start:c_end] = True
        
        # Combina as duas barreiras
        barrier_mask = mask_h | mask_v
        
        # Aplica a vacinação: Onde é Suscetível (S) E está na barreira, vira Vacinado (V)
        # Isso cria a "topografia" intransponível para o vírus.
        target_cells = (grid == S) & barrier_mask
        new_grid[target_cells] = V
        
    return new_grid

def step(grid, params, vac_strategy, vac_params):
    """Executa um passo temporal da simulação."""
    # Extrai parâmetros
    beta = params['beta']   # Infecção (S->E)
    sigma = params['sigma'] # Incubação (E->I)
    gamma = params['gamma'] # Recuperação (I->R)
    
    # Matrizes de vizinhos infectados
    neighbors_i = count_neighbors(grid, I)
    
    # Gera números aleatórios para transições
    rand_matrix = np.random.rand(*grid.shape)
    
    new_grid = grid.copy()
    
    # Máscaras de estado atual
    is_S = (grid == S)
    is_E = (grid == E)
    is_I = (grid == I)
    
    # --- Lógica de Vacinação Reativa (Dinâmica) ---
    # Acontece ANTES da infecção neste passo
    became_v = np.zeros_like(grid, dtype=bool)
    if vac_strategy == "Reativa":
        rate = vac_params['reactive_rate']
        # Se vizinho é I e eu sou S -> chance de virar V
        prob_vac = 1 - (1 - rate) ** neighbors_i
        mask_vac = is_S & (rand_matrix < prob_vac)
        new_grid[mask_vac] = V
        became_v = mask_vac

    # --- Lógica SEIR ---
    
    # S -> E (Infecção)
    # Probabilidade de NÃO ser infectado = (1 - beta)^vizinhos_infectados
    # Probabilidade de SER infectado = 1 - (1 - beta)^vizinhos_infectados
    prob_infection = 1 - (1 - beta) ** neighbors_i
    # Só infecta se era S e não virou V neste turno
    mask_infection = is_S & (~became_v) & (rand_matrix < prob_infection)
    new_grid[mask_infection] = E
    
    # E -> I (Incubação)
    mask_incub = is_E & (rand_matrix < sigma)
    new_grid[mask_incub] = I
    
    # I -> R (Recuperação)
    mask_recov = is_I & (rand_matrix < gamma)
    new_grid[mask_recov] = R
    
    return new_grid

def box_counting(grid, target_state=I):
    """Calcula a dimensão fractal (Box-Counting) para o estado alvo."""
    # Binariza o grid: 1 onde tem o estado alvo, 0 caso contrário
    pixels = (grid == target_state)
    
    if np.sum(pixels) == 0:
        return 0, [] # Sem dimensão se vazio
        
    Lx, Ly = pixels.shape
    L = min(Lx, Ly)
    
    # Tamanhos das caixas (potências de 2)
    scales = np.logspace(np.log10(1), np.log10(L/2), num=10, endpoint=True, base=10).astype(int)
    scales = np.unique(scales) # Remove duplicatas se grid for pequeno
    scales = scales[scales > 0]
    
    counts = []
    
    for scale in scales:
        # Soma blocos de tamanho 'scale'
        # Truque rápido: reshape para dividir em blocos e somar eixos
        # Funciona perfeitamente apenas se grid for múltiplo, então usamos uma aproximação de 'janela deslizante' não sobreposta
        # Para ser robusto e rápido, usamos histograma 2d ou strided tricks, 
        # mas aqui faremos um loop simples otimizado nas fatias
        
        ns = 0
        # Grid reduzido (coarse-graining)
        for x in range(0, Lx, scale):
            for y in range(0, Ly, scale):
                if np.any(pixels[x:x+scale, y:y+scale]):
                    ns += 1
        counts.append(ns)
    
    # Ajuste linear (Regressão) log(N) vs log(1/s)
    # D = - slope
    coeffs = np.polyfit(np.log(scales), np.log(counts), 1)
    dimension = -coeffs[0]
    
    return dimension, list(zip(scales, counts))

def render_grid(grid):
    """Converte grid numérico para imagem RGB."""
    rgb_grid = np.zeros((*grid.shape, 3), dtype=np.uint8)
    for state, color in COLORS.items():
        rgb_grid[grid == state] = color
    return rgb_grid

# --- INTERFACE STREAMLIT ---

st.title("🔬 Modelo Epidemiológico SEIRV 2D + Dimensão Fractal")
st.markdown("""
Este modelo simula a propagação usando um Autômato Celular. 
A **Dimensão de Box-Counting** mede a complexidade espacial dos casos infectados ($I$).
""")

# Sidebar
st.sidebar.header("Parâmetros do Modelo")

st.sidebar.subheader("Configuração da Infecção (Paciente Zero)")

# Slider para quantidade exata (1 a 500 pessoas)
initial_infected_count = st.sidebar.number_input(
    "Qtd. Inicial de Infectados", 
    min_value=1, 
    max_value=1000, 
    value=5,
    step=1
)

# Seletor de Onde Começa
infection_location = st.sidebar.selectbox(
    "Local do Foco Inicial",
    ["Aleatória", "Centro", "Canto (Q1)"],
    index=2, # Padrão sugerido: Canto
    help="Use 'Canto (Q1)' junto com Vacinação 'Cluster' para ver a barreira funcionando."
)
grid_size = st.sidebar.slider("Tamanho do Grid", 50, 200, 100)
n_steps = st.sidebar.slider("Passos de Tempo", 10, 300, 100)

st.sidebar.subheader("Taxas (Probabilidades)")
beta = st.sidebar.slider("Taxa de Infecção (Beta)", 0.0, 1.0, 0.3, help="Probabilidade de S infectar-se por 1 vizinho I")
sigma = st.sidebar.slider("Taxa de Incubação (Sigma)", 0.0, 1.0, 0.2, help="Probabilidade de E virar I")
gamma = st.sidebar.slider("Taxa de Recuperação (Gamma)", 0.0, 1.0, 0.1, help="Probabilidade de I virar R")

st.sidebar.subheader("Estratégia de Vacinação")
vac_strategy = st.sidebar.selectbox("Tipo de Vacinação", ["Nenhuma", "Aleatória", "Cluster", "Reativa"])

vac_params = {'random_rate': 0.0, 'cluster_size': 0.0, 'reactive_rate': 0.0}

if vac_strategy == "Aleatória":
    vac_params['random_rate'] = st.sidebar.slider("% da População Vacinada (Inicial)", 0.0, 1.0, 0.1)
elif vac_strategy == "Cluster":
    vac_params['cluster_size'] = st.sidebar.slider("Tamanho do Cluster (% da largura)", 0.1, 1.0, 0.3)
elif vac_strategy == "Reativa":
    vac_params['reactive_rate'] = st.sidebar.slider("Eficácia da Reação Local", 0.0, 1.0, 0.5, help="Chance de S virar V se tiver vizinho I")

# Botão de Execução
if st.button("Executar Simulação"):
    # Setup
    params = {'beta': beta, 'sigma': sigma, 'gamma': gamma}
    grid = init_grid((grid_size, grid_size), initial_count=initial_infected_count, location_mode=infection_location)
    
    # Aplica vacinação estática inicial
    if vac_strategy in ["Aleatória", "Cluster"]:
        grid = apply_vaccination_strategies(grid, vac_strategy, vac_params)
    
    # Armazenamento de Histórico
    history_counts = {k: [] for k in STATES}
    history_fractal = []
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Placeholders para visualização ao vivo
    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("Visualização do Grid")
        grid_placeholder = st.empty()
    with col2:
        st.subheader("Curvas SEIRV")
        chart_placeholder = st.empty()
        st.subheader("Dimensão Fractal (Infectados)")
        fractal_placeholder = st.empty()

    # Loop Principal
    for t in range(n_steps):
        # Passo do AC
        grid = step(grid, params, vac_strategy, vac_params)
        
        # Coleta Dados
        unique, counts = np.unique(grid, return_counts=True)
        counts_dict = dict(zip(unique, counts))
        for s in STATES:
            history_counts[s].append(counts_dict.get(s, 0))
            
        # Calcula Fractal a cada X passos para performance (ou todo passo se grid pequeno)
        if t % 2 == 0 or t == n_steps - 1:
            fd, _ = box_counting(grid, target_state=I)
            history_fractal.append(fd)
        else:
            history_fractal.append(history_fractal[-1] if history_fractal else 0)

        # Atualiza Interface a cada X frames para não travar navegador
        if t % 5 == 0 or t == n_steps - 1:
            # Imagem
            img = render_grid(grid)
            grid_placeholder.image(img, caption=f"Passo {t}", output_format="PNG", width="stretch")
            
            # Gráfico Linhas
            df_hist = pd.DataFrame(history_counts)
            df_hist.columns = ['S', 'E', 'I', 'R', 'V']
            chart_placeholder.line_chart(df_hist)
            
            # Gráfico Fractal
            fractal_placeholder.line_chart(history_fractal)
            
        progress_bar.progress((t + 1) / n_steps)
        status_text.text(f"Simulando passo {t+1}/{n_steps}...")

    status_text.text("Simulação Concluída!")
    
    # Resultados Finais e Explicação
    st.markdown("---")
    st.subheader("Interpretação da Dimensão Fractal")
    st.info("""
    A dimensão fractal (Box-Counting) dos infectados indica como o vírus preenche o espaço:
    - **D próximo de 0:** Infecção extinta ou muito esparsa.
    - **D próximo de 1:** Frentes de onda lineares ou caminhos finos de transmissão.
    - **D próximo de 2:** O vírus ocupou áreas densas ("manchas" sólidas de infecção).
    
    Oscilações em D sugerem que o vírus está rompendo barreiras de imunidade ou criando novos focos (clusters) isolados.
    """)