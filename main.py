import time
import random
import matplotlib.pyplot as plt
import csv
import os
import copy

def guloso_vmp(D, T, I, a):
    distancia = 0
    ini = a
    rota = [a] 
    I.remove(a)  

    while len(I) > 0:
        x = 100**10
        k = a
        for j in I:
            if j != a:
                x = min(x, D[a][j])
                k = j
        rota.append(k) 
        I.remove(k)  
        distancia += D[a][k] 
        a = k 

    distancia += D[a][ini]  
    rota.append(ini)  

    tempo = calcula_tempo_da_rota(rota, T)

    return rota, distancia, tempo

def calcula_tempo_da_rota(rota, T):
    tempo = 0
    for i in range(len(rota) - 1):
        tempo += T[rota[i]][rota[i + 1]] 
    return tempo

def guloso_vmp_d_t(D, T, I):
    distancia = 0
    tempo = 0
    a = random.choice(I)  
    ini = a
    rota = [a]  
    I.remove(a)  

    while len(I) > 0:
        r = random.random() 
        x = 100**10
        k = a
        for j in I:
            if j != a:
                if r < 0.5:
                    x = min(x, D[a][j])
                else: 
                    x = min(x, T[a][j])
                k = j
        
        rota.append(k)  
        I.remove(k)  
        distancia += D[a][k] 
        tempo += T[a][k]  
        a = k 

    distancia += D[a][ini]
    tempo += T[a][ini]
    rota.append(ini)  

    return rota, distancia, tempo

def gerar_conjunto_s(D, T, I):
    n = len(D)
    S = []

    # Gerar n soluções minimizando distância
    for i in I:
        _, distancia, tempo = guloso_vmp(D, T, copy.deepcopy(I), i)
        S.append((distancia, tempo))

    # Gerar n soluções minimizando tempo
    for i in I:
        _, distancia, tempo = guloso_vmp(D, T, copy.deepcopy(I), i)
        S.append((distancia, tempo))

    # Gerar 6n soluções com algoritmo guloso VMP_D_T
    for _ in range(6 * n):
        _, distancia, tempo = guloso_vmp_d_t(D, T, copy.deepcopy(I))
        S.append((distancia, tempo))

    # Gerar 2n soluções totalmente aleatórias
    for _ in range(2 * n):
        rota = random.sample(I, n)
        distancia = sum(D[rota[i]][rota[i + 1]] for i in range(n - 1))
        distancia += D[rota[-1]][rota[0]]  # Fechar o ciclo
        tempo = sum(T[rota[i]][rota[i + 1]] for i in range(n - 1))
        tempo += T[rota[-1]][rota[0]]  # Fechar o ciclo
        S.append((distancia, tempo))

    return S

def dc(S):
    if len(S) == 1:
        return S
    else:
        S.sort(key=lambda x: x[0])  
        print("S: ", S)

        meio = len(S) // 2
        S1 = S[:meio]
        S2 = S[meio:]  

        M1 = dc(S1) 
        M2 = dc(S2)  

        M = combina(M1, M2)
        return M

def combina(M1, M2):
    M = (M1.copy() + M2.copy())
    M.sort()
    dominado = [True for _ in range(len(M))]

    d0,t0 = M[0]
    for i in range(1,len(M)):
        di,ti = M[i]
        if (d0 < di) and (t0 > ti):
            d0,t0 = di,ti
        else:
            dominado[i] = False

    M = [m for i, m in enumerate(M) if dominado[i]]
    return M

def open_arq(arq):
    with open(arq, 'r') as f:
        lines = f.readlines()
    
    dimension_line = next(line for line in lines if line.startswith("DIMENSION"))
    n = int(''.join(filter(str.isdigit, dimension_line.split(":")[1].strip())))

    start_d = lines.index("D = [\n") + 1
    end_d = start_d + n
    D = [list(map(int, line.strip().strip("[],").split(","))) for line in lines[start_d:end_d]]

    start_t = lines.index("T = [\n") + 1
    end_t = start_t + n
    T = [list(map(int, line.strip().strip("[],").split(","))) for line in lines[start_t:end_t]]

    I = list(range(n))

    return D, T, I

def processar_instancia(nome_arquivo):
    D, T, I = open_arq(nome_arquivo)
    
    inicio_tempo = time.time()
    
    S = gerar_conjunto_s(D, T, I)
    
    solucoes_minimais = dc(S)
    scatterplot_conjunto_s(S,solucoes_minimais, nome_arquivo)
    tempo_cpu = (time.time() - inicio_tempo) * 1000  # Convertendo para ms
    
    return {
        "solucoes_minimais": solucoes_minimais,
        "tempo_cpu": tempo_cpu
    }

def scatterplot_conjunto_s(S, M, nome_arquivo):
    os.makedirs("testes", exist_ok=True)
    # Separar os valores de dR (distância) e tR (tempo) para o conjunto S
    dR_S, tR_S = zip(*S)

    # Separar os valores de dR e tR para as soluções minimais M
    dR_M, tR_M = zip(*sorted(M, key=lambda x: x[0]))  # Ordenar M por dR para o lineplot

    # Criar o scatterplot do conjunto S
    plt.figure(figsize=(10, 6))
    plt.scatter(dR_S, tR_S, color='blue', label="Conjunto S (dR, tR)", alpha=0.6)

    # Adicionar o lineplot das soluções minimais
    plt.plot(dR_M, tR_M, color='red', linestyle='-', marker='o', label="Soluções Minimais (M)")

    # Adicionar destaques para o ponto de menor distância e menor tempo
    plt.scatter(dR_M[0], tR_M[0], color='green', s=100, label="Menor Distância")
    plt.scatter(dR_M[-1], tR_M[-1], color='purple', s=100, label="Menor Tempo")

    # Adicionar rótulos e título
    plt.xlabel("dR (Distância)", fontsize=12)
    plt.ylabel("tR (Tempo)", fontsize=12)
    plt.title("Conjunto S com Soluções Minimais", fontsize=14)
    plt.grid(True)

    # Adicionar legenda
    plt.legend()

    caminho_arquivo = os.path.join("testes", f"conjunto_s_{nome_arquivo}.png")
    plt.savefig(caminho_arquivo)
    plt.close()

def gerar_graficos_e_tabelas(tabelas, diretorio_saida):
    os.makedirs(diretorio_saida, exist_ok=True)

    with open(os.path.join(diretorio_saida, "tabela1.csv"), "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Instância", "No de Soluções Mínimas", "Menor dR", "Menor tR", "Tempo CPU (ms)"])
        for instancia, dados in tabelas.items():
            num_solucoes = len(dados["solucoes_minimais"])
            menor_dr = min(dados["solucoes_minimais"], key=lambda x: x[0])
            menor_tr = min(dados["solucoes_minimais"], key=lambda x: x[1])
            writer.writerow([instancia, num_solucoes, menor_dr, menor_tr, dados["tempo_cpu"]])

    with open(os.path.join(diretorio_saida, "tabela2.csv"), "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Instância", "Conjunto de Soluções Mínimas (Ordenadas por dR)"])
        for instancia, dados in tabelas.items():
            solucoes_ordenadas = sorted(dados["solucoes_minimais"], key=lambda x: x[0])
            writer.writerow([instancia] + [f"({dr}, {tr})" for dr, tr in solucoes_ordenadas])

    for numero_figura, (instancia, dados) in enumerate(tabelas.items(), start=1):
        solucoes_ordenadas = sorted(dados["solucoes_minimais"], key=lambda x: x[0])
        dR, tR = zip(*solucoes_ordenadas)
        
        # Identificar os índices dos pontos de mínimo e máximo
        min_idx = tR.index(min(tR))
        max_idx = tR.index(max(tR))
        
        plt.figure(figsize=(10, 5))
        
        # Plotar a curva principal
        plt.plot(dR, tR, 'o-', label="Soluções Minimais", color="blue")
        
        # Destacar o ponto de máximo
        plt.scatter(dR[max_idx], tR[max_idx], color='red', label="Máximo", zorder=5)
        plt.text(dR[max_idx], tR[max_idx], f"({dR[max_idx]}, {tR[max_idx]})", color="red", fontsize=9)
        
        # Destacar o ponto de mínimo
        plt.scatter(dR[min_idx], tR[min_idx], color='green', label="Mínimo", zorder=5)
        plt.text(dR[min_idx], tR[min_idx], f"({dR[min_idx]}, {tR[min_idx]})", color="green", fontsize=9)
    
        # Configurar o título e os rótulos dos eixos
        plt.title(f"Fig{numero_figura}. {instancia}")
        plt.xlabel("dR (Distância)")
        plt.ylabel("tR (Tempo)")
        
        # Adicionar grade e legenda
        plt.grid(True)
        plt.legend()
        
        # Salvar o gráfico no diretório especificado
        plt.savefig(os.path.join(diretorio_saida, f"grafico_{instancia}.png"))
        plt.close()

if __name__ == "__main__":
    arquivos = ['kroAxkroB100','kroAxkroC100','kroAxkroD100','kroBxkroC100','kroBxkroD100','kroCxkroD100','kroAxkroB150']
    tabelas = {}
    for nome in arquivos:
        tabelas[nome] = processar_instancia(nome+'.txt')
    gerar_graficos_e_tabelas(tabelas, "resultados_saida")