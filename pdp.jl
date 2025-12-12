using HiGHS
using JuMP
# using Graphs
# using GraphPlot
using Plots
using LinearAlgebra

mutable struct Label

    C :: Float64
    R :: Vector{Float64}
    s :: Integer
    V :: Vector{Float64}

    Label(C, R, V) = new(C, copy(R), sum(0 .< V), copy(V))

end

Label(L, n) = Label(Inf, zeros(Int, L), zeros(Int, n))

import Base

"""
Essa funcao implementa os rotulos dominados, ou seja, `x ≤ y` se e somente se `y` for dominado por `x`
"""

function Base.:<=(x :: Label, y :: Label) # essa funcao verifica quem do label e' <= (verifica qual rotulo e dominado pelo outro)

    return (x.C <= y.C) && all(x.R .<= y.R) && (x.s <= y.s) && all(p -> (p[1] == 0) || (p[2] > 0), zip(x.V, y.V))

end

function Base.isequal(x :: Label, y :: Label) # essa funcao verifica que o label e' == (verifica a igualdade entre os dois rotulos)

    # Atencao ao uso de == para numeros reais!!
    return (x.C == y.C) && all(x.R .== y.R) && (x.s == y.s) && all(x.V .== y.V)

end

function Base.copy(λ :: Label) # recebe um objeto λ, onde cria um novo objeto com os mesmos valores de C, R, V e retorna uma nova instancia

    return Label(λ.C, λ.R, λ.V)

end

# funcoes auxiliares
fromvec(v) = Label(v[1], [v[3]], v[4:end])
tovec(λ) = [λ.C; λ.s; λ.R...; λ.V...]

"""
Funcionalidade: aplica o algoritmo de correção de rótulos em um gráfico totalmente conectado para o ESPPRC.


### Entrada:
V - conjunto de vértices do gráfico;
S - matriz nxn onde cada entrada representa o valor objetivo obtido ao viajar ao longo de cada aresta do gráfico;
T - matriz nxn onde T[i,j] corresponde ao tempo usado para viajar do destino de i para a origem de j + tempo da
origem de j para o destino de j;
L - Float64 correspondendo ao limite superior para o recurso de tempo;
W - matriz nx2 contendo a janela de tempo de cada vértice;
Lambda - matriz nx1 onde cada entrada é o conjunto de rótulos do vértice correspondente;
n - número de vértices do gráfico.

### Saida:
   Lambda - matriz nx1 onde cada entrada é o conjunto de rótulos do vértice correspondente;
    bests - valor da função objetivo ótima (o rótulo inteiro)

"""
function ESPPRC(S, T, L, W, Λ, n)
    
    bests = Label(1, n)
    
    i = 1
    
    j = 1
    
    println("ESPPRC: L=$L")

    println("\tS = ", S)
    println("\tT = ", T)

    E = zeros(n)

    E[1] = 1

    F = Set([])
    
    changed = false
    
    #println("ESPPRC: n = $n")

    while sum(E) != 0
        
        if E[i] != 0

            #println("EPPRC: Analisando $i")
                
            for k in 0:(n-2) # i é o vértice que estamos tratando e k varia ao longo de 0 a n-2
                
                j = (i+k) % n + 1  # j é um vértice ao qual podemos estender uma rótula
                                
                for λ in Λ[i]

                    if λ.V[j] == 0 # se a extensão for possível, fazemos isso

                        #println("ESPPRC: extensao para = $j")
                        
                        F = union(F, Set([Extend(λ, i, j, S, T, L, W)]))

                        #for f in F

                            #println("\t $f")

                        #end

                    end
                    
                end
                
                # println("Vertex $j")
                Λ[j], changed, bests = EFF(Λ[j], F, bests)

                F = Set([])

                if changed
                    
                    E[j] = 1
                    
                end
                
            end
            
            E[i] = 0
            
        end
        
        i = i % n + 1
        
    end

    return Λ, bests

end

"""
Funcionalidade: estender um rótulo do vértice i ao vértice j (por certificação em ESPPRC(), a extensão sempre será possível)

### Entrada:
lambda_i - rótulo da matriz (2+n)x1 a ser estendido de i;
i - vértice que contém o rótulo atual;
j - vértice que receberá um novo rótulo;
T - matriz nxn onde T[i,j] corresponde ao tempo usado para viajar do destino de i para a origem de j + tempo da
origem de j para o destino de j;
L - Float64 correspondendo ao limite superior para o recurso de tempo;
W - matriz nx2 contendo a janela de tempo de cada vértice.

### Saída:
lambda_j - matriz (2+n)x1 (rótulo estendido).
"""
function Extend(λ_i, i, j, S, T, L, W)
    
    n = size(T, 1)

    λ_j = copy(λ_i)

    λ_j.R[1] = max(W[j, 1], λ_j.R[1] + T[i,j]) # atualizamos o tempo percorrido e marcamos j como inacessível se lambda_j[2] + T[i,j] < W[j, 1] significa

    # que chegamos na tarefa antes do início da janela de tempo de execução da tarefa, ou seja, devemos esperar
    
    # atualizando a ordem em que a tarefa j foi concluída
    
    order = 1
    
    for k in λ_j.V
        if (k != Inf) && (k > order)
            order = k
        end
    end
    
    λ_j.V[j] = order + 1

    # atualizando o uso de recursos
    
    # a primeira entrada corresponde ao valor da função objetivo (número de tarefas concluídas)

    λ_j.C += S[i,j]
    
    # atualizando vértices inalcançáveis ​​de j agora  
    
    for k in 1:n
        
        #println("\tEXTEND: V[$k]=$(λ_j.V[k])")
        #println("\tEXTEND: $k - $((λ_j.R[1] + T[j,k] > W[k,2])), $((λ_j.R[1] + T[j,k] < W[k,1] - 5)), $(((λ_j.R[1] + T[j,k]) > L))")
        if λ_j.V[k] == 0

            # Levamos em consideração o tempo que um veículo teria que esperar para concluir duas tarefas consecutivas
            
            if (λ_j.R[1] + T[j,k] > W[k,2]) || (λ_j.R[1] + T[j,k] < W[k,1] - 5) || ((λ_j.R[1] + T[j,k]) > L) # Aqui, 5 é o limite imposto para o tempo de marcha lenta do veículo

                λ_j.V[k] = Inf
                
            end

        end
        
    end
    
    λ_j.s = sum(λ_j.V .> 0) # medindo vértices inacessíveis (visitados, impossíveis de visitar ou não tão bons para visitar)  
        
    return λ_j
    
end

"""
Funcionalidade: mantém os conjuntos Lambda totalmente incomparáveis. Verifica qual rótulo pode entrar e aplicar a regra de dominação para remover
rótulos do conjunto.

### Entrada:
Lambda - conjunto de rótulos de um vértice;
F - conjunto de rótulos recém-estendidos (que possivelmente entrarão no Lambda);
bests - melhor valor de função objetivo (o rótulo inteiro) obtido até agora.

### Saída:
tmplambda - define Lambda possivelmente alterado;
changed - booleano que armazena se Lambda foi alterado;
bests - melhor valor de função objetivo (o rótulo inteiro) obtido até agora.
"""
function EFF(Λ :: Set, F :: Set, bests)
        
    dominated = false
    
    changed = false
    
    tmplambda = copy(Λ)
    
    if isempty(Λ) && !isempty(F)
        
        changed = true 
        
        for λ in F
            if λ.C < bests.C
                bests = λ
            end
        end
        
        return F, changed, bests
        
    else
        
        if isempty(F)
            
            changed = false
            
            return Λ, changed, bests
            
        else
            
            for λ_f in F
                
                for λ in tmplambda 

                    if λ <= λ_f && λ != λ_f
                        
                        dominated = true
                        
                        #println("Dominated $(λ_f) by $(λ)")
                        
                    end
                    
                    
                    if (λ == λ_f) || dominated
                        
                        break
                        
                    end
                    
                    
                end
                
                if !dominated && !(λ_f in tmplambda)
                    
                    changed = true
                                        
                    for λ in tmplambda

                        if λ_f <= λ
                            
                            dominated = true 
                            
                            #println("Dominated $(λ) by $(λ_f)")
                            
                            if dominated
                            
                                tmplambda = setdiff(tmplambda, Set([λ]))
                                
                            end
                            
                        end
                        
                    end
                    
                    tmplambda = union!(tmplambda, Set([λ_f]))
                                        
                    if λ_f.C < bests.C
                        
                        bests = λ_f
                    
                    end
                    
                end
                
                dominated = false
            
            end
            
        end
        
    end
    
    return tmplambda, changed, bests
    
end

function simplexrev(A,b,c,base)

    k=0 # contador de iterações

    (m,n)=size(A) # dimensões da matriz do problema
    solucao=vec(zeros(n)) # criando um vetor com zeros para receber a solução

    while k<1000

        B=A[:,base] # obtendo a matriz base

        cb=c[base] # vetor custo das variáveis básicas

        xb=B\b # atualizando xb

        y=B' \ cb # atualizando y, B'y = cb

        # encontrando o vetor com as colunas das variáveis não-básicas
        vnb=vec([1:n;]) # criando um vetor com o número de colunas de 1 a n
        for i=1:length(base) # os correspondentes da base recebem valor 0
           vnb[base[i]]=0
        end
        vnb=findall(x -> x>0,vnb) # os que não são zero (são positivos) formam vnb
        cnb=c[vnb] # cnb é a parte de c correspondente as variáveis não básicas

        D=A[:,vnb] # matriz com as colunas não básicas

        rn =  D'*y - cnb # rn = (zn - cn)' # custo relativo das variáveis não básicas

        

        # verificando a otimalidade
        if maximum(rn)<=0
            solucao[base]=xb
            votimo=c'*solucao   # ou cb'*xb

            # Cálculo do custo reduzido final
            zj = (c[base]' * (B \ A))'  # Zj = CB * B⁻¹ * Aj
            custo_reduzido = zj - c  # Zj - Cj
            
            println("Solução ótima=",solucao)
            println("Valor ótimo=",votimo)
            println("Iterações=",k)
            println("Base=",base)
            println("Zj - Cj (Custo reduzido final) = ", custo_reduzido) # Exibir custo reduzido na solução ótima
            return solucao, custo_reduzido
            
        end

        # encontrar quem entra na base
        aux=findmax(rn); # vnb(aux) que entrará na base
        ientra=vnb[aux[2]];
        colentra=A[:,ientra]; # selecionando a coluna a entrar na base A[:,ientra] ou D[:,aux]

        # y = inv(B)*colentra
        yentra = B\colentra # resolvendo o sistema B*yentra = colentra
        auxentra=findall(x -> x>0,yentra) # obtendo os indices do vetor y que são positivos

        if length(auxentra)==0  # se nenhum y_i é positivo o problema é ilimitado
            println("Problema Ilimitado")
            return
        end

        divi=xb[auxentra]./yentra[auxentra]  # Cálculo do quociente xb(i)/y(i) para y(i)>0
        aux1=findmin(divi) # o índice no qual o mínimo ocorre,aux1[2], é que indicará quem sai da base

        isai=base[auxentra[aux1[2]]] #base[auxentra[aux1[2]]] deve sair da base
        posisai=findall(x -> x==isai,base)
        base[posisai[1]]=ientra # atualizando os índices das colunas que formam a base

        k=k+1 # atualizando o contador
    end

end

"""
Funcionalidade: cria dados de acordo com parâmetros para o ESPPRC.

### Entrada:
n - número de tarefas(certificados) do gráfico (incluindo origem e destino);
L - limite superior do Float64 para o recurso de tempo;
dual - matriz 4x1 incluindo variáveis ​​duais (em ordem: lambda_0, lambda_1, lambda_2 E lambda_3)
C - matriz mxm correspondente à matriz de distância entre as m cidades;
task - matriz(n-2)x2 onde a primeira coluna corresponde à origem e a segunda coluna ao destino da tarefa representada pela i-ésima linha
W - matriz nx2 com intervalo de execução de cada tarefa

### Saída:
S, T, L, W, Lambda, n - saída que será a entrada correspondente do ESPPRC (verifique a documentação do ESPPRC)
"""
function data(n, L, dual, C, task)

    #println("Data")
 
    T = zeros(n,n)

    for i in 1:n

        T[i,i] = 2 * L

    end

    T[2:n-1,1] .= 2 * L

    T[n,:] .= 2 * L

    T[1:n-1,n] .= 0.0

    for i in 2:n-1

        for j in 2:n-1 

            if i!=j

                T[i,j] = C[task[i-1,2], task[j-1,1]] + C[task[j-1,1], task[j-1,2]]

            end

        end

    end

    for j in 2:n-1 

        # T[1,j] = C[task[1,1], task[j-1,1]] + C[task[j-1,1], task[j-1,2]]
        T[1,j] = C[task[j-1,1], task[j-1,2]]

    end

    #display(T)

    Lambda = Array{Any,1}(undef, n)

    label_origin = Label(0, [0.0], zeros(n))

    label_origin.V[1] = 1.0

    Lambda[1] = Set([label_origin])

    for i in 2:n
        Lambda[i] = Set{Label}([])
    end

    # Construcao da matriz de custos para o ESPPRC, considerando 1 como o vertice
    # artifical da origem e n como o vertice artificial do destino. Observe que nao
    # ha estimulo para fazer o vertice n durante o processo, pela forma que foi
    # construida a matriz T. Da mesma forma, nao ha estimulo para uma tarefa ser
    # feita novamente. Utilizamos 0 como um custo desestimulante, dado que queremos
    # minimizar a funcao. (Sera que 1 ou n seria um melhor valor?)

    S = Array{Float64}(undef, n, n)

    S[:, 1] .= 0

    S[:, end] .= 0

    S[1, end] = 0

    S[end, :] .= 0

    for i in 2:n - 1

        for j in i + 1:n - 1

            S[i,j]  = - 1 - dual[j]
            S[j, i] = - 1 - dual[i]

        end

        S[i, i] = 0

    end
    
    for j in 2:n-1

        S[1,j] = -1 - dual[j] - dual[1]

    end

    return S, T, Lambda
end

function tem_caminho_negativo(rotulos)
    for lambda in rotulos
        if lambda.C < 0
            return true
        end
    end
    return false
end

function A_final(C, task, Wu, num_caminhoes, r; MAXIT=10)
    rotulos = [[Label(-1.0, [0.0], [0])]]

    L = maximum(Wu[:, 2]) # correspondendo ao limite superior para o recurso de tempo
    #display(L)
    #display(C)
    println("matriz task")
    #display(task)
    println("matriz Wu")
    #display(Wu)
    #L = maximum(Wu[:, 2]) # correspondendo ao limite superior para o recurso de tempo

    n = r + 2 #criando varariavel artificial
     
    W = Array{Float64}(undef, n, 2)
    W[1,:] = [0, L]
    W[2:end-1,:] .= Wu
    W[end, :] = [0, L]
    println("matriz W")
    #display(W)
    #matriz A
    g = r + 1
    K = Matrix{Float64}(I, g, g)  # matriz identidade g x g
    # construindo R
    primeira_linha = ones(1, g - 1) # a primeira linha da matriz R vair ser uma matriz onde todos os elementos sao 1 de ordem (1, g-1) 
    
    parte_de_baixo = Matrix{Float64}(I, g - 1, g - 1) # a parte de baixo da matriz R vai ser a matriz Identidade de ordem (g-1, g-1)
    
    R = vcat(primeira_linha, parte_de_baixo) # usando o comando vcat pra concatenar as duas parte da matriz R em uma unica matriz
    
    # construindo A
    A = hcat(K, R) # usando o comando hcat pra concatenar a matriz K com a matriz R para formarmos a matriz A
    
    n_original = n
    
    b = vcat(num_caminhoes, ones(size(A,1)-1))

    base = collect(1:size(A,1)) # matriz da base 
    #o comando collect transforma um intervalo em vetor,matriz

    #L = maximum(W) # correspondendo ao limite superior para o recurso de tempo
    eh_otima = false
    salvaguarda = 1
    base_inicial = copy(base)
    solucao = [] 

    while !eh_otima && salvaguarda < MAXIT #tem_caminho_negativo(rotulos[end])

        println("""

        ==============
        Iteração $salvaguarda
        ==============
        """
        )

        salvaguarda += 1
        # Número de colunas e linhas (criar a matriz c)
        (m, n) = size(A)
        # Encontrar índices das colunas fora da base
        vnb = setdiff(collect(1:n), base)  # Colunas que não estão na base
        # Criar c automaticamente
        c = zeros(n)  # Inicializa um vetor de zeros com tamanho n
        for j in vnb
            c[j] = -sum(A[2:end, j])  # Soma dos elementos da coluna sem a primeira linha, com sinal negativo
        end

        #simplex 
        solucao, custo_reduzido = simplexrev(A,b,c,base)
        
        base.=base_inicial
        
        #ESPPRC

        var_dual = custo_reduzido[base_inicial]

        println("Custo arestas:\n\tλ0 = ", var_dual[1])
        println("\tλ = ", var_dual[2:end])

        S, T, Lambda = data(n_original, L, var_dual, C, task)
        # display(T)
        # display(L)
        # display(C)
        # display(S)
        rotulos, melhor = ESPPRC(S, T, L, W, Lambda, n_original)
        #println("ESPPRC devolveu")
        #for i in rotulos[end]
            #   println("\t$i")
        #end
        #println()

        eh_otima=true
        # Encontrar todos os lambdas com C < 0
        lambdas_negativos = [lambda for lambda in rotulos[end] if lambda.C < 0]

        println("\tRotulos a serem adicionados:")
        for λ in lambdas_negativos

            println("\t$λ")

        end

        if !isempty(lambdas_negativos)
            eh_otima = false

            valores_coluna = [lambda.V[1:end-1] for lambda in lambdas_negativos]
            valores_coluna = [map(x -> isinf(x) ? 0.0 : 1.0, v) for v in valores_coluna]

            # Adiciona a coluna 
            for nova_coluna in valores_coluna
            #     if all(!all(A[:, j] .== nova_coluna) for j in 1:size(A, 2))
                A = hcat(A, nova_coluna)
            #     end
            end

            #display(A)
        end

    end 

    return A, solucao, g

end

function rota_da_solucao(solucao, A, g)
    indices_nao_base = collect(g+1:length(solucao))
    indices_ativos = [j for j in indices_nao_base if solucao[j] != 0.0]
    colunas_ativas = A[2:end, indices_ativos]
    
    #println("Índices das colunas fora da base com x ≠ 0: ", indices_ativos)
    #println("Colunas correspondentes da matriz A:")
    #println(colunas_ativas)
    
    return  colunas_ativas
end