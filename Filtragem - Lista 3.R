#Filtragem Adaptativa
#Lista 3 - Algoritmos Recursivos do tipo LMS
#Autor:Felipe Pinto Marinho
#Data:14/02/2022

#Carregando algumas bibliotecas relevantes
library(ggplot2)
library(plotly)
library(pracma)
library(MASS)

#Geração dos dados
x = rnorm(301, mean = 0, sd = 1)
g = rep(0, 300)

for (i in 1:length(x)-1) {
  g[i] = x[i +1] - 1.6 * x[i]
}

dados = cbind(x[-length(x)], g)
dados = as.data.frame(dados)

#Representação gráfica dos sinais
figura_ruido = plot_ly(data = dados, x = 1:300, y = dados[, 1], type = "scatter", marker = list(size = 10, color = "rgba(125, 180, 240, .9)", line = list(color = "rgba(0, 152, 0, .8)", width = 2)), name = NULL) 
figura_ruido = figura_ruido %>% layout(xaxis = list(title = "Tempo discreto (k)", titlefont = list(size = 22), tickfont = list(size = 22)), yaxis = list(title = "Observação", titlefont = list(size = 22), tickfont = list(size = 22))) %>% add_lines(x = 1:300, y = dados[, 1], name = 'Ruído Gaussiano', line = list(width = 4))
figura_entrada = plot_ly(data = dados, x = 1:300, y = dados[, 2], type = "scatter", name = "Entrada do Filtro", line = list(color = "red", width = 4), marker = list(size = 10, color = "rgba(240, 18, 24, .9)", line = list(color = "rgba(152, 0, 0, .8)", width = 2), symbol = 'x'))
figura_entrada = figura_entrada %>% layout(xaxis = list(title = "Tempo discreto (k)", titlefont = list(size = 22), tickfont = list(size = 22)), yaxis = list(title = "Observação", titlefont = list(size = 22), tickfont = list(size = 22)))

subplot(figura_ruido, figura_entrada, nrows = 2, titleX = T, titleY = T, shareX = T)

#Aplicaçao dos algoritmos
d = x[-1]  #Problema de equalização o desejado deve ser o sinal anterior a passagem do canal
R = rbind(c(3.56, -1.6), c(-1.6, 3.56))
P = c(1, 1)

#Função para a determinação do MSE
MSE = function(R, P, variância, w){
  return(variância - 2 * dot(w, P) + t(w) %*% R %*% w)
}

#Implementação dos algoritmos recursivos do tipo LMS
#Gradiente Determinístico
Grad_Det = function(N, R, P, mu, x, d, ordem){
  
  #Inicialização da matriz de parâmetros e do vetor de erros quadráticos
  W = matrix(0, nrow = ordem, ncol = N)
  erro = rep(0, N)

  #Inicialização do vetor de parâmetros
  w0 = rep(0, ordem)
  w = w0
  
  #Laço de iteração
  for (k in 1:N) {
    
    #Determinação do erro
    e = d[k+(ordem-1)] - dot(x[k:(k+(ordem-1))], w)
    erro[k] = MSE(R, P, 1, w)
    
    #Atualização
    w = w - (2 * mu) * (-P + (R %*% w))
    
    #Armazena os coeficientes
    W[, k] = w
  }
  
  #Resultados
  resultados = list(W, erro)
  names(resultados) = c("W", "erro")
  return(resultados)
}

#Algoritmo de Newton
Newton = function(N, R, P, mu, x, d, ordem){
  
  #Inicialização da matriz de parâmetros e do vetor de erros quadráticos
  W = matrix(0, nrow = ordem, ncol = N)
  erro = rep(0, N)
  
  #Inicialização do vetor de parâmetros
  w0 = rep(0, ordem)
  w = w0
  
  #Laço de iteraçao
  for (k in 1:N) {
    
    #Determinação do erro
    e = d[k+(ordem-1)] - dot(x[k:(k+(ordem-1))], w)
    erro[k] = MSE(R, P, 1, w)
    
    #Atualização
    w = w - mu * (w - (ginv(R) %*% P))
    
    #Armazena os coeficientes
    W[, k] = w
    
  }
  
  #Resultados
  resultados = list(W, erro)
  names(resultados) = c("W", "erro")
  return(resultados)
}

#Algoritmo LMS
LMS = function(N, mu, x, d, ordem){
  
  #Inicialização da matriz de parâmetros e do vetor de erros quadráticos
  W = matrix(0, nrow = ordem, ncol = N)
  erro = rep(0, N)
  
  #Inicialização do vetor de parâmetros
  w0 = rep(0, ordem)
  w = w0
  
  #Laço de iteração
  for (k in 1:N) {
    
    #Determinaçao do erro
    e = d[k+(ordem-1)] - dot(x[k:(k+(ordem-1))], w)
    erro[k] = MSE(R, P, 1, w)
    
    #Atualização
    w = w + (2 * mu * e) * x[k:(k+(ordem-1))] 
    
    #Armazena os coeficientes
    W[, k] = w 
  }
  
  #Resultados
  resultados = list(W, erro)
  names(resultados) = c("W", "erro")
  return(resultados)
}

#Algortimo LMS normalizado
LMS_normalizado = function(N, mu, gama, x, d, ordem){
  
  #Inicialização da matriz de parâmetros e do vetor de erros quadráticos
  W = matrix(0, nrow = ordem, ncol = N)
  erro = rep(0, N)
  
  #Inicialização do vetor de parâmetros
  w0 = rep(0, ordem)
  w = w0
  
  #Laço de iteração
  for (k in 1:N) {
    
    #Determinaçao do erro
    e = d[k+(ordem-1)] - dot(x[k:(k+(ordem-1))], w)
    erro[k] = MSE(R, P, 1, w)
    
    #Atualização
    w = w + (mu * e ) * ((x[k:(k+(ordem-1))])/(gama + dot(x[k:(k+(ordem-1))], x[k:(k+(ordem-1))]))) 
    
    #Armazena os coeficientes
    W[, k] = w 
  }
  
  #Resultados
  resultados = list(W, erro)
  names(resultados) = c("W", "erro")
  return(resultados)
}


#Resultados
resultado_grad_det = Grad_Det(200, R, P, 0.2, g, d, 2)
resultado_newton = Newton(200, R, P, 0.2, g, d, 2)
resultado_LMS = LMS(200, 0.1, g, d, 2)
resultado_LMS_normalizado = LMS_normalizado(200, 0.1, 0.01, g, d, 2)

#Representação gráfica
#Gradiente Determinístico
erro = data.frame(resultado_grad_det$erro)
figura_erro_grad = plot_ly(data = erro, x = 1:200, y = erro[, 1], type = "scatter", marker = list(size = 10, color = "rgba(125, 180, 240, .9)", line = list(color = "rgba(0, 152, 0, .8)", width = 2)), name = NULL) 
figura_erro_grad = figura_erro_grad %>% layout(title = "MSE vs. Iteração", xaxis = list(title = "Iteração", titlefont = list(size = 22), tickfont = list(size = 22)), yaxis = list(title = "MSE", titlefont = list(size = 22), tickfont = list(size = 22))) %>% add_lines(x = 1:200, y = erro[,1], name = "Gradiente Determinístico", line = list(width = 4))
figura_erro_grad

#Newton
erro = data.frame(resultado_newton$erro)
figura_erro_newton = plot_ly(data = erro, x = 1:200, y = erro[, 1], type = "scatter", marker = list(size = 10, color = "rgba(240, 180, 140, .9)", line = list(color = "rgba(250, 0, 0, .8)", width = 2), symbol = 'x'), name = NULL) 
figura_erro_newton = figura_erro_newton %>% layout(title = "MSE vs. Iteração", xaxis = list(title = "Iteração", titlefont = list(size = 22), tickfont = list(size = 22)), yaxis = list(title = "MSE", titlefont = list(size = 22), tickfont = list(size = 22))) %>% add_lines(x = 1:200, y = erro[,1], name = "Newton", line = list(width = 4, color = "red"))
figura_erro_newton

#LMS
erro = data.frame(resultado_LMS$erro)
figura_erro_LMS = plot_ly(data = erro, x = 1:200, y = erro[, 1], type = "scatter", marker = list(size = 10, color = "rgba(140, 180, 240, .9)", line = list(color = "rgba(0, 0, 250, .8)", width = 2), symbol = 'cross'), name = NULL) 
figura_erro_LMS = figura_erro_LMS %>% layout(title = "MSE vs. Iteração", xaxis = list(title = "Iteração", titlefont = list(size = 22), tickfont = list(size = 22)), yaxis = list(title = "MSE", titlefont = list(size = 22), tickfont = list(size = 22))) %>% add_lines(x = 1:200, y = erro[,1], name = "LMS", line = list(width = 4, color = "blue"))
figura_erro_LMS

#LMS normamlizado
erro = data.frame(resultado_LMS_normalizado$erro)
figura_erro_LMS_normalizado = plot_ly(data = erro, x = 1:200, y = erro[, 1], type = "scatter", marker = list(size = 10, color = "rgba(40, 80, 40, .9)", line = list(color = "rgba(0, 0, 0, .8)", width = 2), symbol = 'triangle-down'), name = NULL) 
figura_erro_LMS_normalizado = figura_erro_LMS_normalizado %>% layout(title = "MSE vs. Iteração", xaxis = list(title = "Iteração", titlefont = list(size = 22), tickfont = list(size = 22)), yaxis = list(title = "MSE", titlefont = list(size = 22), tickfont = list(size = 22))) %>% add_lines(x = 1:200, y = erro[,1], name = "LMS normalizado", line = list(width = 4, color = "black"))
figura_erro_LMS_normalizado

subplot(figura_erro_grad, figura_erro_newton, figura_erro_LMS, figura_erro_LMS_normalizado, shareX = T, titleX = T, titleY = T, nrows = 2)

#Curvas de nível
w_1 = seq(-50, 50)
w_2 = seq(-50,50)
z = matrix(0, nrow = 101, ncol = 101)
for (i in 1:101) {
  for (j in 1:101) {
    
    z[i, j] = MSE(R,P, 1, c(w_1[i], w_2[j])) 
  }
}

par(mfrow = c(2,2))
w_opt = ginv(R) %*% P
cols = hcl.colors(10, "YlOrRd")

#Gradiente Determinístico
contour(w_1, w_2, z, main = "Gradiente Determinístico", xlab = "w1", ylab = "w2", col = cols)
points(x = w_opt[1], y = w_opt[2], cex = 1.5, pch = 19, col = "black")
points(x = resultado_grad_det$W[1,], y = resultado_grad_det$W[2,], cex = 1.5 ,col = "red")

#Newton
contour(w_1, w_2, z, main = "Algoritmo de Newton", xlab = "w1", ylab = "w2", col = cols)
points(x = w_opt[1], y = w_opt[2], cex = 1.5, pch = 19, col = "black")
points(x = resultado_newton$W[1,], y = resultado_newton$W[2,], cex = 1.5 ,col = "red")

#LMS
contour(w_1, w_2, z, main = "Algoritmo LMS", xlab = "w1", ylab = "w2", col = cols)
points(x = w_opt[1], y = w_opt[2], cex = 1.5, pch = 19, col = "black")
points(x = resultado_LMS$W[1,], y = resultado_LMS$W[2,], cex = 1.5 ,col = "red")

#LMS Normalizado
contour(w_1, w_2, z, main = "Algoritmo LMS Normalizado", xlab = "w1", ylab = "w2", col = cols)
points(x = w_opt[1], y = w_opt[2], cex = 1.5, pch = 19, col = "black")
points(x = resultado_LMS$W[1,], y = resultado_LMS$W[2,], cex = 1.5 ,col = "red")

#Número de condicionamento
kappa(R)

#teste
v = c(1,2)
2 * v
e = d[6] - dot(x[5:6], c(1, 1))
e
(2 * 0.2 * e) * x[5:(5+(2-1))]
x[5:6]
w
x[5:6]/(0.01 + dot(x[5:6], x[5:6]))

#####################################################
#Problema 5
#Geração da entrada
x_ruido = rnorm(1011, mean = 0, sd = 1)
ruido_med = rnorm(1011, mean = 0, sd = sqrt(0.001))

#considerando influência aditiva do ruído de medida sobre a entrada
x_entrada = x_ruido + ruido_med

#Obtenção do sinal desejado com base na função de transferência
y = rep(0, 1000)
for (k in 12:length(x_entrada)) {
  y[k-11] = sum(x_entrada[k-11:k])
}

#Representação gráfica dos sinais de entrada e desejado
par(mfrow = c(2, 1))
plot(x = 1:1011, y = x_entrada, type = "l", lwd = 2, col = "Darkred", xlab = "Tempo discreto (k)", ylab = "Observações", main = "Sinal de entrada do filtro")
plot(x = 1:1000, y = y, type = "l", lwd = 2, col = "Darkblue", xlab = "Tempo discreto (k)", ylab = "Observações", main = "Sinal de saída do sistema")

#Executando o algoritmo LMS para diferentes valores de passo de aprendizagem
R = (x_entrada %*% t(x_entrada))[1:11, 1:11]
R = R/(length(x_entrada) - 1)

P = ccf(x_entrada[12: length(x_entrada)], y)
P = P$acf[25:35, 1, 1]

mu_1 = eigen(R)$value[which.max(eigen(R)$value)]/2
mu_2 = eigen(R)$value[which.max(eigen(R)$value)]/10
mu_3 = eigen(R)$value[which.max(eigen(R)$value)]/50

resultado_LMS_mu_1 = LMS(900, mu_1, x_entrada, y, 11)
resultado_LMS_mu_2 = LMS(900, mu_2, x_entrada, y, 11)
resultado_LMS_mu_3 = LMS(900, mu_3, x_entrada, y, 11)

par(mfrow = c(2, 2))
plot(x = 1:900, y = resultado_LMS_mu_1$erro, type = "l", lwd = 2.5, col = "Darkred", xlab = "Iteração", ylab = "MSE", main = "MSE x Iterações para mu_máx/2")
plot(x = 1:900, y = resultado_LMS_mu_2$erro, type = "l", lwd = 2.5, col = "Darkblue", xlab = "Iteração", ylab = "MSE", main = "MSE x Iterações para mu_máx/10")
plot(x = 1:900, y = resultado_LMS_mu_1$erro, type = "l", lwd = 2.5, col = "orange", xlab = "Iteração", ylab = "MSE", main = "MSE x Iterações para mu_máx/50")
resultado_LMS_mu_2$W

#Determinação do Desajuste (Misadjustment)
Misadjustment = function(R, mu){
  return((mu * sum(diag(R)))/(1- (mu * sum(diag(R)))))
}

#Desajuste Empírico
resultado_LMS_mu_1$erro[900]
resultado_LMS_mu_2$erro[900]
resultado_LMS_mu_3$erro[900]

#Desajuste Teórico
M_mu_1 = Misadjustment(R, mu_1)
M_mu_2 = Misadjustment(R, mu_2)
M_mu_3 = Misadjustment(R, mu_3)


