#Aula teórico-prática de matrizes

#contruindo uma matriz

dados <- c(10,5,5,1,20,5,15,2,30,5,25,3,40,5,35,4)

A <- matrix(data = dados, nrow = 4, ncol = 4) #assim os dados são construídos pelas colunas
A

A1 <- matrix(data = dados, nrow = 4, ncol = 4, byrow = TRUE) #assim os dados são lidos pelas linhas
A1

B <- matrix(c(1,5,9,
              2,8,3,
              10,4,7,
              6,2,5), nrow = 3, ncol = 4) 

C <- matrix(c(1,3,
              5,2,
              8,10,
              4,6), nrow = 2, ncol = 4)

#Traço da matriz

trace <- sum(diag(A))
trace_B <- sum(diag(B)) #o traço de uma matriz retangular não tem significado matricial. O R calcula, mas não pega os valores
                        # corretos. Isso mostra que nem porque o R calcula, será correto.

#Posto
qr(A)
qr(A)$rank

#Saber o posto é saber quais dados são relevantes para a análise também.

qr(B)$rank
qr(C)

A*B #Dimensão não compatível 
B*A #dimensão não compatível por isso que ele não fez a multiplicação

A*B %*% t(C)
C %*% t(B)
B %*% t(B)

#SOma de vetores
v <- c(2,4,6,8)
sum(v)

u <- rep(1, length(v))
t(u) %*% v


#Combinação linear vetores
u <- c(1,2,0,-1)
t(u) %*% v


