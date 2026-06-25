library(asreml)
library(tidyverse)
library(ggpubr)
library(metan)

# ------------------------------------------------------------------------------
# 1. PREPARAÇÃO DOS DADOS
# ------------------------------------------------------------------------------

url_soja <- "https://raw.githubusercontent.com/mauricioaraujj/Pan_African_Trials_Network/refs/heads/main/data/data.csv"

dat <- read.csv(url_soja, sep = ';') |> 
  filter(COUNTRY == "Zambia") |>
  mutate(
    env  = as.factor(env),
    gen  = as.factor(gen),
    rep  = as.factor(rep)
  ) |>
  filter(!is.na(GY))

cat("Ambientes totais:", nlevels(dat$env), "\n")
cat("Genótipos totais:", nlevels(dat$gen), "\n")

# ------------------------------------------------------------------------------
# 2. ÍNDICE AMBIENTAL (Finlay-Wilkinson)
# ------------------------------------------------------------------------------

index <- data.frame(
  index = tapply(dat$GY, dat$env, mean) - mean(dat$GY)
) |> rownames_to_column("env")

dat_grafico_amb <- merge(dat, index, by = "env")

ggplot() +
  geom_point(data = subset(index, index < 0),
             aes(x = env, y = index, colour = "Desfavorável", fill = "Desfavorável"),
             size = 2.5, shape = 25) +
  geom_point(data = subset(index, index > 0),
             aes(x = env, y = index, colour = "Favorável", fill = "Favorável"),
             size = 2.5, shape = 24) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 1.2) +
  labs(y = "Índice Ambiental", x = "Ambiente", colour = "", fill = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(values = c("firebrick","forestgreen")) +
  scale_fill_manual(values = c("firebrick","forestgreen"))

####### ou
index_trajetoria <- dat %>%
  group_by(env, loc, YEAR, RAINFED) %>% 
  summarise(GY_mean = mean(GY, na.rm = TRUE), .groups = "drop") %>%
  mutate(index = GY_mean - mean(dat$GY, na.rm = TRUE)) %>%
  group_by(loc) %>%
  mutate(anos_totais = n_distinct(YEAR)) %>%
  ungroup()

ggplot(index_trajetoria, aes(x = factor(YEAR), y = index, group = loc)) +
    geom_line(data = subset(index_trajetoria, anos_totais > 1), 
            aes(colour = loc), linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, colour = "grey40") +
    geom_point(data = subset(index_trajetoria, anos_totais > 1),
             aes(shape = as.factor(RAINFED)), colour = "black", size = 4, stroke = 1.2) +
    geom_point(data = subset(index_trajetoria, anos_totais == 1),
             aes(shape = as.factor(RAINFED)), colour = "firebrick", size = 4, stroke = 1.2) +
  
  # Cores para as linhas (Locais)
  scale_colour_viridis_d(option = "turbo", name = "Código do Local") +
    scale_shape_manual(values = c(16, 15, 8), name = "Regime Hídrico") + 
  
  labs(x = "Ano de Cultivo", 
       y = "Índice Ambiental",
       title = "Trajetória dos Locais e Impacto do Regime Hídrico",
       subtitle = "") +
  
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())


# ------------------------------------------------------------------------------
# 3. REGRESSÃO DE EBERHART-RUSSELL
# ------------------------------------------------------------------------------

ERm <- metan::ge_reg(
  .data   = dat,          
  env     = "env",
  gen     = "gen",
  rep     = "rep",
  resp    = "GY",
  verbose = FALSE
)

# Selecionar os genótipos com maior média geral (b0)
top_gen <- ERm$GY$regression |>
  arrange(desc(b0)) |>
  slice(1:30) |> #dá para mudar os 30
  pull(GEN) |>
  as.character()

#mantém os genótipos mesmo que ná esteja em todos os ambientes?

ggplot(data = subset(ERm$GY$data, GEN %in% top_gen),
       aes(x = IndAmb, y = Y)) +
  facet_wrap(. ~ GEN) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#10342d") +
  geom_point(aes(colour = ENV), size = 2) +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Índice Ambiental", y = "Produtividade (GY)", colour = "Ambiente") +
  scale_color_viridis_d(option = "turbo") +
  stat_regline_equation()



# ------------------------------------------------------------------------------
# 4. DECOMPOSIÇÃO AMMI
# ------------------------------------------------------------------------------

media_geral <- mean(dat$GY, na.rm = TRUE)

medias_env <- dat |>
  group_by(env) |>
  summarise(efeito_env = mean(GY, na.rm = TRUE) - media_geral, .groups = "drop")

medias_gen <- dat |>
  group_by(gen) |>
  summarise(efeito_gen = mean(GY, na.rm = TRUE) - media_geral, .groups = "drop")

dat_ammi_completo <- dat |>
  tidyr::complete(env, gen, rep) |>
  left_join(medias_env, by = "env") |>
  left_join(medias_gen, by = "gen") |>
  mutate(
    # Se o GY for NA, preenche com o modelo aditivo; caso contrário, mantém o GY real
    GY = ifelse(is.na(GY), media_geral + efeito_env + efeito_gen, GY)
  ) |>
  select(env, gen, rep, GY) # Manter apenas as colunas necessárias

ammi <- metan::performs_ammi(
  .data = dat_ammi_completo,
  env   = "env",
  gen   = "gen",
  rep   = "rep",
  resp  = "GY"
)

ggarrange(plot_scores(ammi, type = 1), plot_scores(ammi, type = 2))

plot(ammi, type = "AMMI1")
plot(ammi, type = "AMMI2")

# ------------------------------------------------------------------------------
# 5. MODELOS MISTOS
# ------------------------------------------------------------------------------

# No modelo inicial, todos os efeitos aleatórios estao com matrizes de covariancias diagonais homogêneas

#######################################
# Matriz de covariâncias dos Resíduos
#######################################

## mod1: Matriz dos resíduos diagonal homogênea (Todos ambientes tem a mesma variância residual e não há correlação residual nem dentro e nem entre ambientes)
mod1 <- asreml(
  fixed    = GY ~ 1,
  random   = ~ env + env:rep + gen + gen:env,
  residual = ~ idv(units),
  data     = dat,
  maxit = 100)
aic1 <- summary(mod1)$aic

## mod2: Matriz dos resíduos bloco diagonal heterogênea (DIAG) (Cada ambiente possui sua própria variância residual e não há correlação residual nem dentro e nem entre ambientes)
mod2 <- asreml(
  fixed    = GY ~ 1,
  random   = ~ env + env:rep + gen + gen:env,
  residual = ~ dsum(~idv(units) | env),
  data     = dat,
  maxit = 100)
aic2 <- summary(mod2)$aic

## Melhor estrutura é a do mod2


#######################################
# Matriz de covariâncias dos blocos
#######################################

## mod3: Agora cada ambiente tem sua propria variancia do bloco (+32 parametros - não melhorou AIC)
mod3 <- asreml(
  fixed    = GY ~ 1,
  random   = ~ env + at(env):rep + gen + gen:env,
  residual = ~ dsum(~idv(units) | env),
  data     = dat,
  maxit = 100)
aic3 <- summary(mod3)$aic

## Melhor continuar com estrutura do mod2





#######################################
# Matriz de covariâncias GxE
#######################################

## mod4: Variância dos genótipos é diferente para cada ambiente (DIAG)
mod4 <- asreml(
  fixed    = GY ~ 1,
  random   = ~ env + env:rep + diag(env):gen,
  residual = ~ dsum(~idv(units) | env),
  data     = dat,
  maxit = 100)
aic4 <- summary(mod4)$aic


## mod5: Variância dos genótipos é diferente para cada ambiente e correlações entre efeitos dos genótipos nos diferentes ambientes é a mesma
mod5 <- asreml(
  fixed    = GY ~ 1,
  random   = ~ env + env:rep + corh(env):gen,
  residual = ~ dsum(~idv(units) | env),
  data     = dat,
  maxit = 100)
aic5 <- summary(mod5)$aic


## mod6: Variância e correlações entre os genótipos de cada ambiente sao diferentes 
mod6 <- asreml(
  fixed    = GY ~ 1,
  random   = ~ env + env:rep + us(env):gen,
  residual = ~ dsum(~idv(units) | env),
  data     = dat,
  maxit = 300)
aic6 <- summary(mod6)$aic

# nao convergiu

## mod7: Variância e correlações entre os genótipos de cada ambiente sao diferentes 
mod7 <- asreml(
  fixed    = GY ~ 1,
  random   = ~ env + env:rep + corgh(env):gen,
  residual = ~ dsum(~idv(units) | env),
  data     = dat,
  maxit = 300)
aic7 <- summary(mod7)$aic

# nao convergiu



# ## mod4: Variância da interação GxE é diferente para cada ambiente (DIAG)
# mod4 <- asreml(
#   fixed    = GY ~ 1,
#   random   = ~ env + env:rep + gen + gen:diag(env),
#   residual = ~ dsum(~idv(units) | env),
#   data     = dat,
#   maxit = 100)
# aic4 <- summary(mod4)$aic



# ------------------------------------------------------------------------------
# 6. COMPARAÇÃO DE MODELOS PELO AIC
# ------------------------------------------------------------------------------

modelos_aic <- data.frame(
  modelo = factor(
    c("CS(G)-ID(R)", "CS(G)-DIAG(R)", "CS(G)-DIAG(B)-DIAG(R)", "DIAG(G)-DIAG(R)",
      "HCS(G)-DIAG(R)"),
    levels = c("CS(G)-ID(R)", "CS(G)-DIAG(R)", "CS(G)-DIAG(B)-DIAG(R)", "DIAG(G)-DIAG(R)",
               "HCS(G)-DIAG(R)") 
  ), AIC = c(aic1, aic2, aic3, aic4, aic5))
  
  ggplot(modelos_aic, aes(y = modelo, x = AIC)) +
    geom_point(size = 3, shape = 15, colour = "#10342d") +
    geom_segment(aes(x = min(AIC) - 50, xend = AIC, y = modelo, yend = modelo),
                 linewidth = 1.5, colour = "#10342d") +
    geom_label(aes(label = round(AIC, 1)), size = 4,
               position = position_nudge(y = 0.3, x = -30)) +
    theme_minimal(base_size = 14) +
    labs(y = "Modelo", x = "AIC",
         title = "Comparação de estruturas de covariância G×A",
         subtitle = "Variável resposta: GY | Dados Zâmbia (Soja)")

  
  #model 5 teve menor IAC
  
  # ------------------------------------------------------------------------------
  # 7. PREDIÇÕES DO MELHOR MODELO 
  # ------------------------------------------------------------------------------
  
  melhor_mod <- mod5 # o model 5 teve o melnor AIC
  
  pred <- predict.asreml(
    object   = melhor_mod,
    classify = "gen:env"
  )$pvals
  
  pred |>
    ggplot(aes(x = reorder(env, -predicted.value),
               y = reorder(gen, predicted.value),
               fill = predicted.value)) +
    geom_tile() +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + 
    scale_fill_gradient(low = "#10342d", high = "#b2bc63") +
    labs(x = "Ambiente", y = "Genótipo", fill = "GY",
         title = "Valores preditos — Genótipo × Ambiente")
  
  
  
  # ------------------------------------------------------------------------------
  # 8. ÍNDICES DE ESTABILIDADE 
  # ------------------------------------------------------------------------------
  
  ge_table <- pred[, c("gen", "env", "predicted.value")] |>
    pivot_wider(names_from = "env", values_from = "predicted.value") |>
    column_to_rownames("gen")
  
  # -- Índice Pi (Lin & Binns) --
  Pi <- apply(ge_table, 1, function(x) {
    sum((x - apply(ge_table, 2, max))^2) / (2 * ncol(ge_table))
  })
  
  df_Pi <- data.frame(gen = names(Pi), Pi = Pi) |>
    arrange(Pi) |>
    slice(1:15) #top 15
  
  ggplot(data = df_Pi, aes(x = reorder(gen, Pi), y = Pi)) +
    geom_col(colour = "black", aes(fill = Pi)) +
    scale_fill_gradient(high = "#10342d", low = "#b2bc63") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none") +
    labs(y = expression(P[i]), x = "Genótipo",
         title = "Top 15 - Índice Pi (Lin & Binns)")
  
  # -- Safety-First Index (Kataoka) --
  alpha <- 0.10
  SF <- apply(ge_table, 1, function(x) mean(x) - qnorm(1 - alpha) * sd(x))
  
  df_SF <- data.frame(gen = names(SF), SF = SF) |>
    arrange(desc(SF)) |>
    slice(1:15)
  
  ggplot(data = df_SF, aes(x = reorder(gen, -SF), y = SF)) +
    geom_col(colour = "black", aes(fill = SF)) +
    scale_fill_gradient(low = "#10342d", high = "#b2bc63") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none") +
    labs(y = "Safety-First Index", x = "Genótipo",
         title = "Top 15 - Safety-First (Kataoka)")
  
  
  
  
  ################
  # Conectividade
  ###############
  
  # ------------------------------------------------------------------------------
  # 9. CONECTIVIDADE DOS GENÓTIPOS ENTRE AMBIENTES
  # ------------------------------------------------------------------------------
  
  table_connect <- table(as.character(df_filter$env), as.character(df_filter$gen))
  
  # Convert to presence/absence
  table_connect[table_connect > 0] <- 1
  table_connect <- as.data.frame(table_connect)
  
  # Percentage of experiments in which each treatment is present
  conect_perc <- tapply(table_connect$Freq, table_connect$Var2, sum) /                      length(unique(df_filter$env)) * 100
  
  conect_perc <- data.frame(
    Genótipo = names(conect_perc),
    Conectividade = conect_perc
  )
  
  
  conect_perc <- conect_perc[
    order(conect_perc$Conectividade, decreasing = TRUE),
  ]
  
  kable(
    conect_perc,
    caption = "Connectivity of genotypes across experiments",
    align = "c"
  )
  
  table_connect
  
  #Presence 
  pres <- mean(table_connect$Freq) * 100
  
  #Absence
  abs <- 100 - pres
  
  
  table_connect$Var2 <- factor(as.character(table_connect$Var2), levels = conect_perc$Genótipo)
  ordem_env <- names(sort(tapply(table_connect$Freq, table_connect$Var1, sum), decreasing = TRUE))
  table_connect$Var1 <- factor(as.character(table_connect$Var1), levels = ordem_env)
  
  ggplot(table_connect, aes(x = Var2, y = Var1, fill = factor(Freq))) +
    geom_tile() +
    scale_fill_manual(
      values = c("white", "#A7C7E7"),
      labels = c("Absence", "Presence")
    ) +
    labs(
      title = "Presence–absence of treatments across experiments",
      x = "Genotypes (Sorted by Connectivity)",
      y = "Environments",
      fill = "Presence"
    ) +
    theme_minimal() + 
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank() 
    )
  
  # Matrix
  M <- xtabs(Freq ~ Var1 + Var2, data = table_connect)
  
  # Shared treatments
  shared <- M %*% t(M)
  
  # Number of treatments per experiment
  n_trat <- rowSums(M)
  
  # Jaccard
  jaccard <- shared
  
  for(i in rownames(shared)) {
    for(j in colnames(shared)) {
      jaccard[i, j] <- shared[i, j] / (n_trat[i] + n_trat[j] - shared[i, j])
    }
  }
  
  diag(jaccard) <- 1
  
  pheatmap(
    jaccard,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "average",
    color = rev(hcl.colors(100, "Blues")),
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 90,
    main = "Experiment connectivity based on shared genotypes"
  )
  
  ## Heatmap com os genótipos
  
  # 1. Matriz de contingência original (Linhas = Ambientes, Colunas = Genótipos)
  M <- xtabs(Freq ~ Var1 + Var2, data = table_connect)
  
  # 2. Shared environments per genotype (INVERTIDO: t(M) %*% M gera Genótipo x Genótipo)
  shared_gen <- t(M) %*% M
  
  # 3. Número de ambientes em que cada genótipo está presente (Agora sim usamos colSums)
  n_env_per_gen <- colSums(M)
  
  # 4. Inicializa a matriz de Jaccard para os Genótipos
  jaccard_gen <- shared_gen
  
  # 5. Loop do Jaccard adaptado para os genótipos (com proteção contra divisões por zero)
  for(i in rownames(shared_gen)) {
    for(j in colnames(shared_gen)) {
      denominador <- n_env_per_gen[i] + n_env_per_gen[j] - shared_gen[i, j]
      
      if (denominador == 0) {
        jaccard_gen[i, j] <- 0
      } else {
        jaccard_gen[i, j] <- shared_gen[i, j] / denominador
      }
    }
  }
  
  # Garante a diagonal unitária (similaridade de um genótipo com ele mesmo é 1)
  diag(jaccard_gen) <- 1
  
  # 6. Plot do Heatmap dos Genótipos
  library(pheatmap)
  pheatmap(
    jaccard_gen,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "average",
    color = rev(hcl.colors(100, "Blues")),
    fontsize_row = 6,  
    fontsize_col = 6,
    angle_col = 90,
    main = "Genotype connectivity based on shared environments"
  )
