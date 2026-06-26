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
    rep  = as.factor(rep)) |>
  filter(!is.na(GY))

cat("Ambientes totais:", nlevels(dat$env), "\n")
cat("Genótipos totais:", nlevels(dat$gen), "\n")

# ------------------------------------------------------------------------------
# 2. ÍNDICE AMBIENTAL (Finlay-Wilkinson)
# ------------------------------------------------------------------------------

```{r}
#| label: fig-trajectory
#| fig-cap: "Environmental Index trajectory of the evaluated locations across years"
#| warning: false
#| message: false

library(tidyverse)
library(ggplot2)
library(ggrepel)

trajectory_index <- dat %>%
  group_by(env, loc, YEAR, RAINFED) %>% 
  summarise(GY_mean = mean(GY, na.rm = TRUE), .groups = "drop") %>%
  mutate(index = GY_mean - mean(dat$GY, na.rm = TRUE)) %>%
  group_by(loc) %>%
  mutate(total_years = n_distinct(YEAR)) %>%
  ungroup()

ggplot(trajectory_index, aes(x = factor(YEAR), y = index, group = loc)) +
  geom_line(data = subset(trajectory_index, total_years > 1), 
            aes(colour = loc), linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, colour = "grey40") +
  
  geom_point(data = subset(trajectory_index, total_years > 1),
             aes(shape = as.factor(RAINFED)), colour = "black", size = 4, stroke = 1.2) +
  geom_point(data = subset(trajectory_index, total_years == 1),
             aes(shape = as.factor(RAINFED)), colour = "firebrick", size = 4, stroke = 1.2) +
  
  ggrepel::geom_text_repel(aes(label = env), colour = "grey50", size = 3, 
                           nudge_x = 0.2,       
                           direction = "y",    
                           segment.color = "grey80") + 
  
  scale_colour_viridis_d(option = "turbo", name = "Location Code") +
  scale_shape_manual(values = c(16, 15, 8), name = "Water Regime") + 
  
  labs(x = "Year", 
       y = "Environmental Index",
       title = "Location Trajectory and Water Regime Impact",
       subtitle = "") +
  
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

```
# ------------------------------------------------------------------------------
# 3. REGRESSÃO DE EBERHART-RUSSELL
# ------------------------------------------------------------------------------

```{r}
#| label: fig-eberhart-russell
#| fig-cap: "Eberhart and Russell (1966) regression analysis"
#| 
# 1. Run the Eberhart-Russell model
ERm <- metan::ge_reg(
  .data   = dat,          
  env     = "env",
  gen     = "gen",
  rep     = "rep",
  resp    = "GY",
  verbose = FALSE
)

# 2. Select the top 30 genotypes with the highest overall mean (b0)
top_gen <- ERm$GY$regression |>
  arrange(desc(b0)) |>
  slice(1:30) |> 
  pull(GEN) |>
  as.character()

# 3. Force genotype G318 into the plot (if not already in the top 30)
top_gen <- unique(c(top_gen, "G318"))

# 4. Plot the regression lines
ggplot(data = subset(ERm$GY$data, GEN %in% top_gen),
       aes(x = IndAmb, y = Y)) +
  facet_wrap(. ~ GEN) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#10342d") +
  geom_point(aes(colour = ENV), size = 2) +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Environmental Index", 
       y = "Grain Yield (GY)", 
       colour = "Environment") +
  scale_color_viridis_d(option = "turbo") +
  stat_regline_equation()
```



# ------------------------------------------------------------------------------
# 4. DECOMPOSIÇÃO AMMI
# ------------------------------------------------------------------------------

library(tidyverse)
library(metan)
library(patchwork)

ammi_bruto <- metan::performs_ammi(
  .data = dat,
  env   = "env",
  gen   = "gen",
  rep   = "rep",
  resp  = "GY"
)

p1_bruto <- metan::plot_scores(ammi_bruto, 
                               type = 1, 
                               size.text.gen = 0,   
                               col.gen = "grey60",  
                               size.text.env = 3.5, 
                               col.env = "#10342d",
                               title = FALSE) +     
  labs(title = "A) AMMI 1: Prod. vs. PC1")

p2_bruto <- metan::plot_scores(ammi_bruto, 
                               type = 2, 
                               size.text.gen = 0, 
                               col.gen = "grey60",
                               size.text.env = 3.5,
                               col.env = "firebrick",
                               title = FALSE) +       
  labs(title = "B) AMMI 2: PC1 vs PC2")

# 4. Juntar os gráficos lado a lado
painel_ammi_bruto <- p1_bruto + p2_bruto

# Mostrar o gráfico final no ecrã
print(painel_ammi_bruto)

# ------------------------------------------------------------------------------
# 5. MODELOS MISTOS
# ------------------------------------------------------------------------------

# No modelo inicial, todos os efeitos aleatórios estao com matrizes de covariancias diagonais homogêneas

#######################################
# Matriz de covariâncias dos Resíduos
#######################################





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
  
  melhor_mod <- mod5 # o model 5 teve o menor AIC
  
  pred <- predict.asreml(object = melhor_mod,classify = "gen:env")$pvals
  
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
  
  
  # G342, G307, G310
  
  varcomp <- summary(mod5)$varcomp
  
  pred <- predict.asreml(object = melhor_mod,classify = "gen:env")$pvals
  pred$reliability <- 1 - (pred$std.error^2/varcomp[1,1])
  
  
  ################
  # Conectividade
  ###############
  
  # ------------------------------------------------------------------------------
  # 9. CONECTIVIDADE DOS GENÓTIPOS ENTRE AMBIENTES
  # ------------------------------------------------------------------------------

table_connect <- table(as.character(df_filter$env), as.character(df_filter$gen))
table_connect[table_connect > 0] <- 1
table_connect <- as.data.frame(table_connect)

df_water <- df_filter %>% 
  distinct(env, RAINFED) %>% 
  rename(Var1 = env) 

table_connect <- merge(table_connect, df_water, by = "Var1", all.x = TRUE)
table_connect$Status <- ifelse(table_connect$Freq == 0, "Absence", as.character(table_connect$RAINFED))
conect_perc <- tapply(table_connect$Freq, table_connect$Var2, sum) / length(unique(df_filter$env)) * 100

conect_perc <- data.frame(
  Genótipo = names(conect_perc),
  Conectividade = conect_perc
)

conect_perc <- conect_perc[order(conect_perc$Conectividade, decreasing = TRUE), ]

kable(
  conect_perc,
  caption = "Connectivity of genotypes across experiments",
  align = "c"
)

table_connect$Var2 <- factor(as.character(table_connect$Var2), levels = conect_perc$Genótipo)
ordem_env <- names(sort(tapply(table_connect$Freq, table_connect$Var1, sum), decreasing = TRUE))
table_connect$Var1 <- factor(as.character(table_connect$Var1), levels = ordem_env)

ggplot(table_connect, aes(x = Var2, y = Var1, fill = Status)) +
  geom_tile() +
  scale_fill_manual(
    values = c(
      "Absence"       = "white", 
      "Irrigation"    = "#1b4332",  
      "Supplementary" = "#40916c",  
      "Rainfed"       = "#95d5b2"   
    )
  ) +
  labs(
    title = "Presence–absence of treatments across experiments by water regime",
    x = "Genotypes",
    y = "Environments",
    fill = "Water Regime"
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank() 
  )

genotipos_alvo <- c("G342", "G138") 
niveis_genotipos <- levels(table_connect$Var2)
rotulos_customizados <- ifelse(niveis_genotipos %in% genotipos_alvo, niveis_genotipos, "")
cores_eixo_x <- ifelse(niveis_genotipos %in% genotipos_alvo, "red", "black")
ggplot(table_connect, aes(x = Var2, y = Var1, fill = Status)) +
  geom_tile() +
  scale_fill_manual(
    values = c(
      "Absence"       = "white", 
      "Irrigation"    = "#1b4332",  
      "Supplementary" = "#40916c",  
      "Rainfed"       = "#95d5b2"   
    )
  ) +
  scale_x_discrete(labels = rotulos_customizados) + 
  labs(
    title = "Presence–absence of treatments across experiments by water regime",
    x = "Genotypes",
    y = "Environments",
    fill = "Water Regime"
  ) +
  theme_minimal() + 
  theme(
  
    axis.text.x = element_text(color = cores_eixo_x, angle = 90, vjust = 0.5, hjust = 1, size = 8, face = "bold"),
    axis.ticks.x = element_line(color = ifelse(rotulos_customizados == "", NA, "gray")), 
    panel.grid = element_blank() 
  )
  
  ###########################
  # Genotypes
  ##########################
  g342 <- dat |> filter(gen == "G342") |> select (env) |> unique()
  g307 <- dat |> filter(gen == "G307") |> select (env) |> unique()
  g310 <- dat |> filter(gen == "G310") |> select (env) |> unique()
  
  
  # Melhores sequeiro
  top5 <- pred |>
    filter(env %in% c("E0103","E0104","E0105","E0106","E0107","E0108",
                      "E0169","E0170","E0171","E0172","E0173","E0174",
                      "E0175","E0176")) |>
    group_by(env) |>
    arrange(desc(predicted.value), .by_group = TRUE) |>
    slice_head(n = 5)
  
  # G342 - 13/14
  # G307 - 7/14
  # G310 - 4/14
  
  #############################
  # Análise modelo 5
  ############################
  
  mod5 <- asreml(
    fixed    = GY ~ 1,
    random   = ~ env + env:rep + corh(env):gen,
    residual = ~ dsum(~idv(units) | env),
    data     = dat,
    maxit = 100)
  aic5 <- summary(mod5)$aic
  
  plot(mod5)
  
  pred <- predict.asreml(object = mod5,classify = "gen:env",sed=TRUE)
  varcomp <- summary(mod5)$varcomp
  
  h2_mod5 <- 1-(mean(pred$sed[upper.tri(pred$sed)]^2)/(2*mean(varcomp[4:36,1])))
  
  
  # ------------------------------------------------------------------------------
  # 10. RELIABILITY DE CADA GENÓTIPO EM CADA ENV
  # ------------------------------------------------------------------------------
  library(dplyr)
  library(asreml)
  
  data <- read.csv(url_soja, sep = ';') |> 
    filter(COUNTRY == "Zambia") |>
    mutate(
      env  = as.factor(env),
      gen  = as.factor(gen),
      rep  = as.factor(rep)) |>
    filter(!is.na(GY))
  
  
  ## Removing outliers
  df <- data|> 
    mutate (GY = case_when(
      GY == 10650.85 ~ NA,
      TRUE ~ GY
    )) |> 
    mutate(PH_R8 = case_when(
      PH_R8 == 0 ~ NA,
      PH_R8 == 0.0 ~ NA,
      PH_R8 == 308.33 ~ NA,
      TRUE ~ PH_R8
    ))
  
  ## Converting columns
  df$YEAR <- as.factor(df$YEAR)
  df$loc <- as.factor(df$loc)
  df$env <- as.factor(df$env)
  df$check <- as.factor(df$check)
  df$rep <- as.factor(df$rep)
  df$gen <- as.factor(df$gen)
  df$RAINFED <- as.factor(df$RAINFED)
  
  df_filter <- df |> filter(!(env %in% c("E0108", "E0175", "E0285")))
  
  # Model 5
  mod5 <- asreml(
    fixed    = GY ~ 1,
    random   = ~ env + env:rep + corh(env):gen,
    residual = ~ dsum(~idv(units) | env),
    data     = df_filter,
    maxit = 100)
  aic5 <- summary(mod5)$aic ; aic5
  
  print(summary(mod5)$varcomp)
  varcomp <- summary(mod5)$varcomp
  
  preds <- predict(mod5, classify = "env:gen", sed = TRUE, pworkspace = "8000mb")
  
  tabela_blups <- preds$pvals
  
  # Calcular a PEV (Prediction Error Variance)
  tabela_blups$PEV <- tabela_blups$std.error^2
  
  print(summary(mod5)$varcomp)
  varcomp <- summary(mod5)$varcomp
  
  sigma2_G_valores <- c(
    "E0103" = varcomp[4,1],  
    "E0104" = varcomp[5,1],  
    "E0105" = varcomp[6,1],   
    "E0106" = varcomp[7,1],
    "E0107" = varcomp[8,1],
    "E0169" = varcomp[9,1],
    "E0170" = varcomp[10,1],
    "E0171" = varcomp[11,1],
    "E0172" = varcomp[12,1],
    "E0173" = varcomp[13,1],
    "E0174" = varcomp[14,1],
    "E0176" = varcomp[15,1],
    "E0224" = varcomp[16,1],
    "E0225" = varcomp[17,1],
    "E0226" = varcomp[18,1],
    "E0227" = varcomp[19,1],
    "E0228" = varcomp[20,1],
    "E0229" = varcomp[21,1],
    "E0230" = varcomp[22,1],
    "E0252" = varcomp[23,1],
    "E0253" = varcomp[24,1],
    "E0254" = varcomp[25,1],
    "E0255" = varcomp[26,1],
    "E0278" = varcomp[27,1],
    "E0279" = varcomp[28,1],
    "E0280" = varcomp[29,1],
    "E0281" = varcomp[30,1],
    "E0282" = varcomp[31,1],
    "E0283" = varcomp[32,1],
    "E0284" = varcomp[33,1]
  )
  
  tabela_blups$sigma2_G <- sigma2_G_valores[tabela_blups$env]
  
  # Reliability = 1 - (PEV / Vgen)
  tabela_blups$Reliability <- 1 - (tabela_blups$PEV / tabela_blups$sigma2_G)
  
  # RESULTADO FINAL
  resultado_final <- tabela_blups[, c("env", "gen", "predicted.value", "std.error", "Reliability")]
  head(resultado_final)
  
  mod5_reliability <- readRDS("reability2.rds")
  
  
  rel_gen <- mod5_reliability |>
    group_by(gen) |>
    summarise(
      mean_reliability = mean(Reliability, na.rm = TRUE),
      n_env = n(),
      .groups = "drop"
    )
  
  rel_env <- mod5_reliability |>
    group_by(env) |>
    summarise(
      mean_reliability = mean(Reliability, na.rm = TRUE),
      n_gen = n(),
      .groups = "drop"
    )
  
