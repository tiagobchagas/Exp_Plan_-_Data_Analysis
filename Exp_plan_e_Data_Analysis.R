# Planejamento Experimental e Análises de dados agronômicos
# Autor: Dsc. José Tiago Barroso Chagas



# Instalação de pacotes ---------------------------------------------------

# Os seguinte pacotes são requeridos para execução deste script
# Favor instalá-los antes de prosseguir para as demais etapas
pacotes=c("psych","metan","ggcorrplot", "lme4", "FieldSimR", "ggrepel",
          "knitr", "kableExtra", "ggpubr")

install.packages(pacotes)

sapply(pacotes, library, character.only=T, logical.return=T)

# Simulando o experimento de Milho ----------------------------------------

library(FieldSimR)

ntraits <- 2 # Number of traits
nenvs <- 3 # Number of environments
nblocks <- c(2, 2, 3) # Number of blocks per environment
block_dir <- "col" # Arrangement of blocks ("side-by-side")
ncols <- c(10, 10, 15) # Number of columns per environment
nrows <- 20 # Number of rows per environment
plot_length <- 8 # Plot length; here in meters (column direction)
plot_width <- 2 # Plot width; here in meters (row direction)

#H2 = 0.3 for grain yield and H2 = 0.5 for plant height 

H2 <- c(0.3, 0.3, 0.3, 0.5, 0.5, 0.5) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)
var <- c(0.086, 0.12, 0.06, 15.1, 8.5, 11.7) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)
#set.seed(2024)

# Calculation of error variances based on the genetic variance and target heritability vectors.
calc_varR <- function(var, H2) {
  varR <- (var / H2) - var
  return(varR)
}

varR <- calc_varR(var, H2)
cat("Variâncias residuais Trait per environment -
    c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3")
round(varR, 2) # Vector of error variances: c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)


spatial_model <- "Bivariate" # Spatial error model.
prop_spatial <- 0.4 # Proportion of spatial trend.

ScorR <- rand_cor_mat(ntraits, min.cor = 0, max.cor = 0.5, pos.def = TRUE)
#Correlação dos erros devido tendências espaciais
round(ScorR, 2)
# extraneous variation
ext_ord <- "zig-zag"
ext_dir <- "row"
prop_ext <- 0.2

EcorR <- rand_cor_mat(ntraits, min.cor = 0, max.cor = 0.5, pos.def = TRUE)
# Correlação dos erros devido causas não conhecidas   (Caminhamento,máquinas e implementos)
round(EcorR, 2)
gv_df <- gv_df_unstr

error_ls <- field_trial_error(
  ntraits = ntraits,
  nenvs = nenvs,
  nblocks = nblocks,
  block.dir = block_dir,
  ncols = ncols,
  nrows = nrows,
  plot.length = plot_length,
  plot.width = plot_width,
  varR = varR,
  ScorR = ScorR,
  EcorR = EcorR,
  RcorR = NULL,
  spatial.model = spatial_model,
  prop.spatial = prop_spatial,
  ext.ord = ext_ord,
  ext.dir = ext_dir,
  prop.ext = prop_ext,
  return.effects = TRUE
)
pheno_df <- make_phenotypes(
  gv_df,
  error_ls$error.df,
  randomise = TRUE
)


# coloca banco de dados real



# Blups --------------------------------------------------------
library(metan)
inspect(pheno_df,plot=TRUE,threshold = 100)

colnames(pheno_df)=c("env","block","col","row","id","rep","GY","PH")
head(pheno_df)
modelb= gamem_met(pheno_df,
                  env = env,
                  gen = id,
                  rep = c("rep"),
                  resp = c("GY","PH"), random = "gen",
                  verbose = TRUE)

modelb$GY$MeansGxE$Trait="GY"
modelb$PH$MeansGxE$Trait="PH"
a=modelb$GY$MeansGxE
b=modelb$PH$MeansGxE
blups=rbind(a,b)

get_model_data(modelb, what = "vcomp")
get_model_data(modelb, "lrt")

plot(modelb,var = 1)
plot(modelb,var = 2)

head(blups)
bmean = blups |>
  pivot_wider(names_from = c("Trait","ENV"),values_from = c("Y"),id_cols = "GEN") |> 
column_to_rownames("GEN") 
head(bmean)
# Análise Fatorial --------------------------------------------------------



library(psych) 

factor_analysis <- principal(bmean, nfactors = length(colnames(bmean)),
                             rotate = "none",scores = TRUE) 

#Autovales
eigenvalues <- round(factor_analysis$values, 5)
eigenvalues


# Variância compartilhada
variancia_compartilhada <- as.data.frame(factor_analysis$Vaccounted) |>  
  slice(1:3)
rownames(variancia_compartilhada) <- c("Autovalores",
                                       "Prop. da Variância",
                                       "Prop. da Variância Acumulada")
# Variância compartilhada pelas variáveis originais para a formação de cada fator
round(variancia_compartilhada, 3) %>%
  kable(caption = "Autovalores") %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = FALSE, 
                font_size = 20)

# Escores fatoriais
round(factor_analysis$weights,4) |> kable(caption = "Escores") |>
  kable_styling(bootstrap_options = "striped", 
                full_width = FALSE, 
                font_size = 12)

# Cálculo dos fatores propriamente ditos
fatores <- as.data.frame(factor_analysis$scores)
round(factor_analysis$scores,4) |> kable(caption = "Fatores") |>
  kable_styling(bootstrap_options = "striped", 
                full_width = FALSE, 
                font_size = 12)
# Cargas fatoriais
round(unclass(factor_analysis$loadings),4) |> kable(caption = "Cargas") |>
  kable_styling(bootstrap_options = "striped", 
                full_width = FALSE, 
                font_size = 12) 
# Cálculo das cargas fatoriais
cargas_fatoriais <- as.data.frame(unclass(factor_analysis$loadings))



# Coeficientes de correlação de Pearson para cada par de fatores (ortogonais)
rho <- cor(as.matrix(fatores))
round(rho, 4)

# cargas_fatoriais dos dois primeiros fatores
P1=cargas_fatoriais[, c(1,2)] |> data.frame() %>%
  rownames_to_column("variáveis") %>%
  ggplot(aes(x = PC1, y = PC2, label = variáveis)) +
  geom_point(color = "darkorchid",
             size = 3) +
  geom_text_repel() +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "orange") +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "orange") +
  expand_limits(x= c(-1.25, 0.25), y=c(-0.25, 1)) +
  theme_bw()


plot(P1)





# Índice de seleção -------------------------------------------------------




library(metan)
head(bmean)
# selecionar os 20% melhores
bmean |>  fai_blup(SI = 20,
                DI = c("max, max, max, min, min, min"),
                UI = c("min, min, min, max, max, max"),
                mineval = 0.7) |> plot()
#classificar todos os genótipos
FAI= bmean |>  fai_blup(SI = 100,
                   DI = c("max, max, max, min, min, min"),
                   UI = c("min, min, min, max, max, max"),
                   mineval = 0.7) 

library(corrplot)
FAI$cormat |> corrplot()


FAI$eigen$eigen.values 
eigenvalues 



# Rotação das cargas fatoriais --------------------------------------------

factor_analysis <- principal(bmean, nfactors = length(colnames(bmean)),
                             rotate = "varimax",scores = TRUE) 

# Cálculo das cargas fatoriais
cargas_fatoriais <- as.data.frame(unclass(factor_analysis$loadings))
cargas_fatoriais

# cargas_fatoriais dos dois primeiros fatores
P2=cargas_fatoriais[, c(1,2)] |> data.frame() %>%
  rownames_to_column("variáveis") %>%
  ggplot(aes(x = RC5 , y = RC3, label = variáveis)) +
  geom_point(color = "darkorchid",
             size = 3) +
  geom_text_repel() +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "orange") +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "orange") +
  expand_limits(x= c(-1.25, 0.25), y=c(-0.25, 1)) +
  theme_bw()
plot(P2)
library(ggpubr)
ggarrange(P1,P2)
