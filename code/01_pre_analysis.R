# negative time points: initation + burn-in samples
library(sjSDM)
library(tidyverse)
library(DALEX)

load("sim0.RData")
data = sim.data[[1]]
coords = sim.data[[2]]
data = cbind(coords, data)
A = sim.data[[3]]
Niche = sim.data[[4]]
data_wo_burnin= data[data$time >= 0, ]


# Spatial snapshot --------------------------------------------------------
# Several patches at on one time-point

# Prepare data
spatial_snapshot = data_wo_burnin %>% filter(time == 1200)
X = spatial_snapshot %>% select(env)
Y = spatial_snapshot[, 6:ncol(spatial_snapshot)] %>% as.matrix()
spatial_coords = spatial_snapshot %>% select(x, y) %>% scale
surviving_species = colSums(Y) > 0

# Train model
spatial_poisson = sjSDM(Y[, surviving_species], 
                        linear(X, ~env + I(env^2)), 
                        spatial = linear(spatial_coords, formula = ~ x + y + x:y + I(x^2) + I(y^2)),
                        family = poisson(),
                        device = 0)

PA = apply(Y[, surviving_species], 1:2, function(p) ifelse(p > 0,1,0))
spatial_PA= sjSDM(PA, 
                  linear(X, ~env + I(env^2)), 
                  spatial = linear(spatial_coords, formula = ~ x + y + x:y + I(x^2) + I(y^2)),
                  family = binomial(),
                  device = 0)


# Prepare outputs - Abundance
beta_env = t(coef(spatial_poisson)$env[[1]])
rownames(beta_env) = colnames(spatial_poisson$data$X)
beta_space = t(coef(spatial_poisson)$spatial[[1]])
rownames(beta_space) = colnames(spatial_poisson$spatial$X)

results_spatial_abundance = list(
  model = spatial_poisson,
  data = list(Y = spatial_poisson$data$Y,
              env = spatial_poisson$data$X,
              space = spatial_poisson$spatial$X),
  beta_env = beta_env,
  beta_space = beta_space,
  cooccurrence = cov2cor(getCov(spatial_poisson)),
  link = function(P) exp(P)
)

# Prepare outputs - PA
beta_env = t(coef(spatial_PA)$env[[1]])
rownames(beta_env) = colnames(spatial_PA$data$X)
beta_space = t(coef(spatial_PA)$spatial[[1]])
rownames(beta_space) = colnames(spatial_PA$spatial$X)

results_spatial_PA = list(
  model = spatial_PA,
  data = list(Y = spatial_PA$data$Y,
              env = spatial_PA$data$X,
              space = spatial_PA$spatial$X),
  beta_env = beta_env,
  beta_space = beta_space,
  cooccurrence = cov2cor(getCov(spatial_PA)),
  link = function(P) 1/(1+exp(-P))
)

results_spatial = list(
  results_abundance = results_spatial_abundance,
  results_PA = results_spatial_PA,
  A = A[surviving_species, surviving_species],
  Niche = Niche[surviving_species, ]
)


# Temporal snapshot -------------------------------------------------------
# One patch over several time-points
time_indices = 1001:1201#seq(201, by = 5, length.out = 201)
temporal_snapshot = (data_wo_burnin %>% filter(patch == 42))[time_indices,]
X = temporal_snapshot %>% select(env)
Y = temporal_snapshot[, 6:ncol(temporal_snapshot)] %>% as.matrix()
surviving_species = colSums(Y) > 0
Y = Y[,surviving_species]
X = cbind(X[2:201,], scale(Y[1:200,]+0.0001))
X = data.frame(X)
colnames(X) = c("env", paste0("sp_", 1:sum(surviving_species)))
form = paste0("~env + I(env^2) + ", paste0( paste0( "sp_", 1:sum(surviving_species) ), collapse = " + ") )

# Train model - Poisson
temporal_poisson = sjSDM(Y[-1,], 
                         linear(X, form), 
                         family = poisson(),
                         device = 0)
# Train model - Abundance
X = temporal_snapshot %>% select(env)
Y = temporal_snapshot[, 6:ncol(temporal_snapshot)] %>% as.matrix()
surviving_species = colSums(Y) > 0
Y = Y[,surviving_species]
PA = apply(Y, 1:2, function(P) ifelse(P>0, 1, 0))
X = cbind(X[2:201,], scale( PA[1:200,] ))
X = data.frame(X)
colnames(X) = c("env", paste0("sp_", 1:sum(surviving_species)))
form = paste0("~env + I(env^2) + ", paste0( paste0( "sp_", 1:sum(surviving_species) ), collapse = " + ") )

temporal_PA = sjSDM(PA[-1,], 
                    linear(X, form), 
                    family = binomial(),
                    device = 0)


# Prepare outputs - Abundance
beta_env = t(coef(temporal_poisson)[[1]])[1:3,]
rownames(beta_env) = colnames(temporal_poisson$data$X)[1:3]
A_species = t(coef(temporal_poisson)[[1]])[-(1:3),]
rownames(A_species) = colnames(temporal_poisson$data$X)[-(1:3)]

results_temporal_abundance = list(
  model = temporal_poisson,
  data = list(Y = temporal_poisson$data$Y,
              env = temporal_poisson$data$X[,1:3],
              biotic = temporal_poisson$data$X[,-(1:3)]),
  beta_env = beta_env,
  A_species = A_species,
  cooccurrence = cov2cor(getCov(temporal_poisson)),
  link = function(P) exp(P)
)

# Prepare outputs - PA
beta_env = t(coef(temporal_PA)[[1]])[1:3,]
rownames(beta_env) = colnames(temporal_PA$data$X)[1:3]
A_species = t(coef(temporal_PA)[[1]])[-(1:3),]
rownames(A_species) = colnames(temporal_PA$data$X)[-(1:3)]

results_temporal_PA = list(
  model = temporal_PA,
  data = list(Y = temporal_PA$data$Y,
              env = temporal_PA$data$X[,1:3],
              biotic = temporal_PA$data$X[,-(1:3)]),
  beta_env = beta_env,
  A_species = A_species,
  cooccurrence = cov2cor(getCov(temporal_PA)),
  link = function(P) 1/(1+exp(-P))
)

results_temporal = list(
  results_abundance = results_temporal_abundance,
  results_PA = results_temporal_PA,
  A = A[surviving_species, surviving_species],
  Niche = Niche[surviving_species, ]
)

# Spatio-temporal  --------------------------------------------------------
# Several patches over several time-points
time_indices = 1001:1101#seq(201, by = 5, length.out = 201)
spatial_temporal = data_wo_burnin[data_wo_burnin$patch %in% 1:10,]
spatial_temporal = spatial_temporal[spatial_temporal$time %in% 1001:1101,]
X = spatial_temporal %>% select(env, time, x, y)
Y = spatial_temporal %>% select(-x, -y, -env)
surviving_species = which(colSums(Y[,-c(1:2)]) > 0, arr.ind = TRUE)
Y = Y[,c(1, 2, surviving_species+2) ]

X_corrected = cbind(X[-c(1:10),], Y[-c(1001:1010),]) %>% select(-time, -patch)
X_corrected = cbind(X_corrected[,c(1:3)], scale(X_corrected[,-c(1:3)]+0.0001))

X_corrected= data.frame(X_corrected)
colnames(X_corrected) = c("env","x", "y" ,paste0("sp_", 1:length(surviving_species)))
form = paste0("~env + I(env^2) + ", paste0( paste0( "sp_", 1:length(surviving_species) ), collapse = " + ") )

coords = X_corrected %>% select(x, y) %>% scale %>% as.data.frame()
colnames(coords) = c("x", "y")
A_abund = A[surviving_species, surviving_species]
Niche_abund = Niche[surviving_species,]
# Train model - Poisson
spatial_temporal_poisson = sjSDM(Y[-c(1:10),] %>% select(-time, -patch) %>% as.matrix(), 
                                 linear(X_corrected, form), 
                                 spatial = linear(coords, ~ x + y + x:y + I(x^2) + I(y^2)),
                                 family = poisson(),
                                 device = 0)

# Train model - PA
X = spatial_temporal %>% select(env, time, x, y)
Y = spatial_temporal %>% select(-x, -y, -env)
PA = apply(Y[,-c(1:2)], 1:2, function(P) ifelse(P>0, 1, 0))
surviving_species = which((colSums(PA) > 0) & (colSums(PA) < nrow(PA)), arr.ind = TRUE)
Y = Y[,c(1, 2, surviving_species+2) ]
PA = cbind(Y[,1:2], apply(Y[,-c(1:2)], 1:2, function(P) ifelse(P>0, 1, 0) ))
colnames(PA) = colnames(Y)

X_corrected = cbind(X[-c(1:10),], PA[-c(1001:1010),]) %>% select(-time, -patch)
X_corrected = cbind(X_corrected[,c(1:3)], scale(X_corrected[,-c(1:3)])+0.0001)

X_corrected= data.frame(X_corrected)
colnames(X_corrected) = c("env","x", "y" ,paste0("sp_", 1:length(surviving_species)))
form = paste0("~env + I(env^2) + ", paste0( paste0( "sp_", 1:length(surviving_species) ), collapse = " + ") )

coords = X_corrected %>% select(x, y) %>% scale %>% as.data.frame()
colnames(coords) = c("x", "y")

spatial_temporal_PA= sjSDM(PA[-c(1:10),] %>% select(-time, -patch) %>% as.matrix(), 
                           linear(X_corrected, form), 
                           spatial = linear(coords, ~ x + y + x:y + I(x^2) + I(y^2)),
                           family = binomial(),
                           device = 0)



# Prepare outputs - Abundance
beta_env = t(coef(spatial_temporal_poisson)$env[[1]])[1:3,]
rownames(beta_env) = colnames(spatial_temporal_poisson$data$X)[1:3]
beta_space = t(coef(spatial_temporal_poisson)$spatial[[1]])
rownames(beta_space) = colnames(spatial_temporal_poisson$spatial$X)
A_species = t(coef(spatial_temporal_poisson)$env[[1]])[-(1:3),]
rownames(A_species) = colnames(spatial_temporal_poisson$data$X)[-(1:3)]

results_temporal_abundance = list(
  model = spatial_temporal_poisson,
  data = list(Y = spatial_temporal_poisson$data$Y,
              env = spatial_temporal_poisson$data$X[,1:3],
              biotic = spatial_temporal_poisson$data$X[,-(1:3)]),
  beta_env = beta_env,
  beta_space = beta_space,
  A_species = A_species,
  cooccurrence = cov2cor(getCov(spatial_temporal_poisson)),
  link = function(P) exp(P)
)


# Prepare outputs - PA
beta_env = t(coef(spatial_temporal_PA)$env[[1]])[1:3,]
rownames(beta_env) = colnames(spatial_temporal_PA$data$X)[1:3]
beta_space = t(coef(spatial_temporal_PA)$spatial[[1]])
rownames(beta_space) = colnames(spatial_temporal_PA$spatial$X)
A_species = t(coef(spatial_temporal_PA)$env[[1]])[-(1:3),]
rownames(A_species) = colnames(spatial_temporal_PA$data$X)[-(1:3)]

results_temporal_abundance = list(
  model = spatial_temporal_PA,
  data = list(Y = spatial_temporal_PA$data$Y,
              env = spatial_temporal_PA$data$X[,1:3],
              biotic = spatial_temporal_PA$data$X[,-(1:3)]),
  beta_env = beta_env,
  beta_space = beta_space,
  A_species = A_species,
  cooccurrence = cov2cor(getCov(spatial_temporal_PA)),
  link = function(P) 1/(1+exp(-P))
)

results_spatial_temporal = list(
  results_abundance = results_temporal_abundance,
  results_PA = results_temporal_PA,
  A_PA = A[surviving_species, surviving_species],
  Niche_PA = Niche[surviving_species, ],
  A_abund = A_abund,
  Niche_abund = Niche_abund
)

save(results_spatial, results_temporal, results_spatial_temporal, file = "results/01_prelim_results.RData")
