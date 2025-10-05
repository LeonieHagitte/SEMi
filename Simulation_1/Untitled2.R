
library(lavaan)
    
    # Flexible SEM data generator with typ_des_models
    generate_data <- function(sample_size = 100,
                              typ_des_models = 1,
                              sample_type = FALSE,
                              loading_offset = 0,
                              moderator_value = 0,
                              effect_size = 0,
                              resid_vars = rep(1, 4),      # residual variances
                              latent_var = 1,              # latent variance
                              intercepts = rep(0, 4)) {    # intercepts of items
      
      # --- sample model type ---
        if (sample_type) {
        typ_des_models <- sample(1:3, size = 1)  # randomly choose 1, 2, or 3
        message("Sampled typ_des_models: ", typ_des_models)
      }
      
      # --- Set loadings based on typ_des_models ---
      if (typ_des_models == 1) {
        loading1 <- 0.7
        loading2 <- 0.7
        loading3 <- 0.7
        loading4 <- 1
      } else if (typ_des_models == 2) {
        loading1 <- loading_offset + moderator_value * effect_size
        loading2 <- 0.7
        loading3 <- 0.7
        loading4 <- 1
      } else if (typ_des_models == 3) {
        loading1 <- loading_offset + moderator_value * effect_size
        loading2 <- loading_offset + moderator_value * effect_size
        loading3 <- 0.7
        loading4 <- 1
      } else {
        stop("typ_des_models must be 1, 2, or 3")
      }
    
      model1_metric <- paste0(
        'eta =~ ', loading1,'*x1 + ', loading2,'*x2 + ', loading3,'*x3 + ', loading4,'*x4\n',
        'x1 ~~ x1\nx2 ~~ x2\nx3 ~~ x3\nx4 ~~ x4\n',
        'eta ~~ eta\n',
        'x1 ~ 0*1\nx2 ~ 0*1\nx3 ~ 0*1\nx4 ~ 0*1\n',
        'eta ~ 0*1\n'
      )
      
simdat <- lavaan::simulateData(model1_metric, sample.nobs = sample_size)
return(simdat)

    }
    
    # Example usage:
    set.seed(123)
    sim_data <- generate_data(typ_des_models = 2, loading_offset = 0.5, 
                              moderator_value = 1, effect_size = 0.2, 
                              sample_size = 10)
head(sim_data)
  