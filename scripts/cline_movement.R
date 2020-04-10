# load libraries
library(ggplot2); library(viridis); library(mgcv); 
library(boot); library(magrittr); library(tidymv); 
library(data.table); library(purrr); library(plyr);
library(cowplot)

# load simulation function (modified from hzam to export cline width paramns)
run_simulation <- function(neutral_loci=3, pop1_range_limit_left=0, pop1_range_limit_right=0.48, pop2_range_limit_left=0.52,
                           pop2_range_limit_right=1, K_half=800, sigma_comp=0.1, range_limit_left=0,  range_limit_right=1, 
                           pref_ratio=0.5, max_generations=10, per_reject_cost=0, meandispersal=0.01, growth_rate=1.05,
                           beginning_columns=1, survival_fitness_method=1, beta=1, male_trait_loci=1, assort_mating=1, 
                           female_trait_loci=1, underdominant_loci=1, same_loci_MTL_UDL=TRUE,
                           mating_trait_dominance=1, dom_coefficient=0.75, hybrid_fitness=1, ngen_fit=5){
  
  pref_SD <- sqrt( -1 / (2 * log(pref_ratio)))     # width of female acceptance curve for male trait
  s_per_locus <- 1 - hybrid_fitness^(1/underdominant_loci)  # loss in fitness due to each heterozygous locus
  
  # each locus has two alleles, so takes up two columns
  # for clarity later, calculate start and end columns for each locus type
  MTL_col_start <- beginning_columns + 1
  MTL_col_end <- beginning_columns + 2*male_trait_loci
  if (assort_mating == 1) {
    FTL_col_start <- MTL_col_start
    FTL_col_end <- MTL_col_end
  }
  if (assort_mating == 0) {
    FTL_col_start = MTL_col_end + 1
    FTL_col_end = MTL_col_end + 2*female_trait_loci
  }
  if (same_loci_MTL_UDL == TRUE) {
    UDL_col_start <- FTL_col_start
    UDL_col_end <- FTL_col_end
  }
  if (same_loci_MTL_UDL == FALSE) {
    UDL_col_start <- FTL_col_end + 1
    UDL_col_end <- FTL_col_end + 2*underdominant_loci
  }
  NL_col_start <- UDL_col_end + 1
  NL_col_end <- UDL_col_end + 2*neutral_loci
  total_loci <- (NL_col_end - beginning_columns)/2
  
  # Set up starting conditions
  
  pop1_starting_N <- round(K_half * (pop1_range_limit_right - pop1_range_limit_left))  # remember that K_half is for each sex, so whole pop is 2K, divided among the two species
  pop2_starting_N <- round(K_half * (pop2_range_limit_right - pop2_range_limit_left))
  starting_N <- pop1_starting_N + pop2_starting_N
  
  # generate the population of females
  pop_matrix_F = matrix(-9, nrow=starting_N, ncol=1+2*total_loci) # create matrix to store population locations and genotypes; columns in this order: location, genotype columns
  # generate starting values for pop1:
  pop_matrix_F[1:pop1_starting_N,1] <- runif(n=pop1_starting_N, min=pop1_range_limit_left, max=pop1_range_limit_right)  # assigns random locations for pop1
  pop_matrix_F[1:pop1_starting_N,2:(1+2*total_loci)] <- 0   # assigns genotypes of pop1
  # generate starting values for pop2:
  pop_matrix_F[(1+pop1_starting_N):starting_N,1] <- runif(n=pop2_starting_N, min=pop2_range_limit_left, max=pop2_range_limit_right)  # assigns random locations for pop2
  pop_matrix_F[(1+pop1_starting_N):starting_N,2:(1+2*total_loci)] <- 1   # assigns genotypes of pop2
  
  # generate the population of males
  pop_matrix_M = matrix(-9, nrow=starting_N, ncol=1+2*total_loci) # create matrix to store population locations and genotypes; columns in this order: location, genotype columns
  # generate starting values for pop1:
  pop_matrix_M[1:pop1_starting_N,1] <- runif(n=pop1_starting_N, min=pop1_range_limit_left, max=pop1_range_limit_right)  # assigns random locations for pop1
  pop_matrix_M[1:pop1_starting_N,2:(1+2*total_loci)] <- 0   # assigns genotypes of pop1
  # generate starting values for pop2:
  pop_matrix_M[(1+pop1_starting_N):starting_N,1] <- runif(n=pop2_starting_N, min=pop2_range_limit_left, max=pop2_range_limit_right)  # assigns random locations for pop2
  pop_matrix_M[(1+pop1_starting_N):starting_N,2:(1+2*total_loci)] <- 1   # assigns genotypes of pop2
  
  # set up locations every 0.001 across range, and calculate density at each when 2K individuals perfectly spaced:
  spaced_locations <- round(seq(from = range_limit_left, to = range_limit_right, length.out = 1001), digits=3)  # the round is needed to ensure only 3 decimal places, avoiding error later
  ind_locations_if_even_at_K <- seq(from = range_limit_left, to = range_limit_right, length.out = 2*K_half)
  get_density_if_even_at_K <- function(focal_location) {
    return(sum(exp(-((ind_locations_if_even_at_K - focal_location)^2)/(2*(sigma_comp^2)))))
  } 
  ideal_densities_at_spaced_locations <- sapply(spaced_locations, get_density_if_even_at_K)
  
  # set up list for cline data
  pop_matrix_C <- c()
  
  # Run cycles of mate choice, reproduction, dispersal
    for (time in 1:max_generations) { # start of generation loop
      print(paste("generation ",time))
      
      # New way of calculating density dependence, using same method as Irwin 2002 AmNat:
      ind_locations_real <- c(pop_matrix_F[,1] , pop_matrix_M[,1])
      get_density_real <- function(focal_location) {
        return(sum(exp(-((ind_locations_real - focal_location)^2)/(2*(sigma_comp^2)))))
      } 
      real_densities_at_spaced_locations <- sapply(spaced_locations, get_density_real)
      # Liou&Price equation, where density dependence based on same K everywhere, but limited growth rate:
      local_growth_rates <- growth_rate*ideal_densities_at_spaced_locations / (ideal_densities_at_spaced_locations + ((real_densities_at_spaced_locations)*(growth_rate - 1)))
      
      if (mating_trait_dominance==1) {  # the dominant case
        fraction_trait_loci_dominant <- function(x) {
          dominant_count <- 0
          # cycle through trait loci and count up the number of loci with a dominant (="1") allele
          for (locus_count in 1:male_trait_loci) {
            locus_columns <- (2*(locus_count-1))+(MTL_col_start:((MTL_col_start)+1))
            if (sum(x[locus_columns]) %in% (1)) {
              dominant_count <- dominant_count + dom_coefficient
            }
            if (sum(x[locus_columns]) %in% (2)) {
              dominant_count <- dominant_count + 1
            }
          }
          return(dominant_count / male_trait_loci)  
        } 
        male_traits <- apply(pop_matrix_M, 1, fraction_trait_loci_dominant) # the dominant case
      } 
      if (mating_trait_dominance==2) {
        male_traits <- apply(pop_matrix_M, 1, function(x) mean(x[MTL_col_start:MTL_col_end]))  # the additive case
      }
      
      if (mating_trait_dominance==1) {  # the dominant case
        fraction_trait_loci_dominant <- function(x) {
          dominant_count <- 0
          # cycle through trait loci and count up the number of loci with a dominant (="1") allele
          for (locus_count in 1:female_trait_loci) {
            locus_columns <- (2*(locus_count-1))+(FTL_col_start:((FTL_col_start)+1))
            if (sum(x[locus_columns]) %in% (1)) {
              dominant_count <- dominant_count + dom_coefficient
            }
            if (sum(x[locus_columns]) %in% (2)) {
              dominant_count <- dominant_count + 1
            }
          }
          return(dominant_count / female_trait_loci)  
        } 
        female_traits <- apply(pop_matrix_F, 1, fraction_trait_loci_dominant) # the dominant case
      } 
      if (mating_trait_dominance==2) {
        female_traits <- apply(pop_matrix_F, 1, function(x) mean(x[FTL_col_start:FTL_col_end]))  # the additive case
      }
      
      # cycle through the mothers, determining numbers of offspring
      daughters_per_mother_list <- vector("list", nrow(pop_matrix_F)) # pop_matrix_daughters <- matrix(nrow=0, ncol=1+2*total_loci)
      sons_per_mother_list <- vector("list", nrow(pop_matrix_F)) # pop_matrix_sons <- matrix(nrow=0, ncol=1+2*total_loci)
      # now cycle through all females and determine male mate, reproduction, and daughter genotypes and locations
      for (mother in 1:nrow(pop_matrix_F))  {
        # Find a male mate
        # rank males in terms of distance from female
        male_dists <- abs(pop_matrix_F[mother,1] - pop_matrix_M[,1])  # calculates vector of dist of each male from focal female
        # pick closest male, determine his male signalling trait
        mate <- 0
        rejects <- 0
        while (mate == 0) {
          focal_male <- which(male_dists == min(male_dists))  # returns row number of male that is closest
          
          # compare male trait with female's trait (preference), and determine
          # whether she accepts; note that match_strength is determined by a
          # Gaussian, with a maximum of 1 and minimum of zero.
          match_strength <- (exp(1) ^ ((-(male_traits[focal_male] - female_traits[mother])^2) / (2 * (pref_SD ^2))))
          
          if (runif(1) < match_strength) {
            # she accepts male, and they reproduce
            father <- focal_male
            mate <- 1
          } else {
            # she rejects male
            # change that male's distance to a very large number (99) so he
            # won't be considered again (just in this bit)
            male_dists[focal_male] <- 99
            rejects <- rejects + 1  # count number of rejects, for imposition of fitness cost for search time
          }
        }
        
        # Reproduce
        
        # determine fitness cost due to mate search time
        search_fitness <- (1-per_reject_cost) ^ rejects    # calculate proportion fitness lost (1-cost) due to mate search time
        
        # determine local growth rate at mother's location:
        location_rounded <- round(pop_matrix_F[mother,1], digits=3)
        local_growth <- local_growth_rates[spaced_locations==location_rounded]
        
        #combine for total fitness:   
        
        reproductive_fitness <- 2*local_growth * search_fitness  # the 2 is because only females, not males, produce offspring
        # now draw the number of offspring from a poisson distribution with a mean of total_fitness
        offspring <- rpois(1, reproductive_fitness) 
        daughters <- matrix(nrow=0, ncol=1+2*total_loci)
        sons <- matrix(nrow=0, ncol=1+2*total_loci)
        if (offspring >= 1)  { # if offspring, generate their location and genotypes
          for (kid in 1:offspring) {
            kid_info <- rep(-9, 1+2*total_loci)
            while (1 == 1)  {	# disperse the kid
              newloc <- pop_matrix_F[mother,1] + rnorm(1, mean=0, sd=meandispersal)
              if ((newloc <= range_limit_right) & (newloc >= range_limit_left))  {
                break
              }
            }
            kid_info[1] <- newloc
            # generate genotypes; for each locus, first column for allele from mother, second for allele from father
            for (locus in 1:total_loci) {   
              # choose allele from mother
              kid_info[2+2*(locus-1)] <- pop_matrix_F[mother,1+2*(locus-1)+sample(2, 1)]
              # choose allele from father
              kid_info[3+2*(locus-1)] <- pop_matrix_M[father,1+2*(locus-1)+sample(2, 1)]
            }
            # determine sex of kid and add to table
            if (runif(1) > 0.5) {
              daughters <- rbind(daughters, kid_info)
              # pop_matrix_daughters <- rbind(pop_matrix_daughters, kid_info) 
            } else {
              sons <- rbind(sons, kid_info)
              # pop_matrix_sons <- rbind(pop_matrix_sons, kid_info)
            }
          } 
        } 
        daughters_per_mother_list[[mother]] <- daughters
        sons_per_mother_list[[mother]] <- sons
      } 
      pop_matrix_daughters <- do.call("rbind", daughters_per_mother_list)
      pop_matrix_sons <- do.call("rbind", sons_per_mother_list)
      
      # determine fitnesses of daughters due to heterozygosity at (option 1) heterozygosity, or (option 2) epistasis 
      if (survival_fitness_method == 1) {
        underdominance_fitness_daughters <- rep(NaN, dim(pop_matrix_daughters)[1])
        for (daughter in 1:nrow(pop_matrix_daughters))  {
          heterozyg_loci <- 0
          for (locus in 1:underdominant_loci) {
            locus_col <- UDL_col_start + 2*(locus-1)
            if (mean(pop_matrix_daughters[daughter,locus_col:(locus_col+1)]) == 0.5 ) { 
              heterozyg_loci <- heterozyg_loci + 1
            }
          }
          underdominance_fitness_daughters[daughter] <- (1-s_per_locus) ^ heterozyg_loci
        }
        
        # determine whether each daughter survives to adulthood
        random_proportions_daughters <- runif(n=length(underdominance_fitness_daughters), min=0, max=1)
        daughters_survive <- underdominance_fitness_daughters > random_proportions_daughters
        
        # determine fitnesses of sons due to heterozygosity at underdominance loci
        underdominance_fitness_sons <- rep(NaN, dim(pop_matrix_sons)[1])
        for (son in 1:nrow(pop_matrix_sons))  {
          heterozyg_loci <- 0
          for (locus in 1:underdominant_loci) {
            locus_col <- UDL_col_start + 2*(locus-1)
            if (mean(pop_matrix_sons[son,locus_col:(locus_col+1)]) == 0.5 ) { 
              heterozyg_loci <- heterozyg_loci + 1
            }
          }
          underdominance_fitness_sons[son] <- (1-s_per_locus) ^ heterozyg_loci
        } 
        # determine whether each son survives to adulthood
        random_proportions_sons <- runif(n=length(underdominance_fitness_sons), min=0, max=1)
        sons_survive <- underdominance_fitness_sons > random_proportions_sons
        
      } else if (survival_fitness_method == 2) {
        HI_daughters <- rowMeans(pop_matrix_daughters[,UDL_col_start:UDL_col_end]) # apply(pop_matrix_daughters, 1, function(x) mean(x[UDL_col_start:UDL_col_end]))
        epistasis_fitness_daughters <- 1 - (1-hybrid_fitness) * (4*HI_daughters*(1-HI_daughters))^beta
        random_proportions_daughters <- runif(n=length(epistasis_fitness_daughters), min=0, max=1)
        daughters_survive <- epistasis_fitness_daughters > random_proportions_daughters
        
        HI_sons <- rowMeans(pop_matrix_sons[,UDL_col_start:UDL_col_end])  # apply(pop_matrix_sons, 1, function(x) mean(x[UDL_col_start:UDL_col_end]))
        epistasis_fitness_sons <- 1 - (1-hybrid_fitness) * (4*HI_sons*(1-HI_sons))^beta
        random_proportions_sons <- runif(n=length(epistasis_fitness_sons), min=0, max=1)
        sons_survive <- epistasis_fitness_sons > random_proportions_sons  
      }
      
      # assign surviving offspring to new adult population
      pop_matrix_F <- pop_matrix_daughters[daughters_survive,]
      pop_matrix_M <- pop_matrix_sons[sons_survive,]
      
      pop_matrix <- list(pop_matrix_F,pop_matrix_M)
      
      # modulo to fit cline params every ngen_fit
      if (time%%ngen_fit == 0) {
        sim_output_1 <- as.data.frame(pop_matrix[[1]])
        sim_output_2 <- as.data.frame(pop_matrix[[2]])
        location <- c(sim_output_1[,1],sim_output_2[,1])
        location <- rep(location, 2)
        neutral_freqs <- as.vector(c(rowMeans(sim_output_1[,(2+(2*male_trait_loci)):ncol(sim_output_1)]), 
                                     rowMeans(sim_output_2[,(2+(2*male_trait_loci)):ncol(sim_output_2)])))
        trait_freqs <- as.vector(c(rowMeans(sim_output_1[,2:(1+(2*male_trait_loci))]), 
                                   rowMeans(sim_output_2[,2:(1+(2*male_trait_loci))])))
        freqs <- c(neutral_freqs, trait_freqs)
        vars <- c(rep("neutral",length(neutral_freqs)), rep("traits",length(trait_freqs)))
        format_output <- cbind.data.frame(location, freqs, vars)
        
        # calculate neutral loci cline width and center
        model_n <- gam(freqs ~ s(location), data=subset(format_output,vars=="neutral"), quasibinomial(link = "logit"),
                       method = "P-ML")
        model_np <- predict_gam(model_n, length_out = 1000)
        model_np$fit <- inv.logit(model_np$fit)
        colnames(model_np) <- c("location", "freqs", "se.fit")
        left_width_margin_n <- max(model_np[model_np$freqs<=0.1,]$location)
        right_width_margin_n <- min(model_np[model_np$freqs>=0.9,]$location)
        dist_from_HI50percent_n <- abs(model_np$freqs-0.5)
        center_n <- model_np[dist_from_HI50percent_n == min(dist_from_HI50percent_n),]$location
        width_n <- right_width_margin_n - left_width_margin_n
        
        # calculate mating trait cline width and center
        model_t <- gam(freqs ~ s(location), data=subset(format_output,vars=="traits"), quasibinomial(link = "logit"),
                       method = "P-ML")
        model_tp <- predict_gam(model_t, length_out = 10000)
        model_tp$fit <- inv.logit(model_tp$fit)
        colnames(model_tp) <- c("location", "freqs", "se.fit")
        left_width_margin_t <- max(model_tp[model_tp$freqs<=0.1,]$location)
        right_width_margin_t <- min(model_tp[model_tp$freqs>=0.9,]$location)
        dist_from_HI50percent_t <- abs(model_tp$freqs-0.5)
        center_t <- model_tp[dist_from_HI50percent_t == min(dist_from_HI50percent_t),]$location
        width_t <- right_width_margin_t - left_width_margin_t
        
        # export df
        neutral_df <- cbind.data.frame("Neutral loci", center_n, width_n, time)
        colnames(neutral_df) <- c("cline_type","cline_center","cline_width", "generation")
        trait_df <- cbind.data.frame("Mating trait loci", center_t, width_t, time)
        colnames(trait_df) <- c("cline_type","cline_center","cline_width", "generation")
        data_f <- rbind(neutral_df, trait_df)
        
        # print table as gut-check
        print(data_f)
        
        pop_matrix_C[[time]] <- data_f
      } 
    } # end of generation loop
  pop_matrix <- list(pop_matrix, pop_matrix_C)
  return(pop_matrix)
}

# function to grab cline fit data from final generation
grab_cline_data <- function(list, nrep=50, ngen=500){
  tmp <- list()
  for(i in 1:length(list)){
    tmp[[i]] <-list[[i]][2]
  }
  tmp <- unlist(tmp,recursive=FALSE)
  tmp <- unlist(tmp,recursive=FALSE)
  tmp <- Filter(Negate(is.null), tmp)
  tmp <- lapply(tmp, function(x) return(x[1:2,])) # this drops extra rows from error
  tmp <- rbindlist(tmp)
  tmp$replicate <- rep(1:nrep, each=(ngen*2)/5)
  tmp <- tmp[tmp$generation==500,]
  return(tmp)
}

# function to grab population data
grab_pop_matrices <- function(list){
  tmp <- list()
  for(i in 1:length(list)){
    tmp[[i]] <-replicates_6[[i]][1]
  }
  tmp <- unlist(tmp,recursive=FALSE)
  tmp <- unlist(tmp,recursive=FALSE)
  obj <- list()
  for(x in 1:length(tmp)){
    print(x)
    if(x%%2==1){
      obj[[x]] <- rbind.data.frame(tmp[[x]],tmp[[x+1]])
    } 
  }
  obj <- Filter(Negate(is.null), obj)
  return(obj)
}


# function to calculate when cline moved left
left_movement <- function(cline_500){
  count <- 0
  for(i in 1:max(cline_500$replicate)){
    tmp <- cline_500[cline_500$replicate==i,]
    if(tmp$cline_center[1]>tmp$cline_center[2]){
      count <- count + 1
    }
  }
  left <- count/max(cline_500$replicate) 
  return(left)
}


# function to calculate when the mating trait cline was more narrow than the genome-wide cline
trait_narrower <- function(cline_500){
  count_2 <- 0
  for(i in 1:max(cline_500$replicate)){
    tmp <- cline_500[cline_500$replicate==i,]
    if(tmp$cline_width[1]>tmp$cline_width[2]){
      count_2 <- count_2 + 1
    }
  }
  narrow <- count_2/max(cline_500$replicate)
  return(narrow)
}

# function to do a wilcox ranked sign test on cline width—"cline  500" = cline output from 500th
# generation (i.e. "clines_1" etc in this script)
wilcox_cline_width <- function(cline_500){
  neut_w <- cline_500[cline_500$cline_type=="Neutral loci",]$cline_width
  pheno_w <- cline_500[cline_500$cline_type=="Mating trait loci",]$cline_width
  df <- cbind.data.frame(neut_w, pheno_w)
  df <- df[!is.infinite(rowSums(df)),]
  neut_w <- df$neut_w
  pheno_w <- df$pheno_w
  wilcox.test(neut_w, pheno_w, conf.int = T, paired=T) # W = 1641, p-value = 0.007102
}

# function to do a wilcox ranked sign test on cline center
wilcox_cline_center <- function(cline_500){
  neut <- cline_500[cline_500$cline_type=="Neutral loci",]$cline_center
  pheno <- cline_500[cline_500$cline_type=="Mating trait loci",]$cline_center 
  wilcox.test(neut, pheno, conf.int = T,paired=T) 
}

cline_stats <- function(cline_500){
  # mean cline center
  mean_center <- ddply(cline_500, "cline_type", summarise, mean_center=mean(cline_center))
  # cline center variance
  sd_center <- ddply(cline_500, "cline_type", summarise, sd_center=sd(cline_center))
  # mean cline width
  mean_width <- ddply(cline_500[cline_500$cline_width<Inf,], "cline_type", summarise, mean_width=mean(cline_width, na.rm = T))
  # cline width variance
  sd_width <- ddply(cline_500[cline_500$cline_width<Inf,], "cline_type", summarise, sd_center=sd(cline_width))
  df <- cbind.data.frame(mean_center, sd_center[,2], mean_width[,2], sd_width[,2])
  colnames(df) <- c("cline_type", "mean_cline_center", "sd_cline_center", "mean_cline_width", "sd_cline_width")
  return(df)
}

# run fifty replicate simulations w/ strong partial dominance, moderate assortative mating
replicates_1 <- list()
for(rep in 1:50){
  replicates_1[[rep]] <- run_simulation(max_generations = 500, dom_coefficient = 0.75, 
                                      pref_ratio=0.5, male_trait_loci = 2)
}

# strong partial dominance, weak assortative mating
replicates_2 <- list()
for(rep in 1:50){
  replicates_2[[rep]] <- run_simulation(max_generations = 500, dom_coefficient = 0.75, 
                                        pref_ratio=0.75, male_trait_loci = 2)
}

# moderate partial dominance, moderate assortative making
replicates_3 <- list()
for(rep in 1:50){
  replicates_3[[rep]] <- run_simulation(max_generations = 500, dom_coefficient = 0.5, 
                                        pref_ratio=0.5, male_trait_loci = 2)
}

# moderate partial dominance, weak assortative mating
replicates_4 <- list()
for(rep in 1:50){
  replicates_4[[rep]] <- run_simulation(max_generations = 500, dom_coefficient = 0.5, 
                                        pref_ratio=0.75, male_trait_loci = 2)
}


# strong partial dominance, moderate assortative mating, one locus 
replicates_5 <- list()
for(rep in 1:50){
  replicates_5[[rep]] <- run_simulation(max_generations = 500, dom_coefficient = 0.75, 
                                        pref_ratio=0.5, male_trait_loci = 1)
}

#  strong partial dominance, weak assortative mating, one locus
replicates_6 <- list()
for(rep in 1:50){
  replicates_6[[rep]] <- run_simulation(max_generations = 500, dom_coefficient = 0.75, 
                                        pref_ratio=0.75, male_trait_loci = 1)
}

# moderate partial dominance, moderate assortative mating, one locus
replicates_7 <- list()
for(rep in 1:50){
  replicates_7[[rep]] <- run_simulation(max_generations = 500, dom_coefficient = 0.5, 
                                        pref_ratio=0.5, male_trait_loci = 1)
}

# moderate partial dominance, weak assortative mating, one locus
replicates_8 <- list()
for(rep in 1:50){
  replicates_8[[rep]] <- run_simulation(max_generations = 500, dom_coefficient = 0.5, 
                                        pref_ratio=0.75, male_trait_loci = 1)
}

# extract cline data—numbers correlate to simulation
clines_1 <- grab_cline_data(replicates_1, 50, 500) 
clines_2 <- grab_cline_data(replicates_2, 50, 500)
clines_3 <- grab_cline_data(replicates_3, 50, 500) 
clines_4 <- grab_cline_data(replicates_4, 50, 500) 
clines_5 <- grab_cline_data(replicates_5, 50, 500) 
clines_6 <- grab_cline_data(replicates_6, 50, 500) 
clines_7 <- grab_cline_data(replicates_7, 50, 500)
clines_8 <- grab_cline_data(replicates_8, 50, 500) 

# write data to file
write.csv(clines_1, "~/Dropbox/wagtails/data/cline_500/clines_1.csv")
write.csv(clines_2, "~/Dropbox/wagtails/data/cline_500/clines_2.csv")
write.csv(clines_3, "~/Dropbox/wagtails/data/cline_500/clines_3.csv")
write.csv(clines_4, "~/Dropbox/wagtails/data/cline_500/clines_4.csv")
write.csv(clines_5, "~/Dropbox/wagtails/data/cline_500/clines_5.csv")
write.csv(clines_6, "~/Dropbox/wagtails/data/cline_500/clines_6.csv")
write.csv(clines_7, "~/Dropbox/wagtails/data/cline_500/clines_7.csv")
write.csv(clines_8, "~/Dropbox/wagtails/data/cline_500/clines_8.csv")

# assess how often mating trait cline moves towards recessive homozygote
left_movement(clines_1) # 0.58
left_movement(clines_2) # 0.58
left_movement(clines_3) # 0.44—meaning rightward movement toward dominant homozygote
left_movement(clines_4) # 0.52
left_movement(clines_5) # 0.72
left_movement(clines_6) # 0.62
left_movement(clines_7) # 0.52
left_movement(clines_8) # 0.5

# assess how often mating trait clines are narrower than genome-wide cline
trait_narrower(clines_1) # 0.6
trait_narrower(clines_2) # 0.66
trait_narrower(clines_3) # 0.62
trait_narrower(clines_4) # 0.64
trait_narrower(clines_5) # 1
trait_narrower(clines_6) # 0.96
trait_narrower(clines_7) # 0.98
trait_narrower(clines_8) # 0.94

# cline statistics
cline_stats(clines_1)

#          cline_type mean_cline_center sd_cline_center mean_cline_width sd_cline_width
# 1      Neutral loci         0.5154811       0.1275102        0.4759701      0.1619554
# 2 Mating trait loci         0.4190672       0.3073818        0.2249119      0.1041093

# cline statistics
cline_stats(clines_2)

#          cline_type mean_cline_center sd_cline_center mean_cline_width sd_cline_width
# 1      Neutral loci         0.4890879       0.1054418        0.4241969      0.1452947
# 2 Mating trait loci         0.4401698       0.2122520        0.2710710      0.1376481

cline_stats(clines_3)

#          cline_type mean_cline_center sd_cline_center mean_cline_width sd_cline_width
# 1      Neutral loci         0.5139281       0.1120579       0.47226796      0.1813014
# 2 Mating trait loci         0.5076960       0.2919870       0.07299187      0.3423503

cline_stats(clines_4)

#          cline_type mean_cline_center sd_cline_center mean_cline_width sd_cline_width
# 1      Neutral loci         0.5137647       0.1352037        0.4411549      0.1468363
# 2 Mating trait loci         0.5063718       0.2486510        0.2756875      0.1464122

cline_stats(clines_5)
#         cline_type mean_cline_center sd_cline_center mean_cline_width sd_cline_width
# 1      Neutral loci         0.5081300       0.1164983       0.43856613     0.14122698
# 2 Mating trait loci         0.4281591       0.1575376       0.09935873     0.05010819

cline_stats(clines_6)
#          cline_type mean_cline_center sd_cline_center mean_cline_width sd_cline_width
# 1      Neutral loci         0.4797799       0.1446467        0.4395531     0.13654915
# 2 Mating trait loci         0.4272766       0.1717500        0.1218575     0.07373175

cline_stats(clines_7)
#         cline_type mean_cline_center sd_cline_center mean_cline_width sd_cline_width
#1      Neutral loci         0.4965957       0.1122563       0.43854900     0.16226342
#2 Mating trait loci         0.4920439       0.1555197       0.09399186     0.04581947

cline_stats(clines_8)
#           cline_type mean_cline_center sd_cline_center mean_cline_width sd_cline_width
# 1      Neutral loci         0.4973354       0.1157160        0.4307719      0.1530064
# 2 Mating trait loci         0.4868155       0.1718089        0.1208503      0.1216942


# see whether distribution of cline centers differs across simulations
wilcox_cline_center(clines_1) # p-value = 0.04892
wilcox_cline_center(clines_2) # p-value = 0.1645
wilcox_cline_center(clines_3) # p-value = 0.9923
wilcox_cline_center(clines_4) # p-value = 0.7137
wilcox_cline_center(clines_5) # p-value = 0.0006325
wilcox_cline_center(clines_6) # p-value = 0.06382
wilcox_cline_center(clines_7) # p-value = 0.8168
wilcox_cline_center(clines_8) # p-value = 0.7574

# see whether distribution of cline widths differs across simulations
wilcox_cline_width(clines_1) # p-value = 1.276e-07
wilcox_cline_width(clines_2) # p-value = 0.0001594
wilcox_cline_width(clines_3) # p-value = 4.712e-07
wilcox_cline_width(clines_4) # p-value = 0.0004383
wilcox_cline_width(clines_5) # p-value = 1.819e-12
wilcox_cline_width(clines_6) # p-value = 1.819e-12
wilcox_cline_width(clines_7) # p-value = 2.274e-13
wilcox_cline_width(clines_8) # p-value = 1.819e-12

# extract population matrices
pop_matrix_1 <- grab_pop_matrices(replicates_1)
pop_matrix_2 <- grab_pop_matrices(replicates_2)
pop_matrix_3 <- grab_pop_matrices(replicates_3)
pop_matrix_4 <- grab_pop_matrices(replicates_4)
pop_matrix_5 <- grab_pop_matrices(replicates_5)
pop_matrix_6 <- grab_pop_matrices(replicates_6)
pop_matrix_7 <- grab_pop_matrices(replicates_7)
pop_matrix_8 <- grab_pop_matrices(replicates_8)

