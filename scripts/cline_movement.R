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


# run fifty replicate simulations
replicates <- list()
for(rep in 1:50){
  replicates[[rep]] <- run_simulation(max_generations = 500)
}

# wrangle population data
pop_matrices <- lapply(replicates, `[[`, 1)

# extract representative example population matrices
num_trait_loci <- 1
#plot_example <- pop_matrices[[47]]
#plot_example <- pop_matrices[[17]]
plot_example <- pop_matrices[[44]]
sim_output_1 <- as.data.frame(plot_example[[1]])
sim_output_2 <- as.data.frame(plot_example[[2]])
location <- c(sim_output_1[,1],sim_output_2[,1])
location <- rep(location, 2)
neutral_freqs <- as.vector(c(rowMeans(sim_output_1[,(2+(2*num_trait_loci)):ncol(sim_output_1)]), 
                             rowMeans(sim_output_2[,(2+(2*num_trait_loci)):ncol(sim_output_2)])))
trait_freqs <- as.vector(c(rowMeans(sim_output_1[,2:(1+(2*num_trait_loci))]), 
                           rowMeans(sim_output_2[,2:(1+(2*num_trait_loci))])))
freqs <- c(neutral_freqs, trait_freqs)
vars <- c(rep("neutral",length(neutral_freqs)), rep("traits",length(trait_freqs)))
data <- cbind.data.frame(location, freqs, vars)
write.csv(data,"~/Dropbox/wagtails/data/displaced_cline_example.csv")

# wrangle cline data
cline_data <- lapply(replicates, `[[`, 2)
cline_data <- unlist(cline_data,recursive=FALSE)
cline_filter <- Filter(Negate(is.null), cline_data)
cline_df <- do.call("rbind", cline_filter)
cline_df$replicate <- rep(1:50, each=200)
write.csv(cline_df,"~/Dropbox/wagtails/data/simulation_cline_parameters.csv")

# extract final generation
cline_500 <- cline_df[cline_df$generation==500,]

# count number of times mating trait cline ended up in recessive pop
count <- 0
for(replicate in 1:max(cline_500$replicate)){
  tmp <- cline_500[cline_500$replicate==replicate,]
  if(tmp$cline_center[1]>tmp$cline_center[2]){
    count <- count + 1
  }
}
count/max(cline_500$replicate) # 64% of runs ended with mating trait cline moving left (towards recessive pop)

# count number of times mating trait cline ended up in recessive pop
count_2 <- 0
for(replicate in 1:max(cline_500$replicate)){
  tmp <- cline_500[cline_500$replicate==replicate,]
  if(tmp$cline_width[1]>tmp$cline_width[2]){
    count_2 <- count_2 + 1
  }
}
count_2/max(cline_500$replicate) # 98% of runs ended with mating trait cline more narrow

# mean cline center
mean_center <- ddply(cline_500, "cline_type", summarise, mean_center=mean(cline_center))

#         cline_type mean_center
# 1 Mating trait loci   0.4030753
# 2      Neutral loci   0.4760450

# cline center variance
sd_center <- ddply(cline_500, "cline_type", summarise, sd_center=sd(cline_center))

#     cline_type sd_center
# 1 Mating trait loci 0.1452830
# 2      Neutral loci 0.1360019

# mean cline width
mean_width <- ddply(cline_500[cline_500$cline_width<Inf,], "cline_type", summarise, mean_width=mean(cline_width, na.rm = T))

#         cline_type mean_width
# 1 Mating trait loci  0.09505139
# 2      Neutral loci  0.44965245

# cline width variance
sd_width <- ddply(cline_500[cline_500$cline_width<Inf,], "cline_type", summarise, sd_center=sd(cline_width))

# cline_type  sd_center
# 1 Mating trait loci 0.05421949
# 2      Neutral loci 0.14954128

# wilcox ranked sign test: for nonparametric, nonnormal, dependent distributions
neut <- cline_500[cline_500$cline_type=="Neutral loci",]$cline_center
pheno <- cline_500[cline_500$cline_type=="Mating trait loci",]$cline_center 
wilcox.test(neut, pheno, conf.int = T,paired=T) 

# data:  neut and pheno
# V = 932, p-value = 0.004539
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   0.02685437 0.12565345
# sample estimates:
#   (pseudo)median 
# 0.07860535 

neut_w <- cline_500[cline_500$cline_type=="Neutral loci",]$cline_width
pheno_w <- cline_500[cline_500$cline_type=="Mating trait loci",]$cline_width
wilcox.test(neut_w, pheno_w, conf.int = T, paired=T) # W = 1641, p-value = 0.007102

# Wilcoxon signed rank test with continuity correction
# 
# data:  neut_w and pheno_w
# V = 1272, p-value = 9.008e-10
# alternative hypothesis: true location shift is not equal to 0
# 0 percent confidence interval:
#   NaN NaN
# sample estimates:
#   midrange 
# Inf 


