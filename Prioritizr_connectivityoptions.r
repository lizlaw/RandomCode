## Tests out different connectivity options for Prioritizr

require(prioritizr)
require(gurobi)
require(Matrix)

# set seed for reproducibility
set.seed(500)
rescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

# load data
data(sim_pu_raster, sim_features, sim_locked_in_raster, sim_locked_out_raster)
pu <- sim_pu_raster^0
pu[] <- (1:100)/100
lockin <- pu*0
lockin[c(1,12,23,34,45,56,67,78,89,100)]<- 1

# boundary length matrix
#bm1 <- boundary_matrix(sim_pu_raster)
bm1 <- boundary_matrix(pu)
bm1[] <- rescale(bm1[])
image(bm1)

#cm1 <- connectivity_matrix(sim_pu_raster, !is.na(sim_pu_raster))
cm1 <- connectivity_matrix(pu, !is.na(pu))
cm1[] <- rescale(cm1[])

#r2 <- !is.na(sim_pu_raster)
#r2[sim_locked_in_raster==1] <- 0
#cm2 <- connectivity_matrix(sim_pu_raster,r2)
cm2 <- connectivity_matrix(pu,lockin)
cm2[] <- rescale(cm2[])

#r3 <- !is.na(sim_pu_raster)
#r3[sim_locked_out_raster==1] <- 0
#cm3 <- connectivity_matrix(sim_pu_raster,r3)
cm3 <- connectivity_matrix(pu,(1-lockin))
cm3[] <- rescale(cm3[])

r4 <- r2
r4[r3==0] <- 0
cm4 <- connectivity_matrix(sim_pu_raster,r4)
cm4[] <- rescale(cm4[])

centroids <- xyFromCell(pu,1:100)
d_matrix <- (1 / (as(dist(centroids), "Matrix") + 1))
d_matrix[] <- rescale(d_matrix[])
image(d_matrix)
dm7 <- d_matrix
dm7[d_matrix < 0.9] <- 0

plot(stack(!is.na(pu), r2,r3,r4))
image(bm1)
image(cm1)
image(cm2)
image(cm3)
image(cm4)

# create minimal problem
p0 <- problem(pu, sim_features) %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_binary_decisions() %>%
  add_default_solver()

# create problem with low boundary penalties
p1 <- p0 %>% add_connectivity_penalties(100, connectivity_data = bm1)
p2 <- p0 %>% add_connectivity_penalties(100, connectivity_data = cm2)
p3 <- p0 %>% add_connectivity_penalties(100, connectivity_data = cm3)
p4 <- p0 %>% add_connectivity_penalties(25, connectivity_data = bm1+cm2)
p5 <- p0 %>% add_connectivity_penalties(100, connectivity_data = bm1+cm3)
p6 <- p0 %>% add_connectivity_penalties(100, connectivity_data = dm7)


s0 <- solve(p0)
s1 <- solve(p1)
s2 <- solve(p2)
s3 <- solve(p3)
s4 <- solve(p4)
s5 <- solve(p5)
s6 <- solve(p6)

plot(stack(lockin, s0,s1,s2,s3, s4,s5, s6))
