
# DONT CHANGE
library(data.table)
library(kableExtra)
set.seed(123456)
n = 5
W1 <- rt(n, df = 5)
pi <- plogis(W1 )
A <- rbinom(n, 1, pi)
mu0 =  0.35 + 0.65*plogis( W1-2 )
mu1 <- mu0 + plogis(2*W1+2) - plogis( W1-2) - 0.349

mu0 =  0.35 + 0.65*plogis( W1-2 )
mu1 <- mu0 + plogis(2*W1+2) - plogis( W1-2) - 0.349
mu <- ifelse(A==1, mu1, mu0)
Y <- rbinom(n, 1, mu)
theta <- W1
betas <- seq(-2,2, length = 100)
risks <- sapply(betas, function(a) {
  theta <- a * theta
  loss <- (mu1 + mu0 +1/pi*(Y - mu))*log(1 + exp(theta)) - (mu1 + A/pi0*(Y - mu1))*theta
  mean(loss, na.rm = T)

})
weights <- round( (mu1 + mu0 + 1/pi*(Y - mu)),2)
outcome <- round((mu1  + A/pi0*(Y - mu1)) / weights,2)
dat <- data.frame(Covariate = round(W1,3), Weight  = weights, Outcome = outcome)
dat[1:5,]
kab <- kableExtra::kable(dat, booktabs = TRUE) %>% kableExtra::kable_styling(latex_options = c("striped", "hold_position", "scale_down"))

kab






# Intro DR-learner vs EP-learner example

library(delayed)
library(sl3)
library(hte3)
library(stringr)
set.seed(12345)
n = 1500
W1 <- rt(n, df = 5)
pi0 <- plogis(W1  )

A <- rbinom(n, 1, pi0)
#mu0 = plogis( W1-2)
#mu1 <- plogis(2*W1+2)
mu0 = plogis( W1-2)
mu1 <- plogis(2*W1+2)
mu0 =  0.35 + 0.65*plogis( W1-2 )
mu1 <- pmax(mu0 + plogis(2*W1+2) - plogis( W1-2) - 0.35, 0.001)
mu0 =  0.35 + 0.65*plogis( W1-2 )
mu1 <- mu0 + plogis(2*W1+2) - plogis( W1-2) - 0.349
mu <- ifelse(A==1, mu1, mu0)
Y <- rbinom(n, 1, mu)
trueCATE <- mu1 - mu0

lrnr_stack <- Stack$new(Lrnr_ranger$new(min.node.size=30,max.depth = 1),Lrnr_ranger$new(min.node.size=30,max.depth = 2),
                        Lrnr_ranger$new(min.node.size=30,max.depth = 3),Lrnr_ranger$new(min.node.size=30,max.depth = 4),
                        Lrnr_ranger$new(min.node.size=30,max.depth = 5),
                        Lrnr_ranger$new(min.node.size=30,max.depth = 6),
                        Lrnr_ranger$new(min.node.size=30,max.depth = 7),Lrnr_ranger$new(min.node.size=30,max.depth = 8))


data <- data.table(W1 = W1, A = A, Y = Y)
hte3_task <- hte3::make_hte3_Task_tx(data, modifiers = "W1", confounders = "W1",
                                     treatment = "A", outcome = "Y",
                                     learner_pi = lrnr_stack,
                                     learner_mu = lrnr_stack)

EY0W_est <- hte3_task$get_nuisance_estimates("mu")[,1]
EY1W_est <- hte3_task$get_nuisance_estimates("mu")[,2]
pA1W_est <- hte3_task$get_nuisance_estimates("pi")[,2]
pseudo_outcome <- EY1W_est - EY0W_est + (2*A - 1) * (1/ ifelse(A==1, pA1W_est, 1 - pA1W_est)) * (Y - ifelse(A==1, EY1W_est, EY0W_est ))

df <- data.frame(X = abs(pseudo_outcome))

df <- data.frame(X = abs(pseudo_outcome))
df <- df[order(df$X), , drop = FALSE]
df$ecdf <- 1 -  ecdf(df$X)(df$X)

ggplot(df, aes(x = X, y = ecdf)) +
  geom_step() +
  scale_x_continuous(limits = c(0, 5), breaks = seq(-6, 6, length = 13)) +
  ggtitle("Empirical RCDF of absolute value of pseudo-outcome values") +
  labs(x = "Pseudo-outcome absolute value", y = "Cumulative Probability") + theme_bw()
ggsave("HistogramDR.pdf", width = 6, height = 6)


ggplot(df, aes(x = X, y = ecdf)) +
  geom_step() + scale_y_continuous(limits = c(0, 0.05)) +   scale_x_continuous() +
  ggtitle("Empirical RCDF of absolute value of pseudo-outcome values") +
  labs(x = "Pseudo-outcome absolute value", y = "Cumulative Probability") + theme_bw()

ggsave("HistogramDRextreme.pdf", width = 6, height = 6)




##########
## GET DR-learner Predictions
##########
# Get DR-learner preds for CV and by tree depth
lrnr_DR <- cross_validate_cate(Lrnr_cate_DR$new(base_learner  = lrnr_stack), hte3_task)
lrnr_DR_cv <- lrnr_DR$lrnr_sl
lrnr_DR_by_depth <- lrnr_DR_cv$fit_object$full_fit$fit_object$learner_fits$Stack
# get predictions
DR_preds_by_depth <- lrnr_DR_by_depth$predict(hte3_task)
DR_preds_cv <- lrnr_DR_cv$predict(hte3_task)
# match predictions to tree depths
lrnr_names <- names(DR_preds_by_depth)
tree_depth <- stringr::str_match(lrnr_names, "_30_([0-9]+)")[,2]
all_depths <- unique(tree_depth)
all_depths_names <- paste0("Max Depth = ", all_depths)
# collect predictions
DR_preds <- data.table(
  preds = c(DR_preds_cv, unlist(DR_preds_by_depth)),
  depth = rep(c("CV-selected Max Depth", all_depths_names), each = nrow(DR_preds_by_depth)),
  Method = "DR-learner",
  Covariate = W1,
  trueCATE = trueCATE
)


##########
## GET T-learner Predictions
##########
# Get T-learner preds for CV and by tree depth
lrnr_T <- cross_validate_cate(Lrnr_cate_T$new(base_learner  = lrnr_stack), hte3_task)
lrnr_T_cv <- lrnr_T$lrnr_sl
lrnr_T_by_depth <- lrnr_T_cv$fit_object$full_fit$fit_object$learner_fits$Stack
# get predictions
T_preds_by_depth <- lrnr_T_by_depth$predict(hte3_task)
T_preds_cv <- lrnr_T_cv$predict(hte3_task)
# match predictions to tree depths
lrnr_names <- names(T_preds_by_depth)
tree_depth <- stringr::str_match(lrnr_names, "_30_([0-9]+)")[,2]
all_depths <- unique(tree_depth)
all_depths_names <- paste0("Max Depth = ", all_depths)
# collect predictions
T_preds <- data.table(
  preds = c(T_preds_cv, unlist(T_preds_by_depth)),
  depth = rep(c("CV-selected Max Depth", all_depths_names), each = nrow(T_preds_by_depth)),
  Method = "T-learner",
  Covariate = W1,
  trueCATE = trueCATE
)

##########
## GET EP-learner Predictions
##########
# Get eP-learner preds for CV and by tree depth
ep_stack <- Stack$new(
  Lrnr_cate_EP$new(base_learner  = lrnr_stack, sieve_num_basis = 1),
  Lrnr_cate_EP$new(base_learner  = lrnr_stack, sieve_num_basis = 2),
  Lrnr_cate_EP$new(base_learner  = lrnr_stack, sieve_num_basis = 3),
  Lrnr_cate_EP$new(base_learner  = lrnr_stack, sieve_num_basis = 4),
  Lrnr_cate_EP$new(base_learner  = lrnr_stack, sieve_num_basis = 5.)
)
# train learners
lrnr_EP <- cross_validate_cate(ep_stack, hte3_task)
lrnr_EP_cv <- lrnr_EP$lrnr_sl
lrnr_EP_stack <- lrnr_EP_cv$fit_object$full_fit$fit_object$learner_fits$Stack
# get predictions
ep_preds_all <- lrnr_EP_stack$predict(hte3_task)
ep_preds_cv <- lrnr_EP_cv$predict(hte3_task)
# extract best pred for each depth (CV over sieves)
lrnr_names <- names(ep_preds_all)
tree_depth <- stringr::str_match(lrnr_names, "_30_([0-9]+)")[,2]
sieve_num <- stringr::str_match(lrnr_names, "EP_([0-9]+)")[,2]
cv_risks <- lrnr_EP$cv_risk
all_depths <- unique(tree_depth)
best_index_by_depth <-c()
for(depth in all_depths) {
  best_index <- as.vector(which(tree_depth == depth & cv_risks ==  min(cv_risks[tree_depth == depth])  ))
  best_index_by_depth <- c(best_index_by_depth, best_index)
}
names(best_index_by_depth) <- all_depths
all_depths_names <- paste0("Max Depth = ", all_depths)
ep_preds_by_depth <- ep_preds_all[, best_index_by_depth, with = FALSE]
# collect predictions
EP_preds <- data.table(
  preds = c(ep_preds_cv, unlist(ep_preds_by_depth)),
  depth = rep(c("CV-selected Max Depth", all_depths_names), each = nrow(ep_preds_by_depth)),
  Method = "EP-learner",
  Covariate = W1,
  trueCATE = trueCATE
)

real_preds <- data.table(preds = trueCATE,
                         depth = rep(c("CV-selected Max Depth", all_depths_names), each = length(trueCATE)),
                         Method = "TRUE",
                         Covariate = W1,
                         trueCATE = trueCATE)


# get all preds and bound in [-1,1]
all_preds <- as.data.frame(rbindlist(list(EP_preds, DR_preds, T_preds , real_preds)))
all_preds$preds <- pmax(all_preds$preds, -1)
all_preds$preds <- pmin(all_preds$preds, 1)
all_preds_sub <- all_preds[  all_preds$depth %in%  paste0("Max Depth = ", c(1,2,4,7)),]




color_palette <- c("#0000FF", "#FF0000", "#339933", "black")
linetype_palette <-c("solid", "dashed", "dotdash", "dotted")
size_palette <- c(0.7,0.65, 0.65, 0.7)
names(color_palette) <- names(linetype_palette) <- names(size_palette) <- c("EP-learner", "DR-learner", "T-learner", "TRUE")



plt <- ggplot(all_preds_sub) +
  geom_line(aes(x = Covariate,
                y = preds, color = Method, linetype = Method ,size = Method) )  + scale_x_continuous(limits=c(-3,3)) + labs(x = "Covariate", y = "CATE")  + scale_y_continuous(limits = c(-1, 1) ) + theme_bw() + facet_grid(~ depth )+ scale_colour_manual(values = color_palette)+ scale_linetype_manual(values=linetype_palette)  +  scale_size_manual(  values=size_palette )  +
  theme(legend.position="none")


ggsave(file = "introDRlearner_bydepth.pdf", width = 9, height = 4)

all_preds_cv <- all_preds[  all_preds$depth %in% c("CV-selected Max Depth"),]

plt <- ggplot(all_preds_cv) +
  geom_line(aes(x = Covariate,
                y = preds, color = Method, linetype = Method ,size = Method) )  + scale_x_continuous(limits=c(-3,3)) + labs(x = "Covariate", y = "CATE")  + theme_bw() + facet_grid(~ depth )+ scale_colour_manual(values = color_palette)+ scale_linetype_manual(values=linetype_palette)  +  scale_size_manual(  values=size_palette )  + theme(legend.margin=margin(t=0, r=1, b=-1.5, l=1, unit="cm")) + scale_y_continuous(limits = c(-0.5, 0.5)) + facet_grid(~ depth )+
  theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(1.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'), #change legend key width
        legend.title =  element_blank(), #change legend title font size
        legend.text = element_text(size=13),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))   + geom_histogram(aes(x=Covariate))  +  geom_density(aes(x=Covariate, y = ..density..), color = NA, fill = "grey", alpha = 0.3)


ggsave(file = "introDRlearner_cv.pdf", width = 9, height =5)


# get MSE
all_preds <- as.data.table(all_preds)
mse <- all_preds[,   mean((preds - trueCATE)^2), by = c("depth", "Method")]
setkey(mse, depth)
mse
