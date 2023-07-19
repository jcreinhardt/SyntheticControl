# Setup
rm(list = ls())
gc()

# Libraries
library(dplyr)
library(Synth)
library(synthdid)

# Load data
data(basque)
df <- apply(basque, 2, function(x) if (is.numeric(x)) as.numeric(x) else x)

# --------------- Replicate Abadie et al. 2003 --------------------------

# Run baseline synthetic control
df_for_synth <- dataprep(foo = basque,
    predictors = c("school.illit",
        "school.prim",
        "school.med",
        "school.high",
        "school.post.high",
        "invest"),
    predictors.op = c("mean"),
    dependent = c("gdpcap"),
    unit.variable = c("regionno"),
    time.variable = c("year"),
    special.predictors = list(
        list("gdpcap", 1960:1969, c("mean")),
        list("sec.agriculture", seq(1961, 1969, 2), c("mean")),
        list("sec.energy", seq(1961, 1969, 2), c("mean")),
        list("sec.industry", seq(1961, 1969, 2), c("mean")),
        list("sec.construction", seq(1961, 1969, 2), c("mean")),
        list("sec.services.venta", seq(1961, 1969, 2), c("mean")),
        list("sec.services.nonventa", seq(1961, 1969, 2), c("mean")),
        list("popdens", 1969, c("mean"))),
    treatment.identifier = 17,
    controls.identifier = c(2:16, 18),
    time.predictors.prior = c(1964:1969),
    time.optimize.ssr = c(1960:1969),
    unit.names.variable = c("regionname"),
    time.plot = c(1955:1997))

# 1. combine highest and second highest
# schooling category and eliminate highest category
df_for_synth$X1["school.high", ] <- df_for_synth$X1["school.high", ] + df_for_synth$X1["school.post.high", ]
df_for_synth$X1 <- as.matrix(df_for_synth$X1[-which(rownames(df_for_synth$X1) == "school.post.high"), ])
df_for_synth$X0["school.high", ] <- df_for_synth$X0["school.high", ] + df_for_synth$X0["school.post.high", ]
df_for_synth$X0 <- df_for_synth$X0[-which(rownames(df_for_synth$X0) == "school.post.high"), ]

# 2. make total and compute shares for the schooling catgeories
lowest <- which(rownames(df_for_synth$X0) == "school.illit")
highest <- which(rownames(df_for_synth$X0) == "school.high")
df_for_synth$X1[lowest:highest, ] <- (100 * df_for_synth$X1[lowest:highest, ]) / sum(df_for_synth$X1[lowest:highest, ])
df_for_synth$X0[lowest:highest, ] <- 100 * scale(df_for_synth$X0[lowest:highest, ], center = FALSE, scale = colSums(df_for_synth$X0[lowest:highest, ]))

sc_results <- synth(data.prep.obj = df_for_synth)

print(synth.tab(
    dataprep.res = df_for_synth,
    synth.res = sc_results))

path.plot(synth.res = sc_results,
    dataprep.res = df_for_synth,
    Ylab = c("real per-capita GDP (1986 USD, thousand)"),
    Xlab = c("year"),
    Ylim = c(0, 13),
    Legend = c("Basque country", "synthetic Basque country"),
)

gaps.plot(synth.res = sc_results,
    dataprep.res = df_for_synth,
    Ylab = c("gap in real per-capita GDP (1986 USD, thousand)"),
    Xlab = c("year"),
    Ylim = c(-1.5, 1.5),
)

# --------------- Outcome-only matching --------------------------
dataprep_df <- dataprep(foo = basque,
    predictors = c("gdpcap"),
    predictors.op = c("mean"),
    dependent = c("gdpcap"),
    unit.variable = c("regionno"),
    time.variable = c("year"),
    treatment.identifier = 17,
    controls.identifier = c(2:16, 18),
    time.predictors.prior = c(1955:1969),
    time.optimize.ssr = c(1955:1969),
    unit.names.variable = c("regionname"),
    time.plot = c(1955:1997))

sc_est <- synth(dataprep_df)

path.plot(synth.res = sc_est,
    dataprep.res = dataprep_df,
    Ylab = c("real per-capita GDP (1986 USD, thousand)"),
    Xlab = c("year"),
    Ylim = c(0, 13),
    Legend = c("Basque country", "synthetic Basque country"),
)

weights <- sc_est$solution.w
y_hat <- matrix(filter(df, !regionno %in% c(1,17)) %>%
        select(gdpcap) %>%
        unlist,
    ncol = 43) %*% weights
y <- filter(df, regionno == 17) %>% select(gdpcap)
mean_squ_err <- (y - y_hat)^2
mse_ratio <- mean(mean_squ_err[16:43, ]) / mean(mean_squ_err[1:15, ])

# ---------- Compare MSEs from treated to controls ----------------
get_mse <- function(dataframe, treat_id, treat_time,
                    outcome, treat_id_true,
                    unit_var, time_var) {

    all_unit_ids <- unique(dataframe[, unit_var])
    control_ids <- setdiff(all_unit_ids, c(treat_id, treat_id_true))

    dataprep_out <- dataprep(foo = dataframe,
                             predictors = outcome,
                             predictors.op = "mean",
                             dependent = outcome,
                             unit.variable = unit_var,
                             time.variable = time_var,
                             treatment.identifier = treat_id,
                             controls.identifier = control_ids,
                             time.predictors.prior = min(dataframe[, time_var]):treat_time,
                             time.optimize.ssr = min(dataframe[, time_var]):treat_time,
                             time.plot = unique(dataframe[, time_var]))

    sc_estimate <- synth(dataprep_out)

    weights <- sc_estimate$solution.w
    y_hat <- matrix(data = dataframe[dataframe[,unit_var] %in% control_ids, outcome],
                    nrow = length(unique(dataframe[, time_var]))) %*% weights
    y <- dataframe[dataframe[, unit_var] == treat_id, outcome]
    mean_squ_err <- (y - y_hat)^2
    mean(mean_squ_err[16:43, ]) / mean(mean_squ_err[1:15, ])
}

sapply(setdiff(2:18, 17),
        get_mse,
        dataframe = filter(df, regionno != 1),
        treat_id_true = 17,
        treat_time = 1970,
        outcome = "gdpcap",
        unit_var = "regionno",
        time_var = "year") %>%
    hist(breaks = 10)



get_mse(dataframe = filter(df, regionno != 1),
        treat_id = 2,
        treat_id_true = 17,
        treat_time = 1970,
        outcome = "gdpcap",
        unit_var = "regionno",
        time_var = "year")



# ---------- Now do the same for California for increased sample size ----------------











