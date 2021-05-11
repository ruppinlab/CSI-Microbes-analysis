
# this example works
vec <- c(.01, .049, .5, .0001, .5, .5, .5)
names(vec) <- c("family1", "genus1", "genus2", "species1", "species2", "species3", "species4")

parent <- c("family1", "family1", "genus1", "genus1", "genus2", "genus2")
child <- c("genus1", "genus2", "species1", "species2", "species3", "species4")
df <- data.frame(parent, child)
mat <- as.matrix(df)
colnames(mat) <- NULL

out <- hFDR.adjust(vec, mat, alpha = 0.05)

# we get an error for the below
vec <- c(.01, .049)
names(vec) <- c("family1", "genus1")

parent <- c("family1")
child <- c("genus1")
df <- data.frame(parent, child)
mat <- as.matrix(df)
colnames(mat) <- NULL

out <- hFDR.adjust(vec, mat, alpha = 0.05)


# we get an error for the below
vec <- c(.01, .049, .05)
names(vec) <- c("family1", "genus1", "species1")

parent <- c("family1", "genus1")
child <- c("genus1", "species1")
df <- data.frame(parent, child)
mat <- as.matrix(df)
colnames(mat) <- NULL

out <- hFDR.adjust(vec, mat, alpha = 0.05)

# we also get an error for this
vec <- c(.01, .001, .05, .001)
names(vec) <- c("family1", "genus1", "species1", "species2")

parent <- c("family1", "genus1", "genus1")
child <- c("genus1", "species1", "species2")
df <- data.frame(parent, child)
mat <- as.matrix(df)
colnames(mat) <- NULL

out <- hFDR.adjust(vec, mat, alpha = 0.05)

# no error for this
vec <- c(.01, .1, .05)
names(vec) <- c("family1", "genus1", "genus2")

parent <- c("family1", "family1")
child <- c("genus1", "genus2")
df <- data.frame(parent, child)
mat <- as.matrix(df)
colnames(mat) <- NULL

out <- hFDR.adjust(vec, mat, alpha = 0.05)
