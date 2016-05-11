##############
context("bsf")
##############

data <- cbind(v1 = c(1,1,1,2,2,3,3,3,3,3),
              v2 = c(1,1,3,3,3,3,3,3,3,3),
              v3 = c(2,2,1,1,1,1,1,3,3,3))
rownames(data) <- paste0("O", 1:10)

test_that("bsf basic properties", {
    expect_that(
        bsf.dist(rep(1:10,each=3), rep(1:10,each=3), 10),
        equals(0)
    )
    expect_that(
        bsf.dist(rep(1:10,each=3), rep(10:1,each=3), 10),
        equals(1)
    )
})

test_that("bsf results without normalization", {
    expect_that(
        bsf.dist(data[,1], data[,2], 3, NULL, 1),
        equals(0.45)
    )
    expect_that(
        bsf.dist(data[,1], data[,3], 3, NULL, 1),
        equals(47/30)
    )
    expect_that(
        bsf.dist(data[,2], data[,3], 3, NULL, 1),
        equals(2.25)
    )
})

test_that("bsf results with normalization", {
    maxD <- sum(abs(1:3-(3:1)))
    expect_that(
        bsf.dist(data[,1], data[,2], 3),
        equals(0.45/maxD)
    )
    expect_that(
        bsf.dist(data[,1], data[,3], 3),
        equals(47/30/maxD)
    )
    expect_that(
        bsf.dist(data[,2], data[,3], 3),
        equals(2.25/maxD)
    )
})

test_that("bsf results with precomputed SF", {
    SF <- abs(row(matrix(NA,3,3))-col(matrix(NA,3,3)))
    expect_that(
        bsf.dist(data[,1], data[,2], 3, SF, 1),
        equals(0.45)
    )
    expect_that(
        bsf.dist(data[,1], data[,3], 3, SF, 1),
        equals(47/30)
    )
    expect_that(
        bsf.dist(data[,2], data[,3], 3, SF, 1),
        equals(2.25)
    )
})


##############
context("qep")
##############

qp <- qep(data)

test_that("summary output", {

    out <- paste("Number of conditions: 3 ",
                 "Number of observations: 10 ",
                 "",
                 "Unique types of binning:",
                 "   conditions",
                 "    1 1 1",
                 "  1 3 2 5",
                 "  2 2 0 2",
                 "  3 5 8 3", sep="\n")

    expect_output(summary(qp), out)
})

test_that("qep format", {
    expect_match(class(qp), "qep")
    expect_equal(dim(qp), c(10,3))
    expect_equal(nrow(qp), 10)
    expect_equal(ncol(qp), 3)
    expect_equal(names(dimnames(qp))[1], "observations")
    expect_equal(names(dimnames(qp))[2], "conditions")
    expect_equal(rownames(qp), paste0("O", 1:10))
    expect_equal(colnames(qp), paste0("v", 1:3))
})

test_that("bsf dist on qep", {
    d <- dist(qp)
    dm <- as.matrix(d)
    maxD <- sum(abs(1:3-(3:1)))
    expect_match(class(d), "dist")
    expect_equal(length(d), 3)
    expect_equal(dm[1,2], 0.45/maxD)
    expect_equal(dm[1,3], 47/30/maxD)
    expect_equal(dm[2,3], 2.25/maxD)
})


##############
context("qmix")
##############

data2 <- cbind(v2.1 = c(1,1,1,2,2,3,3),
               v2.2 = c(1,1,3,3,3,3,3))
rownames(data2) <- c(paste0("O",3:5), paste0("O",7:10))
data3 <- cbind(v3.1 = c(1,1,2,2,2,3,3,3,1,3,1),
               v3.2 = c(1,1,3,2,3,3,3,3,1,3,2),
               v3.3 = c(2,2,1,1,2,1,1,3,1,3,3))
rownames(data3) <- c(paste0("O",4:6), paste0("O",7:14))

qm <- qepmix(list(qep(data), qep(data2), qep(data3)))

v3.m <- data[intersect(rownames(data),rownames(data2)),3]
v4.m <- data2[intersect(rownames(data),rownames(data2)),1]

test_that("bsf dist on qepmix", {    
    expect_equal(bsf.dist(v3.m,v4.m,3,maxD=1), 9/20)
})
