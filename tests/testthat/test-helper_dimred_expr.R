context("Dimred for expression")

test_that("Retrieving dimred_methods", {
  methods <- dynmethods:::list_dimred_methods()
  expect_is(methods, "list")
  expect_gt(length(methods), 0)
  expect_named(methods)
})

tasks <- dyntoy::toy_tasks %>% sample_n(3)

for (taski in seq_len(nrow(tasks))) {
  task <- extract_row_to_list(tasks, taski)
  expr <- log2(task$counts+1)

  expect_error(dynmethods:::dimred(expr, method = "prism_in_the_steets_of_london", ndim = 5))

  methods <- dynmethods:::list_dimred_methods()

  for (method_name in names(methods)) {
    test_that(paste0("Perform dimred ", method_name, " with expression on task ", task$id), {
      method_fun <- methods[[method_name]]
      expect_is(method_fun, "function")

      pdf("/dev/null")
      sink("/dev/null")
      out1 <- method_fun(expr)
      sink()

      expect_is(out1, "matrix")
      expect_identical(rownames(out1), rownames(expr))
      expect_identical(colnames(out1), paste0("Comp", seq_len(ncol(out1))))

      sink("/dev/null")
      out2 <- method_fun(expr, ndim = 2)
      sink()
      expect_is(out2, "matrix")
      expect_identical(rownames(out2), rownames(expr))
      expect_identical(colnames(out2), paste0("Comp", seq_len(ncol(out2))))
      expect_equal(ncol(out2), 2)

      sink("/dev/null")
      out3 <- dynmethods:::dimred(expr, method = method_name, ndim = 5)
      sink()
      expect_is(out3, "matrix")
      expect_identical(rownames(out3), rownames(expr))
      expect_identical(colnames(out3), paste0("Comp", seq_len(ncol(out3))))
      expect_equal(ncol(out3), 5)

      v1 <- as.vector(as.matrix(dist(out1)))
      v2 <- as.vector(as.matrix(dist(out2)))
      v3 <- as.vector(as.matrix(dist(out3)))

      # even though the number of dimensions is different, one would expect the correlation between each execution to be somewhat greater than zero
      expect_true(all(cor(cbind(v1, v2, v3)) > .1))
      dev.off()
    })
  }
}


test_that("Testing process_dimred", {
  sp <- matrix(runif(1:10), nrow = 2)
  new_sp <- dynmethods:::process_dimred(sp, c("A", "B"))
  expect_identical(rownames(new_sp), c("A", "B"))
  expect_identical(colnames(new_sp), paste0("Comp", seq_len(5)))
})
