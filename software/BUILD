cc_test(
  name = "problem_test",
  size = "small",
  srcs = [
    "problem_test.cc",
    "problem.cc",
    "problem.h",
    "mass.h",
    "maybe.h"
    ],
  deps = ["@com_google_googletest//:gtest_main"],
)

cc_test(
  name = "boost_test",
  size = "small",
  srcs = [
    "boost_test.cc",
    ],
  deps = [
    "@boost//:numeric_ublas",
    "@com_google_googletest//:gtest_main"],
)