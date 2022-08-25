cc_library(
  name = "maybe",
  srcs = [
    "maybe.h",
  ],
  hdrs = [
    "maybe.h"
  ],
)

cc_library(
  name = "mass",
  srcs = [
    "mass.h"
  ],
  hdrs = [
    "mass.h"
  ]
)

cc_library(
  name = "spring",
  srcs = [
    "spring.h",
    "spring.cc",
    "problem.h",
  ],
  hdrs = [
    "spring.h"
  ],
  deps = [
    ":mass",
    "@boost//:numeric_ublas"
  ]
)

cc_library(
  name = "problem",
  srcs = [
    "problem.h",
    "problem.cc"
  ],
  hdrs = [
    "problem.h",
  ],
  deps = [
    ":mass",
    ":spring",
    ":maybe",
    "@boost//:numeric_ublas",
    "@com_google_googletest//:gtest_main",
  ]
)

cc_test(
  name = "problem_test",
  size = "small",
  srcs = [
    "problem_test.cc",
    ],
  deps = [
    ":problem",
    "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
  name = "boost_test",
  size = "small",
  srcs = [
    "boost_test.cc",
    ],
  deps = [
    "@boost//:numeric_ublas",
    "@com_google_googletest//:gtest_main"
    ],
)

cc_test(
  name = "spring_test",
  size = "small",
  srcs = [
    "spring_test.cc"
    ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":problem",
    ],
)