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
  ],
  hdrs = [
    "spring.h"
  ],
  deps = [
    ":mass",
    "@boost//:numeric_ublas",
    "@com_google_googletest//:gtest_main"
  ]
)
cc_test(
  name = "spring_test",
  size = "small",
  srcs = [
    "spring_test.cc"
    ],
  deps = [
    ":problem",
    ],
)

cc_library(
  name = "damper",
  srcs = [
    "damper.h",
    "damper.cc",
  ],
  hdrs = [
    "damper.h"
  ],
  deps = [
    ":mass",
    "@boost//:numeric_ublas",
    "@com_google_googletest//:gtest_main"
  ]
)
cc_test(
  name = "damper_test",
  size = "small",
  srcs = [
    "damper_test.cc"
    ],
  deps = [
    ":problem",
    ],
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
    ":damper",
    ":maybe",
    "@boost//:numeric_ublas",
    "@boost//:numeric_odeint",
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

cc_library(
  name = "bounded",
  srcs = [
    "bounded.cc"
  ],
  hdrs = [
    "bounded.h"
  ],
  deps = [
    ":maybe",
    "@com_google_googletest//:gtest_main"
  ]
)
cc_test(
  name = "bounded_test",
  size = "small",
  srcs = [
    "bounded_test.cc"
    ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":bounded",
    ],
)

cc_library(
  name = "utils",
  srcs = [
    "utils.h",
    "utils.cc"
  ],
  hdrs = [
    "utils.h"
  ],
  deps = [
    ":maybe",
    ":bounded"
  ]
)
cc_test(
  name = "utils_test",
  size = "small",
  srcs = [
    "utils_test.cc"
    ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":utils",
    ],
)

cc_library(
  name = "problem_description",
  srcs = [
    "problem_description.cc"
  ],
  hdrs = [
    "problem_description.h"
  ],
  deps = [
    ":problem",
    ":utils",
    ":bounded"
  ]
)
cc_test(
  name = "problem_description_test",
  size = "small",
  srcs = [
    "problem_description_test.cc"
    ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":problem_description",
    ],
)

cc_library(
  name = "creature",
  srcs = [
    "creature.h",
    "creature.cc"
  ],
  hdrs = [
    "creature.h"
  ],
  deps = [
    ":bounded",
    ":utils",
    "@com_google_googletest//:gtest_main",
  ]
)
cc_test(
  name = "creature_test",
  size = "small",
  srcs = [
    "creature_test.h",
    "creature_test.cc"
    ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":creature",
    ],
)

cc_library(
  name = "evolution",
  srcs = [
    "evolution.h",
    "evolution.tcc",
  ],
  hdrs = [
    "evolution.h",
  ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":creature",
    ":maybe",
  ]
)
cc_test(
  name = "evolution_test",
  size = "small",
  srcs = [
    "creature_test.h",
    "evolution_test.cc"
    ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":evolution",
    ],
)

cc_library(
  name = "problem_creature",
  srcs = [
    "problem_creature.h",
    "problem_creature.cc",
  ],
  hdrs = [
    "problem_creature.h",
  ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":creature",
    ":problem_description",
    ":bounded",
    ":utils",
    ":pimodel",
  ]
)
cc_test(
  name = "problem_creature_test",
  size = "small",
  srcs = [
    "problem_creature_test.cc"
    ],
  deps = [
    "@com_google_googletest//:gtest_main",
    ":problem_creature",
    ":evolution"
    ],
)

cc_library(
  name = "polynomial",
  srcs = [
    "polynomial.h",
    "polynomial.cc",
  ],
  hdrs = [
    "polynomial.h",
  ],
  deps = [
    "@boost//:math",
    "@boost//:lexical_cast",
    ":maybe",
    "@com_google_googletest//:gtest_main",
    ":utils"
  ]
)

cc_test(
  name = "polynomial_test",
  srcs = [
    "polynomial_test.cc",
  ],
  deps = [
    ":polynomial",
    "@com_google_googletest//:gtest_main",
    "@boost//:numeric_ublas",
    ":utils"
  ]
)

cc_library(
  name = "model",
  srcs = [
    "model.h",
    "model.cc",
  ],
  hdrs = [
    "model.h",
  ],
  deps = [
    ":maybe",
    ":utils",
    "@com_google_googletest//:gtest_main"
  ]
)

cc_test(
  name = "model_test",
  srcs = [
    "model_test.cc"
  ],
  deps = [
    ":model",
    "@com_google_googletest//:gtest_main",
    ":utils"
  ]
)

cc_library(
  name = "pimodel",
  srcs = [
    "pimodel.h",
    "pimodel.cc",
  ],
  hdrs = [
    "pimodel.h",
  ],
  deps = [
    ":bounded",
    ":maybe",
    ":model",
    ":polynomial",
    ":problem_description",
    "@com_google_googletest//:gtest_main",
    "@boost//:numeric_ublas",
  ]
)

cc_test(
  name = "pimodel_test",
  size = "enormous",
  srcs = [
    "pimodel_test.cc"
  ],
  deps = [
    ":pimodel",
    ":problem",
    "@com_google_googletest//:gtest_main",
    ":utils"
  ]
)

cc_binary(
  name = "pimodel_binary_test",
  srcs = [
    "pimodel_binary_test.cc"
  ],
  deps = [
    ":pimodel",
    ":utils"
  ]
)

cc_binary(
  name = "experiment_1",
  srcs = [
    "experiment_1.cc"
  ],
  deps = [
    ":bounded",
    ":evolution",
    ":utils",
    ":problem_creature"
  ]
)
