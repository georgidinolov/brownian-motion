cc_binary(
	name = "2d-brownian-motion-data-sets",
	srcs = ["2d-brownian-motion-data-sets.cpp"],
	includes = ["2DBrownianMotionPath.hpp"],
	deps = [":2d-brownian-motion"],
)

cc_binary(
	name = "1d-brownian-motion-path-test",
	srcs = ["1d-brownian-motion-path.cpp"],
	deps = [":1d-brownian-motion"],
)

cc_binary(
	name = "2d-brownian-motion-path-test",
	srcs = ["2d-brownian-motion-path.cpp"],
	deps = [":2d-brownian-motion"],
)

cc_library(
	name = "1d-brownian-motion",
	srcs = ["1DBrownianMotionPath.cpp"],
	hdrs = ["1DBrownianMotionPath.hpp"],
	visibility = ["//visibility:public"],
	linkopts = ["-lm"],
)

cc_library(
	name = "2d-brownian-motion",
	srcs = ["2DBrownianMotionPath.cpp"],
	hdrs = ["2DBrownianMotionPath.hpp"],
	visibility = ["//visibility:public"],
	linkopts = ["-lm"],
)