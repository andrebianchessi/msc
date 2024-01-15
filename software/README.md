## Requirements
- [Bazel](https://bazel.build/)
Last tested working version (`bazel version` output):
```
Build label: 7.0.0
Build target: @@//src/main/java/com/google/devtools/build/lib/bazel:BazelServer
Build time: Mon Dec 11 16:51:49 2023 (1702313509)
Build timestamp: 1702313509
Build timestamp as int: 1702313509
```

## Run tests
- `bazel test //:all`

To view the outputs:
`bazel test --test_output=all //:all`