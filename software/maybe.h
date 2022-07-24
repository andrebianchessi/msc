// Struct to represent functions returns that might have failed.
// Maybe.isError = true indicates that the function returned an error,
// and Maybe.errMsg has the corresponding error message.
// If the function has no return value, we should use Maybe<Void> so
// that we still can return isError and errMsg.

#pragma once
#include <string>

typedef void* Void;

template <class c>
struct Maybe {
    c* val;
    bool isError = false;
    std::string errMsg = "";
};