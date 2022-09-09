// Struct to represent functions returns that might have failed.
// Maybe.isError = true indicates that the function returned an error,
// and Maybe.errMsg has the corresponding error message.
// If the function has no return value, we should use Maybe<Void> so
// that we still can return isError and errMsg.

#pragma once
#include <string>

class Void {
   public:
    Void() {}
};

template <class c>
class Maybe {
   public:
    c val;
    bool isError;
    std::string errMsg;
    Maybe() {
        this->isError = false;
        this->errMsg = "";
        this->val = c();
    }
    Maybe<c>& operator=(const Maybe<c>& rhs) {
        this->isError = rhs.isError;
        this->errMsg = rhs.errMsg;
        this->val = rhs.val;
        return *this;
    }
};