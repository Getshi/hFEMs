#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <assert.h>
#include <stdarg.h>
#include <iostream>

#ifndef DEBUG_VERBOSE_LEVEL
// #ifdef NDEBUG // release mode
#define DEBUG_VERBOSE_LEVEL 3  // define from compiler
#endif

#define DECLARE_UNUSED(x) ((void)x);

namespace Debug {

// convenience type extensions
// std::ostream& operator<< (std::ostream& stream, const Magnum::Vector3& val);
// std::ostream& operator<< (std::ostream& stream, const Magnum::Vector2& val);

// printf wrapper
void logf(const char *format, ...);
void warningf(const char *format, ...);
void errorf(const char *format, ...);

// space-separated outstream
template <typename... Args>
void log(Args &&... args) {
#if DEBUG_VERBOSE_LEVEL >= 3
  std::cout << "Log: ";
  int dummy[] = {0, (std::cout << std::forward<Args>(args) << " ", 0)...};
  DECLARE_UNUSED(dummy)
  std::cout << std::endl;
#endif
}

template <typename... Args>
void warning(Args &&... args) {
#if DEBUG_VERBOSE_LEVEL >= 2
  std::cout << "Warning: ";
  int dummy[] = {0, (std::cout << std::forward<Args>(args) << " ", 0)...};
  DECLARE_UNUSED(dummy)
  std::cout << std::endl;
#endif
}

template <typename... Args>
void error(Args &&... args) {
#if DEBUG_VERBOSE_LEVEL >= 1
  std::cerr << "Error: ";
  int dummy[] = {0, (std::cerr << std::forward<Args>(args) << " ", 0)...};
  DECLARE_UNUSED(dummy)
  std::cerr << std::endl;
#endif
}

}  // namespace Debug

#endif /* end of include guard: _DEBUG_H_ */
