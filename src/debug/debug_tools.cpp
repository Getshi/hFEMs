#include "debug_tools.h"

// convenience type extensions
// std::ostream& Debug::operator<< (std::ostream& stream, const Magnum::Vector3&
// val) {
//   stream << "("<< val[0] << " " << val[1] << " " << val[2] << ")";
//   return stream;
// }
// std::ostream& Debug::operator<< (std::ostream& stream, const Magnum::Vector2&
// val) {
//   stream << "("<< val[0] << " " << val[1] << ")";
//   return stream;
// }

void Debug::logf(const char *format, ...) {
#if DEBUG_VERBOSE_LEVEL >= 3
  std::cout << "Log: ";
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
#endif
}

void Debug::warningf(const char *format, ...) {
#if DEBUG_VERBOSE_LEVEL >= 2
  std::cout << "Warning: ";
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
#endif
}

void Debug::errorf(const char *format, ...) {
#if DEBUG_VERBOSE_LEVEL >= 1
  std::cerr << "Error: ";
  va_list args;
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
#endif
}
